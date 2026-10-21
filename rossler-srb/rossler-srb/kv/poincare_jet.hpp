// N-fold Poincare map with jets: one revolution (t from 0 to 2 pi) at a time, with exact
// power-of-two rescaling of the jet components before every integration step.
//
// The field is 2 pi-periodic in t, so revolution k can be integrated on [0, 2pi] (the end time
// is the interval enclosure of 2 pi; the enclosure at that end time contains the state at the
// exact time 2 pi).
//
// Rescaling.  The jet components of slot k are stored as  2^{-E[k]} * (true derivative w.r.t.
// p_k), second-order component m as 2^{-F[m]} * (true second derivative).  The field then uses
// dc/dp_CSLOT -> 2^{-E[CSLOT]} (cseed) and the factors hlam[m] = 2^{E[i]+E[j]-F[m]} in the
// second-order product rule (see rossler_jet.hpp), which makes the stored jet satisfy exactly
// the rescaled variational equations.  Multiplication of affine forms by powers of two is exact.

#ifndef POINCARE_JET_HPP
#define POINCARE_JET_HPP

#include <cmath>
#include <climits>
#include <vector>
#include <chrono>
#include <kv/constants.hpp>
#include "rossler_jet.hpp"
#include "mv_driver.hpp"

namespace pj {

namespace ub = boost::numeric::ublas;
typedef kv::interval<double> itv;
typedef kv::affine<double> aff;

inline double magI(const itv& x) { return std::max(std::fabs(x.lower()), std::fabs(x.upper())); }

template <int NG, int NH, int CSLOT> struct JetState {
	typedef rj::RosslerJet<NG, NH, CSLOT> Field;
	ub::vector<aff> x;       // the affine state (scaled jet)
	int E[NG > 0 ? NG : 1];  // scale exponents of the first-order slots
	int F[NH > 0 ? NH : 1];  // scale exponents of the second-order components
	JetState() : x(Field::N) {
		for (int k = 0; k < NG; k++) E[k] = 0;
		for (int m = 0; m < NH; m++) F[m] = 0;
	}

	// true (unscaled) components as affine forms (multiplication by powers of two is exact)
	aff r() const { return x(0); }
	aff z() const { return x(1); }
	aff c() const { return x(2); }
	aff rg(int k) const { return x(Field::IG + k) * std::ldexp(1.0, E[k]); }
	aff zg(int k) const { return x(Field::IG + NG + k) * std::ldexp(1.0, E[k]); }
	aff rh(int m) const { return x(Field::IH + m) * std::ldexp(1.0, F[m]); }
	aff zh(int m) const { return x(Field::IH + NH + m) * std::ldexp(1.0, F[m]); }

	// exponent of the stored magnitude (INT_MIN if the components are exactly zero)
	static int expo(const aff& u, const aff& v) {
		double m = std::max(magI(to_interval(u)), magI(to_interval(v)));
		return (m > 0 && std::isfinite(m)) ? std::ilogb(m) : INT_MIN;
	}
	// Rescaling rule (exact powers of two).  Kv's step-size control makes small steps when some
	// component has large Taylor coefficients, so components are only ever scaled DOWN:
	//   E_k = max(0, exponent of |true slot k|),
	//   F_m = max(E_i + E_j, exponent of |true second derivative m|),
	// which keeps the stored components below 2 in magnitude and the factors cseed and
	// hlam (multiplying the forcing terms of the variational equations) at most 1.
	void rescale() {
		int En[NG > 0 ? NG : 1];
		for (int k = 0; k < NG; k++) {
			int d = expo(x(Field::IG + k), x(Field::IG + NG + k));
			En[k] = (d == INT_MIN) ? 0 : std::max(0, d + E[k]);
			double f = std::ldexp(1.0, E[k] - En[k]);
			x(Field::IG + k) = x(Field::IG + k) * f;
			x(Field::IG + NG + k) = x(Field::IG + NG + k) * f;
			E[k] = En[k];
		}
		for (int m = 0; m < NH; m++) {
			int d = expo(x(Field::IH + m), x(Field::IH + NH + m));
			int base = E[rj::HI[m]] + E[rj::HJ[m]];
			int Fn = (d == INT_MIN) ? base : std::max(base, d + F[m]);
			double f = std::ldexp(1.0, F[m] - Fn);
			x(Field::IH + m) = x(Field::IH + m) * f;
			x(Field::IH + NH + m) = x(Field::IH + NH + m) * f;
			F[m] = Fn;
		}
	}
	// field parameters matching the current scaling
	void setField(Field& f) const {
		f.cseed = (CSLOT >= 0) ? std::ldexp(1.0, -E[CSLOT >= 0 ? CSLOT : 0]) : 1.0;
		for (int m = 0; m < NH; m++) f.hlam[m] = std::ldexp(1.0, E[rj::HI[m]] + E[rj::HJ[m]] - F[m]);
	}
};

struct RevInfo {
	ub::vector<itv> hull;   // hull of (r, z, c) after the revolution
	long steps;
	double seconds;
};

// N revolutions. Returns false if the integration failed.
template <int NG, int NH, int CSLOT>
bool iterate(const itv& a, const itv& b, JetState<NG, NH, CSLOT>& s, int N, int order,
             mv::Stats& st, std::vector<RevInfo>* info = nullptr, int verbose = 0) {
	typedef rj::RosslerJet<NG, NH, CSLOT> Field;
	Field f(a, b);
	const itv twopi = kv::constants<itv>::pi() * 2.;
	// rescaling before every step (exact; see above)
	auto pre = [&]() {
		s.rescale();
		s.setField(f);
	};
	for (int k = 0; k < N; k++) {
		long st0 = st.steps;
		auto t0 = std::chrono::steady_clock::now();
		if (!mv::integrate(f, s.x, itv(0.), twopi, order, st, verbose >= 2 ? verbose : 0, pre)) return false;
		double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
		if (info) {
			RevInfo ri;
			ri.hull.resize(3);
			for (int i = 0; i < 3; i++) ri.hull(i) = to_interval(s.x(i));
			ri.steps = st.steps - st0;
			ri.seconds = sec;
			info->push_back(ri);
		}
		if (verbose) std::cout << "  revolution " << k + 1 << ": steps " << st.steps - st0 << ", " << sec
		                       << " s, maxnum " << aff::maxnum() << std::endl;
	}
	s.rescale();
	return true;
}

} // namespace pj

#endif
