// Long-time integration with affine sets and the mean-value form, built on Kv.
//
// This is the algorithm of kv::ode_maffine / kv::odelong_maffine (Kv 0.4.62) with a different
// step-size strategy.  Kv's own driver lets the hull integration (step 3 below) choose the step;
// Kv's step-size heuristic then reacts to the width of the hull of the set and makes very small
// steps once the set has been stretched by the dynamics.  Here the step is proposed by the
// integration of the centre point and then (adaptively) shortened until the hull integration
// succeeds.  Each step:
//
//   1. I = interval hull of the affine set x,  c = mid(I) (a point, as thin intervals);
//   2. h0 = step chosen by kv::ode (automatic step size) for the centre point; propose
//      t1 = t + kappa h0;
//   3. D  = enclosure of D_x phi(t1; t, I) over the hull  (kv::ode on autodif, fixed step);
//      if this fails, t1 is moved closer to t (kappa is halved) and 3. is redone;
//   4. fc = enclosure of phi(t1; t, c)                    (kv::ode, interval, fixed step);
//   5. x <- fc + D (x - c)   in affine arithmetic.
//
// Correctness (mean value theorem, row by row): for every point x of the set,
// phi_i(x) = phi_i(c) + D_x phi_i(xi_i) (x - c) with xi_i in I, and D_x phi_i(xi_i) lies in row i
// of D.  Affine arithmetic keeps the dependence on the noise symbols of x, so the joint
// relation between initial and final points is preserved.  After the step the noise symbols
// created during the step are collected into one new symbol per component
// (kv::epsilon_reduce2, as in Kv).  The step size only affects efficiency, never correctness:
// every step is a verified Kv enclosure (kv::ode returns 2 only if the Taylor remainder was
// verified on the a-priori enclosure for the whole step).
//
// The rigorous one-step enclosures (Taylor method with verified remainder, kv::ode) are Kv's.

#ifndef MV_DRIVER_HPP
#define MV_DRIVER_HPP

#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <kv/ode-maffine.hpp>

namespace mv {

namespace ub = boost::numeric::ublas;
typedef kv::interval<double> itv;
typedef kv::affine<double> aff;

struct Stats {
	long steps = 0, shrinks = 0;
};

template <class F, class V>
int try_ode(F& f, V& x, const itv& t, itv& t1, const kv::ode_param<double>& p) {
	try {
		return kv::ode(f, x, t, t1, p);
	} catch (std::domain_error&) {
		return 0;
	}
}

struct NoPre {
	void operator()() const {}
};

// integrate the affine set x from t = start to t = end (thin or interval end time).
// returns true on success (end reached); false if a step could not be verified.
// pre() is called before every step (used for the exact power-of-two rescaling of the jet
// components, which may also change parameters of f; f is therefore taken by reference).
template <class F, class Pre = NoPre>
bool integrate(F& f, ub::vector<aff>& x, const itv& start, const itv& end, int order, Stats& st,
               int verbose = 0, Pre pre = Pre()) {
	const int n = x.size();
	kv::ode_param<double> pa, pf;
	pa.set_order(order).set_autostep(true).set_restart_max(10);
	pf.set_order(order).set_autostep(false);
	itv t = start;
	double kappa = 0.5;
	while (true) {
		pre();
		ub::vector<itv> I(n), c(n);
		for (int i = 0; i < n; i++) {
			I(i) = to_interval(x(i));
			c(i) = mid(I(i));
		}
		// 2. step proposed by the centre point
		itv t1 = end;
		{
			ub::vector<itv> tmp = c;
			int r = try_ode(f, tmp, t, t1, pa);
			if (r == 0) {
				if (verbose) std::cout << "mv: point integration failed at t = " << t << std::endl;
				return false;
			}
		}
		double h0 = mid(t1) - t.upper();
		bool last;
		ub::matrix<itv> D;
		ub::vector<itv> fc;
		int tries = 0;
		while (true) {
			double tn = t.upper() + kappa * h0;
			if (tn >= end.lower()) { t1 = end; last = true; }
			else { t1 = itv(tn); last = false; }
			if (!(t1.lower() > t.upper())) return false;
			// 3. Jacobian over the hull, fixed step
			ub::vector<kv::autodif<itv> > Iad = kv::autodif<itv>::init(I);
			itv t1b = t1;
			int r2 = try_ode(f, Iad, t, t1b, pf);
			if (r2 == 2) {
				// 4. centre point, same fixed step
				fc = c;
				t1b = t1;
				int r3 = try_ode(f, fc, t, t1b, pf);
				if (r3 == 2) {
					ub::vector<itv> vi;
					kv::autodif<itv>::split(Iad, vi, D);
					break;
				}
			}
			st.shrinks++;
			kappa *= 0.5;
			if (++tries > 20) {
				if (verbose) std::cout << "mv: step failed at t = " << t << std::endl;
				return false;
			}
		}
		if (tries == 0) kappa = std::min(1.0, kappa * 1.25);
		// 5. mean value form in affine arithmetic
		int maxnum_save = aff::maxnum();
		ub::vector<aff> xc(n), res(n);
		for (int i = 0; i < n; i++) xc(i) = x(i) - c(i);
		for (int i = 0; i < n; i++) {
			res(i) = aff(fc(i));
			for (int j = 0; j < n; j++) res(i) += D(i, j) * xc(j);
		}
		kv::epsilon_reduce2(res, maxnum_save);
		x = res;
		st.steps++;
		if (verbose >= 2 && st.steps % verbose == 0)
			std::cout << "mv: step " << st.steps << " t = " << mid(t1) << " h " << mid(t1) - t.upper()
			          << " kappa " << kappa << " maxnum " << aff::maxnum() << std::endl;
		if (last) return true;
		t = t1;
	}
}

} // namespace mv

#endif
