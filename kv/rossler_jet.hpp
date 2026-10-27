// Rossler system  x' = -(y+z),  y' = x + a y,  z' = b + z (x - c)
// in polar coordinates  x = r cos(theta), y = r sin(theta), reparametrised by the polar angle.
//
//   theta' = (x^2 + a x y + y^2 + y z)/r^2 = 1 + a C S + S z / r,      C = cos(theta), S = sin(theta)
//   r'     = (x x' + y y')/r              = a r S^2 - C z
//   z'     = b + z (r C - c)
//   dr/dtheta = r (a r S^2 - C z) / D,   dz/dtheta = r (b + z (r C - c)) / D,   D = r (1 + a C S) + S z  (= r theta').
//
// Independent variable t = theta + pi/2, so C = cos(theta) = sin t and S = sin(theta) = -cos t.
// The section Sigma = {x = 0, x' > 0} is {t = 0 mod 2 pi} as long as theta' > 0 (i.e. D > 0)
// along the orbit: at t = 0 mod 2pi we have x = 0, y = -r and x' = r theta' > 0, while the other
// zeros of x (t = pi mod 2pi) have x' = -r theta' < 0.  Hence the Poincare map P_c in the section
// coordinates (y,z) = (-r,z) is the time-2pi map of this 2pi-periodic system, and P_c^N is the
// time-2piN map.  D > 0 is guaranteed rigorously: Kv evaluates the field (including the division
// by D) on the a-priori enclosure of every Taylor step; an interval division by an interval
// containing 0 throws, which makes the step (and in the end the proof) fail.  D = r - z > 0 at
// t = 0 (checked by the caller) and D != 0 on every step give D > 0 by continuity.
//
// Jet.  The state also carries derivatives of (r,z) with respect to NG parameters p_k of the
// initial condition: the initial point is X + sum_k p_k d_k ((y,z) coordinates, d_k = "seed"
// vectors chosen by the caller), and in slot CSLOT (if CSLOT >= 0) also c = c0 + cseed p_CSLOT
// (cseed = 1 unless the driver has rescaled the parameters by powers of two, see below).
// The seeds are entered as initial values of the derivative components.  NH = 0: first
// derivatives only.  NH = 1: also the second derivative w.r.t. (p_0,p_0).  NH = 3: second
// derivatives w.r.t. (p_0,p_0), (p_0,p_1), (p_1,p_1).
// The right-hand side of the jet is obtained by evaluating the field in forward-mode jet
// arithmetic (exact chain rule), i.e. it is exactly the first / second variational equation.
//
// State layout (dimension n = 3 + 2 NG + 2 NH):
//   0: r   1: z   2: c
//   3 .. 3+NG-1           : r_k        3+NG .. 3+2NG-1         : z_k
//   3+2NG .. 3+2NG+NH-1   : r_(ij)     3+2NG+NH .. 3+2NG+2NH-1 : z_(ij)
//
// Written from the mathematics for this independent check; no code from the CAPD program.

#ifndef ROSSLER_JET_HPP
#define ROSSLER_JET_HPP

#include <boost/numeric/ublas/vector.hpp>
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>

namespace rj {

namespace ub = boost::numeric::ublas;
typedef kv::interval<double> itv;

// index pairs of the stored second derivatives: h[0]=(0,0), h[1]=(0,1), h[2]=(1,1)
static const int HI[3] = {0, 0, 1};
static const int HJ[3] = {0, 1, 1};

// forward-mode jet: value, NG first derivatives, NH second derivatives
template <class T, int NG, int NH> struct Jet {
	T v;
	T g[NG > 0 ? NG : 1];
	T h[NH > 0 ? NH : 1];
};

template <class T, int NG, int NH> inline Jet<T, NG, NH> jadd(const Jet<T, NG, NH>& x, const Jet<T, NG, NH>& y) {
	Jet<T, NG, NH> r;
	r.v = x.v + y.v;
	for (int k = 0; k < NG; k++) r.g[k] = x.g[k] + y.g[k];
	for (int k = 0; k < NH; k++) r.h[k] = x.h[k] + y.h[k];
	return r;
}
template <class T, int NG, int NH> inline Jet<T, NG, NH> jsub(const Jet<T, NG, NH>& x, const Jet<T, NG, NH>& y) {
	Jet<T, NG, NH> r;
	r.v = x.v - y.v;
	for (int k = 0; k < NG; k++) r.g[k] = x.g[k] - y.g[k];
	for (int k = 0; k < NH; k++) r.h[k] = x.h[k] - y.h[k];
	return r;
}
// multiplication by a quantity that does not depend on the parameters p_k
template <class T, int NG, int NH, class S> inline Jet<T, NG, NH> jscal(const Jet<T, NG, NH>& x, const S& s) {
	Jet<T, NG, NH> r;
	r.v = x.v * s;
	for (int k = 0; k < NG; k++) r.g[k] = x.g[k] * s;
	for (int k = 0; k < NH; k++) r.h[k] = x.h[k] * s;
	return r;
}
// Scaled storage.  The driver keeps the components of moderate size by exact power-of-two
// rescalings: g[k] = 2^{-E_k} (d/dp_k), h[m] = 2^{-F_m} (d^2/dp_i dp_j), (i,j) = (HI[m],HJ[m]).
// Then the product and quotient rules below hold with the factor
// lam[m] = 2^{E_i + E_j - F_m} in front of the terms that are products of two first
// derivatives (lam = 1 for the unscaled jet).  Linear operations are unaffected.
//
// product rule: (xy)_k = x_k y + x y_k,  (xy)_ij = x_ij y + x y_ij + lam (x_i y_j + x_j y_i)
template <class T, int NG, int NH> inline Jet<T, NG, NH> jmul(const Jet<T, NG, NH>& x, const Jet<T, NG, NH>& y,
                                                              const double* lam) {
	Jet<T, NG, NH> r;
	r.v = x.v * y.v;
	for (int k = 0; k < NG; k++) r.g[k] = x.g[k] * y.v + x.v * y.g[k];
	for (int k = 0; k < NH; k++) {
		int i = HI[k], j = HJ[k];
		r.h[k] = x.h[k] * y.v + x.v * y.h[k] + (x.g[i] * y.g[j] + x.g[j] * y.g[i]) * lam[k];
	}
	return r;
}
// q = x / y :  x = q y  =>  q_k = (x_k - q y_k)/y,  q_ij = (x_ij - q y_ij - lam (q_i y_j + q_j y_i))/y
template <class T, int NG, int NH> inline Jet<T, NG, NH> jdiv(const Jet<T, NG, NH>& x, const Jet<T, NG, NH>& y,
                                                              const double* lam) {
	Jet<T, NG, NH> r;
	r.v = x.v / y.v;
	for (int k = 0; k < NG; k++) r.g[k] = (x.g[k] - r.v * y.g[k]) / y.v;
	for (int k = 0; k < NH; k++) {
		int i = HI[k], j = HJ[k];
		r.h[k] = (x.h[k] - r.v * y.h[k] - (r.g[i] * y.g[j] + r.g[j] * y.g[i]) * lam[k]) / y.v;
	}
	return r;
}

// The field.  a, b are intervals containing 0.2; c is a state variable (c' = 0).
template <int NG, int NH, int CSLOT> struct RosslerJet {
	static const int N = 3 + 2 * NG + 2 * NH;
	static const int IG = 3;              // first index of r_k
	static const int IH = 3 + 2 * NG;     // first index of r_(ij)
	itv a, b;
	// dc/dp_CSLOT in the scaled storage (= 2^{-E_CSLOT}); 1 when unscaled.  The driver
	// rescales the stored jet components by powers of two to keep them of moderate size
	// (Kv's step-size control makes tiny steps when some components are large).
	double cseed = 1.0;
	// factors lam[m] of the scaled second-order product rule (see jmul); 1 when unscaled.
	double hlam[3] = {1.0, 1.0, 1.0};
	RosslerJet(const itv& a_, const itv& b_) : a(a_), b(b_) {
		static_assert(NH == 0 || NH == 1 || NH == 3, "NH must be 0, 1 or 3");
		static_assert(NH == 0 || NG >= (NH == 3 ? 2 : 1), "second derivatives need the slots");
		static_assert(CSLOT < NG, "CSLOT out of range");
	}

	template <class T> ub::vector<T> operator()(const ub::vector<T>& x, const T& t) {
		ub::vector<T> y(N);
		T C = sin(t);       // cos(theta)
		T S = -cos(t);      // sin(theta)
		T zero = x(2) * 0.;
		Jet<T, NG, NH> r, z;
		r.v = x(0); z.v = x(1);
		for (int k = 0; k < NG; k++) { r.g[k] = x(IG + k); z.g[k] = x(IG + NG + k); }
		for (int k = 0; k < NH; k++) { r.h[k] = x(IH + k); z.h[k] = x(IH + NH + k); }

		T K1 = a * C * S + 1.;    // 1 + a C S
		T K2 = a * S * S;         // a S^2
		Jet<T, NG, NH> D = jadd(jscal(r, K1), jscal(z, S));        // r (1 + a C S) + S z
		Jet<T, NG, NH> A = jsub(jscal(r, K2), jscal(z, C));        // a r S^2 - C z
		Jet<T, NG, NH> B = jscal(r, C);                            // r C - c   (c_k = [k == CSLOT])
		B.v = B.v - x(2);
		if (CSLOT >= 0) B.g[CSLOT] = B.g[CSLOT] - cseed;
		Jet<T, NG, NH> Bz = jmul(z, B, hlam);                      // z (r C - c)
		Bz.v = Bz.v + b;                                           // b + z (r C - c)
		Jet<T, NG, NH> fr = jdiv(jmul(r, A, hlam), D, hlam);
		Jet<T, NG, NH> fz = jdiv(jmul(r, Bz, hlam), D, hlam);

		y(0) = fr.v;
		y(1) = fz.v;
		y(2) = zero;
		for (int k = 0; k < NG; k++) { y(IG + k) = fr.g[k]; y(IG + NG + k) = fz.g[k]; }
		for (int k = 0; k < NH; k++) { y(IH + k) = fr.h[k]; y(IH + NH + k) = fz.h[k]; }
		return y;
	}
};

} // namespace rj

#endif
