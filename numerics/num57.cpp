// Non-rigorous (double precision) numerics for homoclinic tangencies of the period-1 saddle p_c
// of the Rossler return map near c = 5.7 (a = b = 0.2).  Uses CAPD's double Taylor integrator
// (order 20, default tolerance 1e-18) on the 4D extended field (c' = 0), section x = 0, x' > 0.
//
// Construction (same as tang1/tang2/tang3.py of the c ~ 4.507 proof):
//   X0(s) = p_c + s e_u(c)          (s > 0: positive branch of W^u_loc)
//   eta(s) = P^K2(X0(s))            (K2 even, first fold C of W^u reached at s = s_C)
//   g_n(s,c) = u-coordinate of P^(K2+n)(X0(s)) - p_c in the eigenbasis E(c)
//   G_n(c) = extremum of g_n(.,c) near s_C ;  tangency  <=>  G_n(c) = 0.
//
// Modes:
//   num57 fp c
//   num57 scan c0 c1 nc K2 nmax             proxy: U_n(c) = g_n(s_C(c),c), n = 1..nmax
//   num57 refine cA cB n K2                 Brent root of G_n in [cA,cB]; prints proof data
//   num57 orbit c K2 nmax                   critical orbit and U_n at one c
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>
#include <functional>
#include <cstdlib>
#include "capd/capdlib.h"
using namespace capd;
using namespace std;

static double FLIP = 1;   // -1: use -e_u (other branch of W^u_loc, K2 odd)
static DMap* VF; static DOdeSolver* SOL; static DCoordinateSection* SEC; static DPoincareMap* PM;

struct V2 { double y, z; };
static V2 P(const V2& q, double c) {
  DVector u(4); u[0] = 0; u[1] = q.y; u[2] = q.z; u[3] = c;
  DVector r = (*PM)(u);
  return {r[1], r[2]};
}
// returns image, 2x2 derivative D, and d/dc
static V2 PD(const V2& q, double c, double D[2][2], double Dc[2], double* T = nullptr) {
  DVector u(4); u[0] = 0; u[1] = q.y; u[2] = q.z; u[3] = c;
  DMatrix M(4, 4); double t = 0;
  DVector r = (*PM)(u, M, t);
  DMatrix DP = PM->computeDP(r, M, t);
  for (int i = 0; i < 2; i++) { for (int k = 0; k < 2; k++) D[i][k] = DP[i + 1][k + 1]; Dc[i] = DP[i + 1][3]; }
  if (T) *T = t;
  return {r[1], r[2]};
}
static V2 Pn(V2 q, double c, int n) { for (int k = 0; k < n; k++) q = P(q, c); return q; }

struct FP {
  double c; V2 p; double D[2][2], Dc[2]; double lu, ls; double E[2][2], Ei[2][2]; double dp[2]; double T;
};
static FP fixedPoint(double c, V2 q) {
  FP F; F.c = c;
  double D[2][2], Dc[2];
  for (int it = 0; it < 50; it++) {
    V2 r = PD(q, c, D, Dc);
    double a = D[0][0] - 1, b = D[0][1], cc = D[1][0], d = D[1][1] - 1, det = a * d - b * cc;
    double fy = r.y - q.y, fz = r.z - q.z;
    double dy = -(d * fy - b * fz) / det, dz = -(-cc * fy + a * fz) / det;
    q.y += dy; q.z += dz;
    if (fabs(dy) + fabs(dz) < 1e-15) break;
  }
  F.p = q;
  PD(q, c, F.D, F.Dc, &F.T);
  double a = F.D[0][0], b = F.D[0][1], cc = F.D[1][0], d = F.D[1][1];
  double tr = a + d, det = a * d - b * cc, disc = sqrt(tr * tr / 4 - det);
  double l1 = tr / 2 + disc, l2 = tr / 2 - disc;
  F.lu = fabs(l1) > fabs(l2) ? l1 : l2; F.ls = det / F.lu;
  // eigenvectors, unit length; e_u with positive y-component, e_s with negative z-component
  auto evec = [&](double l, double v[2]) {
    double x1 = b, y1 = l - a, x2 = l - d, y2 = cc;
    if (hypot(x1, y1) > hypot(x2, y2)) { v[0] = x1; v[1] = y1; } else { v[0] = x2; v[1] = y2; }
    double n = hypot(v[0], v[1]); v[0] /= n; v[1] /= n;
  };
  double eu[2], es[2]; evec(F.lu, eu); evec(F.ls, es);
  if (eu[0] < 0) { eu[0] = -eu[0]; eu[1] = -eu[1]; }
  if (es[1] > 0) { es[0] = -es[0]; es[1] = -es[1]; }
  if (FLIP < 0) { eu[0] = -eu[0]; eu[1] = -eu[1]; }
  F.E[0][0] = eu[0]; F.E[1][0] = eu[1]; F.E[0][1] = es[0]; F.E[1][1] = es[1];
  double dE = F.E[0][0] * F.E[1][1] - F.E[0][1] * F.E[1][0];
  F.Ei[0][0] = F.E[1][1] / dE; F.Ei[0][1] = -F.E[0][1] / dE; F.Ei[1][0] = -F.E[1][0] / dE; F.Ei[1][1] = F.E[0][0] / dE;
  // dp/dc = -(D - I)^{-1} Dc
  double A = F.D[0][0] - 1, B = F.D[0][1], C = F.D[1][0], Dd = F.D[1][1] - 1, dt = A * Dd - B * C;
  F.dp[0] = -(Dd * F.Dc[0] - B * F.Dc[1]) / dt; F.dp[1] = -(-C * F.Dc[0] + A * F.Dc[1]) / dt;
  return F;
}
static void toUV(const FP& F, const V2& q, double& U, double& V) {
  double dy = q.y - F.p.y, dz = q.z - F.p.z;
  U = F.Ei[0][0] * dy + F.Ei[0][1] * dz; V = F.Ei[1][0] * dy + F.Ei[1][1] * dz;
}
static V2 X0(const FP& F, double s) { return {F.p.y + s * F.E[0][0], F.p.z + s * F.E[1][0]}; }

// Brent minimisation of f on [a,b] (golden section + parabolic), returns argmin
static double brentMin(const function<double(double)>& f, double a, double b, double tol, double& fx) {
  const double cg = 0.3819660112501051;
  double x = a + cg * (b - a), w = x, v = x, fw, fv, d = 0, e = 0;
  fx = f(x); fw = fv = fx;
  for (int it = 0; it < 200; it++) {
    double xm = 0.5 * (a + b), tol1 = tol * fabs(x) + 1e-300, tol2 = 2 * tol1;
    if (fabs(x - xm) <= tol2 - 0.5 * (b - a)) break;
    if (fabs(e) > tol1) {
      double r = (x - w) * (fx - fv), q = (x - v) * (fx - fw), p = (x - v) * q - (x - w) * r;
      q = 2 * (q - r); if (q > 0) p = -p; q = fabs(q);
      double etemp = e; e = d;
      if (fabs(p) >= fabs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x)) { e = (x >= xm) ? a - x : b - x; d = cg * e; }
      else { d = p / q; double u = x + d; if (u - a < tol2 || b - u < tol2) d = (xm - x >= 0) ? tol1 : -tol1; }
    } else { e = (x >= xm) ? a - x : b - x; d = cg * e; }
    double u = (fabs(d) >= tol1) ? x + d : x + (d >= 0 ? tol1 : -tol1);
    double fu = f(u);
    if (fu <= fx) { if (u >= x) a = x; else b = x; v = w; fv = fw; w = x; fw = fx; x = u; fx = fu; }
    else { if (u < x) a = u; else b = u;
      if (fu <= fw || w == x) { v = w; fv = fw; w = u; fw = fu; } else if (fu <= fv || v == x || v == w) { v = u; fv = fu; } }
  }
  return x;
}
// Brent root finding
static double brentRoot(const function<double(double)>& f, double a, double b, double tol) {
  double fa = f(a), fb = f(b);
  if (fa * fb > 0) { cerr << "brentRoot: no sign change " << fa << " " << fb << endl; return NAN; }
  double c = a, fc = fa, d = b - a, e = d;
  for (int it = 0; it < 200; it++) {
    if (fb * fc > 0) { c = a; fc = fa; d = b - a; e = d; }
    if (fabs(fc) < fabs(fb)) { a = b; b = c; c = a; fa = fb; fb = fc; fc = fa; }
    double tol1 = 2 * 1e-16 * fabs(b) + 0.5 * tol, xm = 0.5 * (c - b);
    if (fabs(xm) <= tol1 || fb == 0) return b;
    if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
      double s = fb / fa, p, q;
      if (a == c) { p = 2 * xm * s; q = 1 - s; }
      else { double qq = fa / fc, r = fb / fc; p = s * (2 * xm * qq * (qq - r) - (b - a) * (r - 1)); q = (qq - 1) * (r - 1) * (s - 1); }
      if (p > 0) q = -q; p = fabs(p);
      if (2 * p < std::min(3 * xm * q - fabs(tol1 * q), fabs(e * q))) { e = d; d = p / q; } else { d = xm; e = d; }
    } else { d = xm; e = d; }
    a = b; fa = fb;
    b += (fabs(d) > tol1) ? d : (xm > 0 ? tol1 : -tol1);
    fb = f(b);
  }
  return b;
}

// locate s_C: first local minimum of y(P^(K2+1)(X0(s))) for s in (0, smax]
static double findSC(const FP& F, int K2, double guess) {
  double c = F.c;
  auto y1 = [&](double s) { return Pn(X0(F, s), c, K2 + 1).y; };
  double lo, hi;
  if (guess > 0) { lo = guess * 0.97; hi = guess * 1.03; }
  else {
    double smax = 3.0 / pow(fabs(F.lu), K2);
    int M = 120; double prev = y1(smax / M), cur; int found = -1; vector<double> ys(M + 1); ys[1] = prev;
    for (int k = 2; k <= M; k++) { ys[k] = y1(k * smax / M); }
    for (int k = 2; k < M; k++) if (ys[k] < ys[k - 1] && ys[k] <= ys[k + 1]) { found = k; break; }
    if (found < 0) { cerr << "findSC: no fold found" << endl; return NAN; }
    lo = (found - 1) * smax / M; hi = (found + 1) * smax / M; (void)cur;
  }
  double fx;
  return brentMin(y1, lo, hi, 1e-13, fx);
}

// G_n(c): extremum of g_n near s_C; returns G, s*, sign (+1 max, -1 min)
struct GRes { double G, s, sign, V; };
static GRes Gn(double c, int n, int K2, V2 pguess, double sguess, FP* Fout = nullptr) {
  FP F = fixedPoint(c, pguess);
  if (Fout) *Fout = F;
  double sC = findSC(F, K2, sguess);
  auto g = [&](double s) { double U, V; toUV(F, Pn(X0(F, s), c, K2 + n), U, V); return U; };
  double h = sC * 2e-3;
  double gm = g(sC), gl = g(sC - h), gr = g(sC + h);
  double sign = (gm > 0.5 * (gl + gr)) ? 1 : -1;   // max if curve bends down
  double fx;
  double s = brentMin([&](double x) { return -sign * g(x); }, sC - h, sC + h, 1e-12, fx);
  double U, V; toUV(F, Pn(X0(F, s), c, K2 + n), U, V);
  return {U, s, sign, V};
}

int main(int argc, char* argv[]) {
  cout.precision(16);
  VF = new DMap("par:a,b;var:x,y,z,c;fun:-(y+z),x+a*y,b+z*(x-c),0;");
  VF->setParameter("a", 0.2); VF->setParameter("b", 0.2);
  SOL = new DOdeSolver(*VF, 20);
  SEC = new DCoordinateSection(4, 0);
  PM = new DPoincareMap(*SOL, *SEC, poincare::MinusPlus);
  string mode = argc > 1 ? argv[1] : "fp";
  V2 pg = {-8.381, 0.0296};
  if (mode == "fp") {
    double c = atof(argv[2]);
    FP F = fixedPoint(c, pg);
    cout << "c=" << c << " p=(" << F.p.y << "," << F.p.z << ") lu=" << F.lu << " ls=" << F.ls << " T=" << F.T << endl;
    cout << "E=[[" << F.E[0][0] << "," << F.E[0][1] << "],[" << F.E[1][0] << "," << F.E[1][1] << "]]  dp/dc=(" << F.dp[0] << "," << F.dp[1] << ")" << endl;
  } else if (mode == "scan") {
    double c0 = atof(argv[2]), c1 = atof(argv[3]); int nc = atoi(argv[4]), K2 = atoi(argv[5]), nmax = atoi(argv[6]);
    double sg = -1; V2 p = pg;
    cout << "# c  lu  sC  U_1..U_nmax  (U_n = u-coord of P^(K2+n)(X0(s_C)) - p)" << endl;
    for (int i = 0; i <= nc; i++) {
      double c = c0 + (c1 - c0) * i / nc;
      FP F = fixedPoint(c, p); p = F.p;
      double sC = findSC(F, K2, sg);
      if (std::isnan(sC)) { sg = -1; continue; }
      // continuation guard: if the fold jumped to the bracket edge, redo the full search
      if (sg > 0 && (fabs(sC / sg - 1) > 0.029)) sC = findSC(F, K2, -1);
      sg = sC;
      V2 q = Pn(X0(F, sC), c, K2);
      cout << c << " " << F.lu << " " << sC;
      for (int n = 1; n <= nmax; n++) { q = P(q, c); double U, V; toUV(F, q, U, V); cout << " " << U; }
      cout << endl;
    }
  } else if (mode == "refine") {
    double cA = atof(argv[2]), cB = atof(argv[3]); int n = atoi(argv[4]), K2 = atoi(argv[5]);
    if (argc > 6) FLIP = atof(argv[6]);
    FP F0 = fixedPoint(0.5 * (cA + cB), pg);
    double sg = findSC(F0, K2, -1);
    auto f = [&](double c) { return Gn(c, n, K2, F0.p, sg).G; };
    double cs = brentRoot(f, cA, cB, 1e-15);
    FP F; GRes R = Gn(cs, n, K2, F0.p, sg, &F);
    double h = 1e-7;
    double dG = (f(cs + h) - f(cs - h)) / (2 * h);
    double dG2 = (f(cs + h / 10) - f(cs - h / 10)) / (2 * h / 10);
    auto g = [&](double s) { double U, V; toUV(F, Pn(X0(F, s), cs, K2 + n), U, V); return U; };
    double ds = R.s * 1e-4;
    double gss = (g(R.s + ds) - 2 * g(R.s) + g(R.s - ds)) / (ds * ds);
    double gss2 = (g(R.s + ds / 3) - 2 * g(R.s) + g(R.s - ds / 3)) / (ds * ds / 9);
    cout << "n=" << n << " K2=" << K2 << " N=" << K2 + n << " flip=" << FLIP << endl;
    cout << setprecision(17);
    cout << "c0=" << cs << "  c0-5.7=" << cs - 5.7 << endl;
    cout << "G(c0)=" << R.G << "  ext=" << (R.sign > 0 ? "max" : "min") << "  s*=" << R.s << endl;
    cout << "dG/dc=" << dG << " (h/10: " << dG2 << ")   g''=" << gss << " (" << gss2 << ")" << endl;
    cout << "V(tip)=" << R.V << endl;
    cout << "p=(" << F.p.y << "," << F.p.z << ") lu=" << F.lu << " ls=" << F.ls << " det=" << F.lu * F.ls << " T=" << F.T << endl;
    cout << "E=[[" << F.E[0][0] << "," << F.E[0][1] << "],[" << F.E[1][0] << "," << F.E[1][1] << "]]" << endl;
    cout << "dp/dc=(" << F.dp[0] << "," << F.dp[1] << ")" << endl;
    // critical orbit
    V2 q = Pn(X0(F, R.s), cs, K2);
    cout << "C=(" << q.y << "," << q.z << ")  orbit y:";
    for (int k = 1; k <= n; k++) { q = P(q, cs); cout << " " << q.y; }
    cout << endl;
  } else if (mode == "orbit") {
    double c = atof(argv[2]); int K2 = atoi(argv[3]), nmax = atoi(argv[4]);
    FP F = fixedPoint(c, pg);
    double sC = findSC(F, K2, -1);
    V2 q = Pn(X0(F, sC), c, K2);
    cout << "c=" << c << " p.y=" << F.p.y << " lu=" << F.lu << " sC=" << sC << " C=(" << q.y << "," << q.z << ")" << endl;
    for (int n = 1; n <= nmax; n++) { q = P(q, c); double U, V; toUV(F, q, U, V); cout << n << " y=" << q.y << " z=" << q.z << " U=" << U << " V=" << V << endl; }
  }
  return 0;
}
