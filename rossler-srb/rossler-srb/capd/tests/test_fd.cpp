// Non-rigorous cross-check of the derivatives used by the proof programs: the first and second
// derivatives of the N-th return map P_c^N (CAPD computeDP, with return-time corrections, for the
// 4D extended field with c as a state variable) are compared with central finite differences of
// P_c^N itself, at the point X0 = p~ + s_m e_u used in the proof, in the directions e_u, e_s and c.
//
// usage: ./test_fd cmid py pz E00 E01 E10 E11 sm N
//   (defaults: the band-merging case J_0 of Proposition 6.1)
#include <iostream>
#include <iomanip>
#include "capd/capdlib.h"
using namespace capd;
using namespace std;
typedef interval I;

static IC2PoincareMap* PM;
static double MAXW = 0;   // largest width of the point enclosures of P^N (noise level of the FD)

struct Pt { double v[2]; };
static Pt Pn(double y, double z, double c, int N) {   // midpoint of a rigorous enclosure of P_c^N(y,z)
  IVector u(4); u[0] = 0.; u[1] = y; u[2] = z; u[3] = c;
  C2Rect2Set s(u);
  IMatrix dF(4, 4); IHessian H(4, 4); I T;
  IVector P = (*PM)(s, dF, H, T, N);
  Pt r; r.v[0] = P[1].mid().leftBound(); r.v[1] = P[2].mid().leftBound();
  MAXW = std::max(MAXW, std::max(P[1].rightBound() - P[1].leftBound(), P[2].rightBound() - P[2].leftBound()));
  return r;
}

int main(int argc, char* argv[]) {
  cout << std::scientific; cout.precision(9);
  double cmid = 4.50698934631197, py = -6.938698731191516, pz = 0.03606348866252088;
  double E00 = 0.9999996961880596, E01 = -0.8077389541129416, E10 = 7.795022697815469e-4, E11 = -0.5895403141503822;
  double sm = 5.3682155586e-4; int N = 15;
  if (argc > 9) { cmid = atof(argv[1]); py = atof(argv[2]); pz = atof(argv[3]); E00 = atof(argv[4]); E01 = atof(argv[5]);
                  E10 = atof(argv[6]); E11 = atof(argv[7]); sm = atof(argv[8]); N = atoi(argv[9]); }
  IMap vf("par:a,b;var:x,y,z,c;fun:-(y+z),x+a*y,b+z*(x-c),0;");
  vf.setParameter("a", I(2) / I(10)); vf.setParameter("b", I(2) / I(10));
  IC2OdeSolver solver(vf, 20);
  ICoordinateSection section(4, 0);
  IC2PoincareMap pm(solver, section, poincare::MinusPlus);
  PM = &pm;

  // directions (in (y, z, c)) and two base points: X0 = p~ + s_m e_u (the fold, used in the proof)
  // and X1 = p~ + (s_m/2) e_u (off the fold, where D P^N e_u is not small)
  double dir[3][3] = {{E00, E10, 0.}, {E01, E11, 0.}, {0., 0., 1.}};
  const char* nm[3] = {"e_u", "e_s", "c"};
  double worst = 0;
  for (int base = 0; base < 2; base++) {
    double ss = (base == 0) ? sm : 0.5 * sm;
    double x0[3] = {py + ss * E00, pz + ss * E10, cmid};
    IVector u(4); u[0] = 0.; u[1] = x0[0]; u[2] = x0[1]; u[3] = x0[2];
    C2Rect2Set s(u);
    IMatrix dF(4, 4); IHessian Hf(4, 4); I T;
    IVector P = pm(s, dF, Hf, T, N);
    IMatrix DP(4, 4); IHessian D2P(4, 4);
    pm.computeDP(P, dF, Hf, DP, D2P, T);
    auto D1 = [&](int i, int d) { I r = 0.; for (int a = 0; a < 3; a++) r += DP[i + 1][a + 1] * dir[d][a]; return r.mid().leftBound(); };
    auto D2 = [&](int i, int d1, int d2) {
      I r = 0.;
      for (int a = 0; a < 3; a++) for (int b = 0; b < 3; b++)
        r += (a == b ? I(2.) : I(1.)) * D2P(i + 1, a + 1, b + 1) * dir[d1][a] * dir[d2][b];
      return r.mid().leftBound(); };
    auto F = [&](double h1, int d1, double h2, int d2) {
      return Pn(x0[0] + h1 * dir[d1][0] + h2 * dir[d2][0], x0[1] + h1 * dir[d1][1] + h2 * dir[d2][1],
                x0[2] + h1 * dir[d1][2] + h2 * dir[d2][2], N); };
    MAXW = 0; (void)F(0., 0, 0., 0);
    const double eps = std::max(2e-15, MAXW);   // noise level of the point values
    cout << "N = " << N << ", c = " << cmid << ", base point " << (base == 0 ? "X0 (fold)" : "X1 (off fold)")
         << " = (" << x0[0] << ", " << x0[1] << ")" << endl;
    for (int d = 0; d < 3; d++) for (int i = 0; i < 2; i++) {
      double jd = D1(i, d), best = 1e300, bestfd = 0, besth = 0;
      for (double h : {1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11}) {
        Pt fp = F(h, d, 0., d), fm = F(-h, d, 0., d);
        double fd = (fp.v[i] - fm.v[i]) / (2 * h);
        if (fabs(fd - jd) < best) { best = fabs(fd - jd); bestfd = fd; besth = h; }
      }
      double noise = eps / besth;
      bool resolved = fabs(jd) > 100 * noise;
      double rel = best / std::max(fabs(jd), 1e-300);
      if (resolved) worst = std::max(worst, rel);
      cout << "  D P^N[" << nm[d] << "]_" << i << ": computeDP " << setw(17) << jd << "  FD " << setw(17) << bestfd
           << "  (h=" << besth << ")  rel.diff " << setw(10) << rel << "  rel.noise " << setw(10) << noise / std::max(fabs(jd), 1e-300)
           << (resolved ? "" : "  [not resolved]") << endl;
    }
    for (int d1 = 0; d1 < 3; d1++) for (int d2 = d1; d2 < 3; d2++) for (int i = 0; i < 2; i++) {
      double jd = D2(i, d1, d2), best = 1e300, bestfd = 0, besth = 0;
      for (double h : {1e-4, 1e-5, 1e-6, 1e-7, 1e-8}) {
        Pt pp = F(h, d1, h, d2), pmn = F(h, d1, -h, d2), mp = F(-h, d1, h, d2), mm = F(-h, d1, -h, d2);
        double fd = (pp.v[i] - pmn.v[i] - mp.v[i] + mm.v[i]) / (4 * h * h);
        if (fabs(fd - jd) < best) { best = fabs(fd - jd); bestfd = fd; besth = h; }
      }
      double noise = eps / (besth * besth);
      bool resolved = fabs(jd) > 100 * noise;
      double rel = best / std::max(fabs(jd), 1e-300);
      if (resolved) worst = std::max(worst, rel);
      cout << "  D2P^N[" << nm[d1] << "," << nm[d2] << "]_" << i << ": computeDP " << setw(17) << jd << "  FD " << setw(17) << bestfd
           << "  (h=" << besth << ")  rel.diff " << setw(10) << rel << "  rel.noise " << setw(10) << noise / std::max(fabs(jd), 1e-300)
           << (resolved ? "" : "  [not resolved]") << endl;
    }
  }
  cout << "largest relative difference over all resolved entries (value > 100 x noise): " << worst << endl;
  cout << "(rel.noise = (width of the point enclosures of P^N)/h^k, relative to the value; the truncation error of the\n"
          " finite differences is not estimated, so very large higher derivatives can also cause disagreement)" << endl;
  return 0;
}
