#include <iostream>
#include "capd/capdlib.h"
using namespace capd; using namespace std;
int main(){
  cout.precision(16);
  IMap vf("par:a,b;var:x,y,z,c;fun:-(y+z),x+a*y,b+z*(x-c),0;");
  vf.setParameter("a",interval(2)/interval(10)); vf.setParameter("b",interval(2)/interval(10));
  IC2OdeSolver solver(vf,20);
  ICoordinateSection section(4,0);
  IC2PoincareMap pm(solver,section,poincare::MinusPlus);
  IVector u{interval(0.),interval(-6.9),interval(0.0365),interval(4.5)};
  C2Rect2Set s(u);
  IMatrix dF(4,4); IHessian H(4,4); interval T;
  IVector P=pm(s,dF,H,T,1);
  IMatrix DP(4,4); IHessian D2P(4,4);
  pm.computeDP(P,dF,H,DP,D2P,T);
  cout<<"P="<<P<<"\nDP="<<DP<<endl;
  for(int i=1;i<3;i++) for(int j=1;j<4;j++) for(int k=1;k<4;k++) cout<<"D2P("<<i<<","<<j<<","<<k<<")="<<D2P(i,j,k)<<endl;
}
