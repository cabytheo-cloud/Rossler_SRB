#include <iostream>
#include "capd/capdlib.h"
using namespace capd; using namespace std;
typedef interval I;
int main(){
  cout.precision(16); cout<<boolalpha;
  I x(-2.,-1.), y(-1., 1.), z(0.5, 2.);
  cout<<"x<0 "<<(x<0.)<<"  y<0 "<<(y<0.)<<"  y>0 "<<(y>0.)<<"  z>0 "<<(z>0.)<<"  x<-1 "<<(x< -1.)<<endl;
  cout<<"subset(y,sym) "<<subset(I(-0.5,0.5),y)<<" "<<subset(I(-1.5,0.5),y)<<endl;
  IVector a{I(0.,1.)}, b{I(-0.1,1.1)}, b2{I(0.,1.1)};
  cout<<"subsetInterior "<<subsetInterior(a,b)<<" "<<subsetInterior(a,b2)<<endl;
  // derivative of P^2 vs composition
  IMap vf("par:a,b;var:x,y,z,c;fun:-(y+z),x+a*y,b+z*(x-c),0;");
  vf.setParameter("a",I(2)/I(10)); vf.setParameter("b",I(2)/I(10));
  IC2OdeSolver solver(vf,20); ICoordinateSection sec(4,0); IC2PoincareMap pm(solver,sec,poincare::MinusPlus);
  IVector u{I(0.),I(-5.5),I(0.0372),I(4.507)};
  auto jet=[&](IVector v,int n, IMatrix& DP, IHessian& D2P){ C2Rect2Set s(v); IMatrix dF(4,4); IHessian H(4,4); I T; IVector P=pm(s,dF,H,T,n); pm.computeDP(P,dF,H,DP,D2P,T); return P;};
  IMatrix D1(4,4),D2(4,4),D12(4,4); IHessian H1(4,4),H2(4,4),H12(4,4);
  IVector P1=jet(u,1,D1,H1);
  IVector P1m=midVector(P1);
  IVector P2=jet(P1m,1,D2,H2);
  IVector P12=jet(u,2,D12,H12);
  cout<<"P^2 direct "<<P12<<"\nP(P) "<<P2<<endl;
  cout<<"D(P^2) direct\n"<<D12<<"\nproduct\n"<<D2*D1<<endl;
}
