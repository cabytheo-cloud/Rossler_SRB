// Verify P(W) subset int(W) for the DECIMAL box W = [-10.7,-2.7] x [0.028,0.034]
// (a = b = 0.2, c = 5.7; Rossler, section x = 0, x' > 0). Decimal endpoints are enclosed by
// exact rational interval arithmetic; inputs cover W outwardly, the inclusion test is inward.
#include <iostream>
#include "capd/capdlib.h"
using namespace capd; using namespace std; typedef interval I;
int main(){
  cout.precision(17); cout<<boolalpha;
  IMap vf("par:a,b,c;var:x,y,z;fun:-(y+z),x+a*y,b+z*(x-c);");
  vf.setParameter("a",I(2)/I(10)); vf.setParameter("b",I(2)/I(10)); vf.setParameter("c",I(57)/I(10));
  IOdeSolver solver(vf,20); ICoordinateSection sec(3,0); IPoincareMap pm(solver,sec,poincare::MinusPlus);
  I YL=I(-107)/I(10), YR=I(-27)/I(10), ZL=I(28)/I(1000), ZR=I(34)/I(1000);
  double yl=YL.leftBound(), yr=YR.rightBound(), zl=ZL.leftBound(), zr=ZR.rightBound(); // outer box
  int Nb=400; bool ok=true;
  double ylo=1e9, yhi=-1e9, zlo=1e9, zhi=-1e9;
  for(int i=0;i<Nb;i++){
    I y = I(yl) + (I(yr)-I(yl))*I(i,i+1)/I(Nb);
    if(i==Nb-1) y = I(y.leftBound(), yr);
    if(i==0) y = I(yl, y.rightBound());
    IVector x{I(0.),y,I(zl,zr)};
    C0HOTripletonSet s(x);
    IVector P=pm(s);
    ylo=std::min(ylo,(double)P[1].leftBound()); yhi=std::max(yhi,(double)P[1].rightBound());
    zlo=std::min(zlo,(double)P[2].leftBound()); zhi=std::max(zhi,(double)P[2].rightBound());
    ok = ok && P[1].leftBound()>YL.rightBound() && P[1].rightBound()<YR.leftBound()
            && P[2].leftBound()>ZL.rightBound() && P[2].rightBound()<ZR.leftBound();
  }
  cout<<"P(W) in ["<<ylo<<","<<yhi<<"]x["<<zlo<<","<<zhi<<"]\nP(W) subset int W (decimal W): "<<ok<<endl;
}
