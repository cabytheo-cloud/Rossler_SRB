import numpy as np
from scipy.optimize import minimize_scalar, brentq
from pmap import P
from fp import fixed_point
def setup(c):
    p,D=fixed_point(c, q=np.array([-6.9,0.0365]) if c<4.8 else None)
    w,V=np.linalg.eig(D); i=np.argmax(abs(w))
    eu=V[:,i].real; es=V[:,1-i].real
    if eu[0]<0: eu=-eu
    E=np.column_stack([eu,es]); Ei=np.linalg.inv(E)
    return p,w[i].real,w[1-i].real,E,Ei
def G(c,M=None,verbose=False):
    p,lu,ls,E,Ei=setup(c)
    # critical point of 1D map along attractor: points P(p + t*eu) lie on W^u
    def F(t):  # y of P^2(p+t eu)?? use curve eta(t)=P(p + t*eu)/ just param by y on W^u approx: point (y, z) with z from W^u
        q=p+t*E[:,0]
        return P(q,c)[0]
    # minimize y-image over t (critical point on W^u_loc-extension, approx linear)
    r=minimize_scalar(F,bracket=(0.5,1.4,2.5),tol=1e-12)
    tC=r.x; qC=p+tC*E[:,0]
    v=qC
    for j in range(3): v=P(v,c)
    uv=Ei@(v-p)
    if verbose: print("c",c,"p",p,"lu",lu,"ls",ls,"tC",tC,"v1",P(qC,c),"v3 uv",uv)
    return uv[0]
if __name__=="__main__":
    for c in np.linspace(4.40,4.60,11):
        print(c,G(c,verbose=False))
