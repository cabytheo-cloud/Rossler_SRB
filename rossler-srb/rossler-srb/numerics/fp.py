import numpy as np, sys
from pmap import P
def attr_pt(c,n=30):
    q=np.array([-6.,0.03])
    for k in range(n): q=P(q,c)
    return q
def fixed_point(c,q=None,per=1):
    if q is None:
        q0=attr_pt(c)
        ys=np.linspace(-11,-2.5,60); vals=[]
        for y in ys:
            try: vals.append(P(np.array([y,q0[1]]),c)[0]-y)
            except Exception: vals.append(np.nan)
        vals=np.array(vals); i=np.where(np.sign(vals[:-1])*np.sign(vals[1:])<0)[0]
        q=np.array([ys[i[0]],q0[1]])
    for it in range(40):
        Q=q.copy(); D=np.eye(2)
        for k in range(per):
            Q,Dk,_=P(Q,c,var=True); D=Dk@D
        F=Q-q
        dq=np.linalg.solve(D-np.eye(2),-F); q=q+dq
        if np.linalg.norm(dq)<1e-14: break
    return q,D
if __name__=="__main__":
    for c in [4.0,4.2,4.4,4.6,5.0,5.7]:
        q,D=fixed_point(c)
        w,V=np.linalg.eig(D)
        print(c,q,w,np.linalg.det(D)); print("  vecs",V.T)
