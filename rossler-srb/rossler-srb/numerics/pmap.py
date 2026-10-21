import numpy as np
from scipy.integrate import solve_ivp
a=b=0.2
def rhs(t,u,c):
    x,y,z=u[0],u[1],u[2]
    return [-y-z, x+a*y, b+z*(x-c)]
def rhs_var(t,u,c):
    x,y,z=u[:3]; V=u[3:].reshape(3,3)
    J=np.array([[0,-1,-1],[1,a,0],[z,0,x-c]])
    return np.concatenate([[-y-z,x+a*y,b+z*(x-c)],(J@V).ravel()])
def ev(t,u,c): return u[0]
ev.direction=1; ev.terminal=True
def P(q,c,var=False,rtol=1e-13,atol=1e-15):
    """Poincare map on x=0, xdot>0. q=(y,z)."""
    y,z=q
    u0=np.array([0.,y,z])
    if var:
        u0=np.concatenate([u0,np.eye(3).ravel()]); f=rhs_var
    else: f=rhs
    # skip initial crossing: integrate a bit first
    s=solve_ivp(f,(0,0.5),u0,args=(c,),method='DOP853',rtol=rtol,atol=atol)
    u1=s.y[:,-1]
    s=solve_ivp(f,(0.5,100),u1,args=(c,),method='DOP853',rtol=rtol,atol=atol,events=ev)
    ue=s.y_events[0][0]; T=s.t_events[0][0]
    if not var: return ue[1:3]
    Phi=ue[3:].reshape(3,3); fv=np.array(rhs(0,ue[:3],c))
    D=Phi-np.outer(fv,Phi[0])/fv[0]
    return ue[1:3], D[1:,1:], T
