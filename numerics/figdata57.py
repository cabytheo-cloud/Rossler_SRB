import numpy as np, json
from pmap import P
P3=json.load(open('near57/params_n16a.json'))
c0=float(P3['cmid']); K2=P3['K2']; n=P3['n']; sm=float(P3['smid'])
E=np.array([[float(P3['E00']),float(P3['E01'])],[float(P3['E10']),float(P3['E11'])]]); Ei=np.linalg.inv(E)
def fixed_point(c,q):
    for it in range(30):
        Q,D,_=P(q,c,var=True); dq=np.linalg.solve(D-np.eye(2),-(Q-q)); q=q+dq
        if np.linalg.norm(dq)<1e-15: break
    return q
p0=fixed_point(c0,np.array([float(P3['py']),float(P3['pz'])]))
# --- fig 1: orbit on the attractor and the critical orbit at c0
q=np.array([-6.,0.03])
for k in range(50): q=P(q,c0)
pairs=[]
for k in range(1500):
    q2=P(q,c0); pairs.append((q[0],q2[0])); q=q2
x=p0+sm*E[:,0]
for k in range(K2): x=P(x,c0)
crit=[x.tolist()]
for k in range(n): x=P(x,c0); crit.append(x.tolist())
print("p0",p0,"C",crit[0],"v16",crit[-1],"v16-p",np.array(crit[-1])-p0)
json.dump({"pairs":pairs,"crit":crit,"p":p0.tolist(),"c0":c0},open('paper/figs/fig1_57.json','w'))
# --- fig 2: unfolding near c0
res={}
for dc in [-5e-12,0.0,5e-12]:
    c=c0+dc; pc=fixed_point(c,p0)
    out=[]
    for s in np.linspace(sm-1.5e-9,sm+1.5e-9,31):
        y=pc+s*E[:,0]
        for k in range(K2+n): y=P(y,c)
        out.append((float(s-sm),float((Ei@(y-pc))[0])))
    res[str(dc)]=out; print("dc",dc,"max g",max(o[1] for o in out))
json.dump(res,open('paper/figs/fig2_57.json','w'))
