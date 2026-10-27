import numpy as np, json
from pmap import P
from tang2 import make_g
c0=4.50698934631197
# --- Fig 1 data: orbit on the section and critical orbit
p,lu,ls,E,Ei,eta,g=make_g(c0)
q=np.array([-6.,0.037])
for k in range(50): q=P(q,c0)
Y=[]
for k in range(2500):
    q=P(q,c0); Y.append(q[0])
s_star=5.3682155586e-04
C=eta(s_star); crit=[C]
for j in range(3): crit.append(P(crit[-1],c0))
# 1D map curve: images of points on W^u (attractor) parametrized by y: use eta(s) for s in range to cover [v2,v1]?
# simpler: take P of points (y, z) with z on attractor via short local curve of W^u: P^{12}(p + s eu) for s in (0, s_max)
ys=np.linspace(min(Y),max(Y),400)
# attractor z as function of y: from orbit points
Yarr=np.array(Y); 
out={"orbit_y":Y,"crit":[list(map(float,x)) for x in crit],"p":list(map(float,p))}
# map curve: for each y in grid, find attractor z by interpolation over orbit pairs
zs=[]
q=np.array([-6.,0.037])
pts=[]
for k in range(50): q=P(q,c0)
for k in range(400):
    q2=P(q,c0); pts.append((q[0],q[1],q2[0])); q=q2
out["pairs"]=pts
json.dump(out,open("paper/figs/fig1.json","w"))
print("fig1 done", p, [c.tolist() for c in crit])
# --- Fig 2 data: unfolding g(s,c)
res={}
for dc in [-5e-5,0.0,5e-5]:
    c=c0+dc
    p,lu,ls,E,Ei,eta,g=make_g(c)
    # s*(c) ~ s_star + ds*/dc * dc, ds*/dc ~ -1.46e-3
    sc=s_star-1.46e-3*dc
    ss=np.linspace(sc-4e-6,sc+4e-6,33)
    res[str(dc)]=[(float(s-s_star),float(g(s))) for s in ss]
    print("dc",dc,"done")
json.dump(res,open("paper/figs/fig2.json","w"))
