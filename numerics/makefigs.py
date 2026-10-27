import json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"serif","mathtext.fontset":"cm","font.size":9,
    "axes.edgecolor":"#8a8985","axes.linewidth":0.6,"xtick.color":"#52514e","ytick.color":"#52514e",
    "axes.labelcolor":"#0b0b0b","xtick.major.width":0.6,"ytick.major.width":0.6})
BLUE,ORANGE,AQUA,GRAY,INK="#2a78d6","#eb6834","#1baf7a","#b9b8b3","#0b0b0b"
d=json.load(open("paper/figs/fig1.json"))
pairs=np.array(d["pairs"]); crit=np.array(d["crit"]); p=d["p"]
fig,ax=plt.subplots(figsize=(4.6,3.6))
lo,hi=-8.9,-3.7
ax.plot([lo,hi],[lo,hi],color=GRAY,lw=0.8,zorder=1)
ax.scatter(pairs[:,0],pairs[:,2],s=1.2,color=BLUE,lw=0,zorder=2)
# cobweb of critical orbit
xs=crit[:,0]
cx,cy=[xs[0]],[xs[0]]
path_x=[xs[0]]; path_y=[lo]
for k in range(3):
    path_x+= [xs[k],xs[k+1]]; path_y+=[xs[k+1],xs[k+1]]
ax.plot([xs[0],xs[0]],[lo,xs[1]],color=ORANGE,lw=1.2,zorder=3)
for k in range(1,3):
    ax.plot([xs[k-1],xs[k]],[xs[k],xs[k]],color=ORANGE,lw=1.2,zorder=3)
    ax.plot([xs[k],xs[k]],[xs[k],xs[k+1]],color=ORANGE,lw=1.2,zorder=3)
ax.plot([xs[2],xs[3]],[xs[3],xs[3]],color=ORANGE,lw=1.2,zorder=3)
ax.scatter([p[0]],[p[0]],s=28,color=INK,zorder=5,edgecolor="white",linewidth=1.2)
ax.annotate(r"$p$",(p[0],p[0]),xytext=(8,-10),textcoords="offset points",color=INK)
ax.scatter([xs[0]],[xs[1]],s=16,color=ORANGE,zorder=5,edgecolor="white",linewidth=1.0)
ax.annotate(r"$(C,\,v_1)$",(xs[0],xs[1]),xytext=(6,-11),textcoords="offset points",color=INK)
ax.annotate(r"$v_2$",(xs[2],xs[2]),xytext=(-14,4),textcoords="offset points",color=INK)
ax.annotate(r"$v_3=P^3(C)\in W^s(p)$",(xs[3],xs[3]),xytext=(-70,22),textcoords="offset points",color=INK,
            arrowprops=dict(arrowstyle="-",color="#8a8985",lw=0.6))
ax.set_xlim(lo,hi); ax.set_ylim(lo,hi)
ax.set_xlabel(r"$y_n$"); ax.set_ylabel(r"$y_{n+1}$")
ax.set_aspect("equal")
for s in ["top","right"]: ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig("paper/figs/returnmap.pdf"); fig.savefig("paper/figs/returnmap.png",dpi=160)
# ---- fig 2
d2=json.load(open("paper/figs/fig2.json"))
fig,ax=plt.subplots(figsize=(4.6,3.0))
ax.axhline(0,color=GRAY,lw=0.8,zorder=1)
labels={"-5e-05":(r"$c=c_*-5\cdot10^{-5}$",AQUA),"0.0":(r"$c\approx c_*$",INK),"5e-05":(r"$c=c_*+5\cdot10^{-5}$",ORANGE)}
for k,(lab,col) in labels.items():
    a=np.array(d2[k]); ax.plot(a[:,0]*1e6,a[:,1]*1e4,color=col,lw=1.6,zorder=3)
    ax.annotate(lab,(a[-1,0]*1e6,a[-1,1]*1e4),xytext=(5,0),textcoords="offset points",color=INK,fontsize=8,va="center")
ax.annotate(r"$\{g=0\}$: $\Gamma(s,c)\in W^s_{\rm loc}(p_c)$",(7.15,0),xytext=(0,4),textcoords="offset points",ha="right",color="#52514e",fontsize=8)
ax.set_xlim(-4.5,7.2)
ax.set_xlabel(r"$s-s_*\ \ (\times 10^{-6})$"); ax.set_ylabel(r"$g(s,c)\ \ (\times 10^{-4})$")
for s in ["top","right"]: ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig("paper/figs/unfolding.pdf"); fig.savefig("paper/figs/unfolding.png",dpi=160)
print("ok")
