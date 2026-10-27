import json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"serif","mathtext.fontset":"cm","font.size":9,
    "axes.edgecolor":"#8a8985","axes.linewidth":0.6,"xtick.color":"#52514e","ytick.color":"#52514e",
    "axes.labelcolor":"#0b0b0b","xtick.major.width":0.6,"ytick.major.width":0.6})
BLUE,ORANGE,AQUA,GRAY,INK="#2a78d6","#eb6834","#1baf7a","#b9b8b3","#0b0b0b"
d=json.load(open("paper/figs/fig1_57.json"))
pairs=np.array(d["pairs"]); crit=np.array(d["crit"]); p=d["p"]; xs=crit[:,0]
fig,ax=plt.subplots(figsize=(4.6,3.9))
lo,hi=-10.9,-2.8
ax.plot([lo,hi],[lo,hi],color=GRAY,lw=0.8,zorder=1)
ax.scatter(pairs[:,0],pairs[:,1],s=1.0,color=BLUE,lw=0,zorder=2)
for k in range(len(xs)-1):
    if k>0: ax.plot([xs[k-1],xs[k]],[xs[k],xs[k]],color=ORANGE,lw=0.7,alpha=0.75,zorder=3)
    ax.plot([xs[k],xs[k]],[xs[k],xs[k+1]] if k>0 else [xs[1],xs[1]],color=ORANGE,lw=0.7,alpha=0.75,zorder=3)
ax.plot([xs[-2],xs[-1]],[xs[-1],xs[-1]],color=ORANGE,lw=0.7,alpha=0.75,zorder=3)
ax.scatter([xs[0]],[xs[1]],s=16,color=ORANGE,zorder=5,edgecolor="white",linewidth=1.0)
ax.annotate(r"$(C,\,v_1)$",(xs[0],xs[1]),xytext=(6,-11),textcoords="offset points",color=INK)
ax.scatter([p[0]],[p[0]],s=28,color=INK,zorder=6,edgecolor="white",linewidth=1.2)
ax.annotate(r"$p$",(p[0],p[0]),xytext=(-12,6),textcoords="offset points",color=INK)
ax.annotate(r"$v_{16}=P^{16}(C)\in W^s(p)$",(p[0],p[0]),xytext=(18,-34),textcoords="offset points",color=INK,
            arrowprops=dict(arrowstyle="-",color="#8a8985",lw=0.6))
ax.set_xlim(lo,hi); ax.set_ylim(lo,hi); ax.set_aspect("equal")
ax.set_xlabel(r"$y_n$"); ax.set_ylabel(r"$y_{n+1}$")
for s in ["top","right"]: ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig("paper/returnmap57.pdf"); fig.savefig("paper/figs/returnmap57.png",dpi=160)
d2=json.load(open("paper/figs/fig2_57.json"))
fig,ax=plt.subplots(figsize=(4.6,3.0))
ax.axhline(0,color=GRAY,lw=0.8,zorder=1)
labels={"-5e-12":(r"$c=c_0-5\cdot10^{-12}$",AQUA),"0.0":(r"$c=c_0$",INK),"5e-12":(r"$c=c_0+5\cdot10^{-12}$",ORANGE)}
for k,(lab,col) in labels.items():
    a=np.array(d2[k]); ax.plot(a[:,0]*1e9,a[:,1]*1e7,color=col,lw=1.6,zorder=3)
    ax.annotate(lab,(a[-1,0]*1e9,a[-1,1]*1e7),xytext=(5,0),textcoords="offset points",color=INK,fontsize=8,va="center")
ax.annotate(r"$\{g=0\}$: $\Gamma(s,c)\in W^s_{\rm loc}(p_c)$",(2.75,0),xytext=(0,4),textcoords="offset points",ha="right",color="#52514e",fontsize=8)
ax.set_xlim(-1.7,2.8)
ax.set_xlabel(r"$s-s_m\ \ (\times 10^{-9})$"); ax.set_ylabel(r"$g(s,c)\ \ (\times 10^{-7})$")
for s in ["top","right"]: ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig("paper/unfolding57.pdf"); fig.savefig("paper/figs/unfolding57.png",dpi=160)
print("ok")
