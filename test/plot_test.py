import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

testpoints = []
rows = ['r1','r2','r3','r4','r5','r6']
cols = ['c1','c2','c3','c4','c5','c6']
pags = ['p1','p2','p3','p4','p5','p6']
ords = ['o3','o5','o7','o9']
legendlines = []
ocolors = ['r', 'g', 'b', 'k']
for c in ocolors:
    legendlines.append(Line2D([0],[0], color=c,lw=4))
eps = [10**float(i) for i in range(4,-16,-1)]
errs = []
fig,ax = plt.subplots(constrained_layout=True)
for i, r in enumerate(rows):
    errs.append([])
    for j, c in enumerate(cols):
        errs[-1].append([])
        for l, o in enumerate(ords):
            errs[-1][-1].append([])
            with open(''.join(["./outs/err",r,c,o,".txt"]),'r') as f:
                for line in f:
                    errs[-1][-1][-1].append(float(line))
            ax.plot(eps,errs[-1][-1][-1],label=''.join([r,c,o]),alpha=0.5,color=ocolors[l])

errs2 = []
fig2,ax2 = plt.subplots(constrained_layout=True)
for i, r in enumerate(rows):
    errs2.append([])
    for j, c in enumerate(cols):
        errs2[-1].append([])
        for k, p in enumerate(pags):
            errs2[-1][-1].append([])
            for l, o in enumerate(ords):
                errs2[-1][-1][-1].append([])
                with open(''.join(["./outs/err2",r,c,p,o,".txt"]),'r') as f:
                    for line in f:
                        errs2[-1][-1][-1][-1].append(float(line))
                ax2.plot(eps,errs2[-1][-1][-1][-1],label=''.join([r,c,o,p]),alpha=0.5,color=ocolors[l])
ax.set_yscale('log')
ax2.set_yscale('log')
ax.set_xscale('log')
ax2.set_xscale('log')
ax.invert_xaxis()
ax2.invert_xaxis()
ax.grid()
ax2.grid()
ax.set_title("Jacobian error")
ax2.set_title("Hessian error")
ax.set_xlabel("$\epsilon$")
ax2.set_xlabel("$\epsilon$")
ax.set_ylabel("error")
ax2.set_ylabel("error")
ax.set_xticks(eps[::2])
ax2.set_xticks(eps[::2])
legendlabels = ["3rd order", "5th order", "7th order", "9th order"]
ax.legend(legendlines, legendlabels)
ax2.legend(legendlines, legendlabels)
plt.show()
