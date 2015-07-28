import matplotlib.pyplot as plt
import TLL as tll
import numpy as np

plt.rcParams['mathtext.fontset'] = 'stix'

u0 = 1.3
upi = .5
ls, gs, corr = tll.dvrg(u0, upi)
lmax = ls[-1]

g1 = np.array(gs[0][0])
g2 = np.array(gs[1][0])
g3 = np.array(gs[2][0])
g4 = np.array(gs[3][0])
cdws = g2 - g3 - 2.0*g1
cdwb = g2 + g3 - 2.0*g1
sdws = g2 + g3
sdwb = g2 - g3
ss = -g1 - g2
st = g1 - g2

phasegs = [(cdws,'cdws','red'),
           (cdwb,'cdwb','orange'),
           (sdws,'sdws','black'),
           (sdwb,'sdwb','green'),
           (ss,'ss','blue'),
           (st,'st','purple')]
fig = plt.figure()
ax1 = fig.add_subplot(211)
for g in gs:
    ax1.plot(ls, g[0], label=g[1], color=g[2], linewidth=2,
             alpha=.5, linestyle='--')
ax1.set_xlim(lmax*.9,lmax*1.025)
#ax1.set_ylim(-.1,.6)
ax1.set_title(r'$U_0=$'+str(u0)+r'$U_A\ \ \ U_{\pi}=$'+str(upi)+r'$U_A$',
              fontsize=16)
ax1.set_xlabel(r'RG Step $l$')
ax1.set_ylabel('Couplings')
ax1.legend(loc='best', fontsize=12)


ax2 = fig.add_subplot(212)
for cor in phasegs:
    cfunc = np.absolute(cor[0])
    ax2.plot(ls, cfunc, label=cor[1], color=cor[2],
             markersize=3, alpha=.5, marker='o')
ax2.set_ylabel(r'$poop$')
ax2.set_xlabel(r'RG Step $l$')
ax2.set_xlim(lmax*.999,lmax)
ax2.legend(loc='best', fontsize=12)

fname = str(u0)+'_'+str(upi)+'_corrs.pdf'
plt.savefig(fname)
plt.show()