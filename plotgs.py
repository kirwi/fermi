import matplotlib.pyplot as plt
import TLL as tll
import numpy as np

plt.rcParams['mathtext.fontset'] = 'stix'

u0 = 1.5
upi = .5
ls, gs, corr = tll.dvrg(u0, upi)
lmax = ls[-1]

fig = plt.figure()
ax1 = fig.add_subplot(111)
mcount = 1
for g in gs:
    skip = mcount*2
    ax1.plot(ls, g[0], color=g[2], alpha=.5, linestyle='--')
    ax1.scatter(ls[::skip], g[0][::skip], label=g[1], color=g[2], marker='o',
                alpha=.5)
    mcount += 1
ax1.set_xlim(0,200)
ax1.set_title(r'$U_0=$'+str(u0)+r'$U_A\ \ \ U_{\pi}=$'+str(upi)+r'$U_A$',
              fontsize=16)
ax1.set_xlabel(r'RG Step $l$')
ax1.set_ylabel('g-ology')
ax1.legend(loc='best', fontsize=12)

#ax2 = fig.add_subplot(212)
#mcount = 1
#for cor in corr:
#    skip = mcount * 5
#    cfunc = np.absolute(cor[0])
#    ax2.plot(ls, cfunc, color=cor[2], linestyle='--', alpha=.5)
#    ax2.scatter(ls[::skip], cfunc[::skip], color=cor[2], label=cor[1],
#             marker='o', alpha=.5)
#    mcount += 1
#ax2.set_xlim(0, 1)
#ax2.set_ylabel(r'$|g_{\delta}|$')
#ax2.set_xlabel(r'RG Step $l$')
#ax2.legend(loc='best', fontsize=12)

plt.savefig('stss.pdf')
plt.show()