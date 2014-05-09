import pylab as plt
import TLL as tll
import numpy as np

ls, gs, corr = tll.dvrg(1.5, 1.)
plt.figure(1)
plt.subplot(211)
for g in gs:
    plt.plot(ls, np.absolute(g[0]), label=g[1], color=g[2])
plt.legend()
plt.subplot(212)
for cor in corr:
    plt.semilogy(ls, np.absolute(cor[0]), label=cor[1], color=cor[2])
plt.legend()

xpnts = tll.tll(1., 1.)
plt.figure(2)
plt.xlim([5,20])
plt.ylim([0,.1])
for tup in xpnts:
    lst = tup[0]
    plt.plot(range(1,20), map(lambda x: (1./x)**(lst[-1]), range(1,20)), 
             label=tup[1], color=tup[2], marker=tup[3])
plt.legend()
plt.show()