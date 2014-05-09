import pylab as plt
import TLL as tll

ls, gs, corr = tll.dvrg(1.5,1.0)
plt.figure(1)
plt.subplot(211)
for g in gs:
    plt.plot(ls, g[0], label=g[1], color=g[2], marker=g[3])

plt.subplot(212)
for cor in cor:
    plt.plot(ls, cor[0], label=cor[1], color=cor[2], marker=cor[3])

plt.show()