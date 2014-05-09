import pylab as plt
from numpy import arange
import TLL as tll

# Create a dictionary that maps tuples in phase space to tuples containing
# (label, color, marker) for plotting.
phase = {}
for uzr in arange(0., 1., 1./10.):
    for upi in arange(0., 1., 1./10.):
        ls, g, corr = tll.dvrg(uzr,upi)
        tup = corr[-1]
        phase[(uzr,upi)] = tup[1:]

for uzr in arange(0., 1., 1./10.):
    for upi in arange(1.+1./10., 2.+1./10., 1./10.):
        ls, g, corr = tll.dvrg(uzr,upi)
        tup = corr[-1]
        phase[(uzr,upi)] = tup[1:]
    
for uzr in arange(1., 2.+1./10., 1./10.):
    for upi in arange(1.+1./10., 2.+1./10., 1./10.):
        ls, g, corr = tll.dvrg(uzr,upi)
        tup = corr[-1]
        phase[(uzr,upi)] = tup[1:]

for uzr in arange(1.+1./10., 2.+1./10., 1./10.):
    for upi in arange(0., 1., 1./10.):
        xpnt = tll.tll(uzr,upi)
        tup = xpnt[0]
        phase[(uzr,upi)] = tup[1:]

# Partition the tuples into lists by (label, color, marker) for plotting in
# groups.
plots = [(filter(lambda x: phase[x]==val, phase), val)
         for val in list(set(phase.values()))]
print plots

#plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#plt.rc('text', usetex=True)

plt.figure()
plt.title(r'Phase Diagram of TLL Region')
plt.xlabel(r'$\frac{U_0}{U}$')
plt.ylabel(r'$\frac{U_{\pi}}{U}$', rotation=0)
plt.axis([0.0,2.0,0.0,2.0])
for tup in plots:
    plt.plot(*zip(*tup[0]), label=tup[1][0], color=tup[1][1],
             marker=tup[1][2], alpha=.45, linestyle='None',
             markersize=12, markeredgewidth=1.5)
#a = arange(0., 1.05, 1./20.)
#b = arange(0., 2.1, 1./20.)
#plt.plot(b, [1. for i in b], 'k--', markersize=8, alpha=.75, linewidth=5)
#plt.plot(b[-21:], b[-21:], 'k--', markersize=8, alpha=.75, linewidth=5)
#plt.plot([1. for i in a], a, 'k--', markersize=8, alpha=.75, linewidth=5)
plt.legend(loc=2, shadow=True, prop={'size':12})
plt.show()