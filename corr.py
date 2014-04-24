from rg1d import *
import pylab as plt
from numpy import *

class Soln:
    def __init__(self, soln, color, label):
        self.soln = soln
        self.color = color
        self.label = label

def sols(u_zr, u_pi):
    g1 = 1.0 - u_pi
    g2 = 1.0 - u_zr
    g3 = g1
    iconds = array( [g1, g2, g3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ] )
    phase = Rgflows()
    lpts,flows = phase.flows(iconds)
    flows = delete(absolute(flows), -1, 1)
    lpts = delete(lpts, -1)
    gs = [Soln(flows[0], 'red', '|g1|'),
          Soln(flows[1], 'green', '|g2|'),
          Soln(flows[2], 'blue', '|g3|'),
          Soln(2.0*flows[1]-flows[0], 'yellow', '|gc|')]
    corrs = [Soln(flows[3], 'cyan', 'cdw+'),
             Soln(flows[4], 'magenta', 'cdw-'),
             Soln(flows[5], 'yellow', 'sdw+'),
             Soln(flows[6], 'black', 'sdw-'),
             Soln(flows[7], 'orange', 'ss'),
             Soln(flows[8], 'pink', 'st')]
    cor_sort = sorted(corrs, key=lambda obj: obj.soln[-1])
    return gs, cor_sort

uhash = {}
for uzr in linspace(0.0,2.0,40):
    for upi in linspace(0.0,2.0,40):
        gs, cor_sort = sols(uzr,upi)
#        plt.figure(1)
#        plt.subplot(211)
#        plt.title('|g| vs. RG-step')
#        for g in gs:
#            plt.plot(lpts, g.soln, color=g.color, label=g.label)
#            plt.legend(loc=1)
#
#        plt.subplot(212)
#        plt.title('ln(Correlation) vs. RG-step')
#        for cor in cor_sort[-3:]:
#            plt.semilogy(lpts, cor.soln, color=cor.color, label=cor.label)
#            ymax = max(cor.soln[-1] for cor in cor_sort[-3:])
#            plt.ylim(ymax*2.0, ymax/2.0)
#            plt.legend(loc=2)
#
#        plt.clf()
        toplot = cor_sort[-1]
        uhash[(uzr,upi)] = toplot.label

cdwp = [ key for key in uhash if uhash[key]=='cdw+' ]
cdwm = [ key for key in uhash if uhash[key]=='cdw-' ]
sdwp = [ key for key in uhash if uhash[key]=='sdw+' ]
sdwm = [ key for key in uhash if uhash[key]=='sdw-' ]
ss = [ key for key in uhash if uhash[key]=='ss' ]
st = [ key for key in uhash if uhash[key]=='st' ]

plt.plot(*zip(*cdwp), color='red', marker='o', linestyle='None', alpha=0.5,
         label='cdw+')
plt.plot(*zip(*cdwm), color='orange', marker='s', linestyle='None', alpha=0.5,
         label='cdw-')
plt.plot(*zip(*sdwp), color='yellow', marker='p', linestyle='None', alpha=0.5,
         label='sdw+')
plt.plot(*zip(*sdwm), color='green', marker='D', linestyle='None', alpha=0.5,
         label='sdw-')
plt.plot(*zip(*ss), color='blue', marker='^', linestyle='None', alpha=0.5,
         label='ss')
plt.plot(*zip(*st), color='purple', marker='h', linestyle='None', alpha=0.5,
         label='st')
plt.xlabel('Uo/U')
plt.ylabel('Upi/U')
plt.title('Phase Diagram')
plt.legend()
plt.show()
