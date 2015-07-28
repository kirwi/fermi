from rg1d import *
import pylab as plt
from numpy import *

class Soln:
    def __init__(self, soln, label, color):
        self.soln = soln
        self.label = label
        self.color = color

def sols(u_zr, u_pi):
    g1 = 1.0 - u_pi
    g2 = 1.0 - u_zr
    g3 = 1.0 - u_pi
    iconds = array([g1, g2, g3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ], float)
    phase = Rgflows()
    lpts,flows = phase.flows(iconds)
    flows = delete(absolute(flows), -1, 1)
    lpts = delete(lpts, -1)
    gs = array([Soln(flows[0], '|g1|', 'red'),
          Soln(flows[1], '|g2|', 'orange'),
          Soln(flows[2], '|g3|', 'blue'),
          Soln(2.0*flows[1]-flows[0], '|gc|')])
    corrs = [Soln(flows[3], 'cdw+'),
             Soln(flows[4], 'cdw-'),
             Soln(flows[5], 'sdw+'),
             Soln(flows[6], 'sdw-'),
             Soln(flows[7], 'ss'),
             Soln(flows[8], 'st')]
    cor_sort = array(sorted(corrs, key=lambda obj: obj.soln[-1]))
    return gs, cor_sort

uhash = {}
for uzr in linspace(0.0,2.0,40):
    for upi in linspace(0.0,2.0,40):
        gs, cor_sort = sols(uzr,upi)
        toplot = cor_sort[-1]
        uhash[(uzr,upi)] = toplot.label

