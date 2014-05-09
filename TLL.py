# This program produces a phase diagram for a 1D-2D interacting fermi lattice
# model at half filling. The phase diagram is a 2-dimensional space, and each
# point in space corresponds to a unique set of values to start the RG-flow
# process. The RG-flow process can be done differentially until one of the
# couplings diverges. The phase space points correspond to different initail
# conditions to the ode solver.

import numpy as np
from scipy.integrate import odeint, ode

# dg/du = glgy(r,u) where g and glgy are 1 X 4 column vectors
def glgy(r,u):
    g1, g2, g3, g4 = r
    ufctr = 1. / (1. - u)**2
    gg1 = ufctr * -g1**2
    gg2 = ufctr * 0.5 * (-g1**2 + g3**2)
    gg3 = ufctr * g3 * (2.0*g2 - g1)
    gg4 = 0.0
    return np.array([gg1, gg2, gg3, gg4], float)

# dg/du = corrs where g and corrs are 1 X 10 vectors for solutions with
# correlations in regions where the couplings diverge
def corrs(u, r):
    g1, g2, g3, g4 = r[:4]
    ufctr = 1. / (1. - u)**2
    hg1, hg2, hg3, hg4 = glgy(r[:4],u)
    hcdws = ufctr * (g2 - g3 - 2.0*g1)
    hcdwb = ufctr * (g2 + g3 - 2.0*g1)
    hsdws = ufctr * (g2 + g3)
    hsdwb = ufctr * (g2 - g3)
    hss = ufctr * (-g1 - g2)
    hst = ufctr * (g1 - g2)
    return np.array([hg1,hg2,hg3,hg4,hcdws,hcdwb,hsdws,hsdwb,hss,hst], float)

# Compute ode for couplings in the tll region (no couplings diverge). Solve
# for various bosonized critical exponents and sort by their ode solution.
def tll(uzr, upi):
    r = np.array([1.-upi, 1.-uzr, 1.-upi, 1.-uzr], float)
    lpts = np.linspace(0., 1., 1000)
    gfuncs = odeint(glgy,r,lpts)
    g1, g2, g3, g4 = gfuncs[:,0], gfuncs[:,1], gfuncs[:,2], gfuncs[:,3]
    g2c, g2s = 0.5 * (3.0*g2 - g1), -0.5 * (g1 + g2)
    g4c, g4s = g4, -g4
    gc = np.sqrt((4.0*np.pi + g4c - g2c)/(4.0*np.pi + g4c + g2c))
    gs = np.sqrt((4.0*np.pi + g4s - g2s)/(4.0*np.pi + g4s + g2s))
    cdw = gc + gs
    cdw4 = 4.0 * gc
    sdwz = gc + gs
    sdwxy = gc + 1.0/gs
    sing = 1.0/gc + gs
    trip = 1.0/gc + 1.0/gs
    out = sorted([(cdw,'cdw','red','h'),
                  (cdw4,'cdw4','orange','h'),
                  (sdwz,'sdwz','yellow','h'),
                  (sdwxy,'sdwxy','green','h'),
                  (sing,'sing/trip0','blue','h'),
                  (trip,'trip+-','purple','h')],
                 key=lambda tup: np.absolute(tup[0][-1]))
    return out

# Compute ode for couplings in the non tll regions (some couplings diverge).
# Solve for the RG-flow of couplings and correlation functions, sorting the
# latter by ode solution.
def out(t,y):
    if np.isnan(y): return -1
    else: return None

def dvrg(uzr, upi):
    r0 = np.array([1.-upi, 1.-uzr, 1.-upi, 1.-uzr,
                  0., 0., 0., 0., 0., 0.], float)
    du = .001
    hfuncs = ode(corrs).set_integrator('dopri5')
    hfuncs.set_initial_value(r0)
    hfuncs.set_solout(out)
    uf = 1.
    ts = []
    sols = []

    while hfuncs.successful() and hfuncs.t < uf:
        ts.append(1. / (1-hfuncs.t))
        sols.append(hfuncs.y)        
        hfuncs.integrate(hfuncs.t+du)

    g1, g2, g3, g4 = [map(lambda x: x[i],sols) for i in range(4)]
    cdws, cdwb, sdws, sdwb, ss, st = [map(lambda x: x[i],sols) 
                                      for i in range(4,10)]
    coups = [(g1,'g1','red','^'),(g2,'g2','orange','^'),
             (g3,'g3','yellow','^'),(g4,'g4','green','^')]
    chis =  sorted([(cdws,'cdws','red','^'),
                    (cdwb,'cdwb','orange','^'),
                    (sdws,'sdws','yellow','^'),
                    (sdwb,'sdwb','green','^'),
                    (ss,'ss','blue','^'),
                    (st,'st','purple','^')], 
                   key=lambda p: np.absolute(p[0][-1]))
    return ts, coups, chis