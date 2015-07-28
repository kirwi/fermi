from numpy import *

class Ueff:

    def __init__(self, w=0.001):
        self.w = w
        self.ehash = {}
        self.mu = None

    def energy(self, x, y):
        if (x,y) in self.ehash:
            return self.ehash[(x,y)]
        else:
            self.ehash[(x,y)] = -2.0 * (cos(x) + cos(y))
            return self.ehash[(x,y)]

    def shape(self, x, y, z, k):
        if ((self.energy(x,y) > self.mu) and (self.energy(x+k,z) < self.mu)):
            return -1.0
        if ((self.energy(x,y) < self.mu) and (self.energy(x+k,z) > self.mu)):
            return 1.0
        return 0.0

    def __call__(self, k):
        res = 0.0
        bzone = linspace(-pi, pi, 200)
        for x in bzone:
            for y in bzone[-100:]:
                for z in bzone[-100:]:
                    res += (self.shape(x,y,z,k) *
                            (self.energy(x,y) - self.energy(x+k,z)) /
                            ((self.energy(x,y) - self.energy(x+k,z))**2 +
                            self.w**2))
        return 4.0 * res / 200.0**3

class Rgflows:

    # Returns an array containing inhomogeneous part of the RG equations to
    # solve. dg/dl = f(g,l) where both f and g are numpy arrays.
    def f(self, u, r):
        g1, g2, g3 = r[:3]
        f1 = -g1**2
        f2 = 0.5 * (g3**2 - g1**2)
        f3 = g3 * (2.0*g2 - g1)
        fcorr_cp = -2.0*g1 + g2 + g3
        fcorr_cm = -2.0*g1 + g2 - g3
        fcorr_sp = g2 + g3
        fcorr_sm = g2 - g3
        fcorr_ss = -g1 - g2
        fcorr_st = g1 - g2
        return array([f1, f2, f3, fcorr_cp, fcorr_cm, fcorr_sp, fcorr_sm,
                       fcorr_ss, fcorr_st], float)

    # Performs Runge-Kutta solution across the input array r and returns a
    # tuple (lpts, out), numpy arrays containing the RG points and updated
    # solutions respectively. The Runge-Kutta process stops when any of the
    # couplings has "diverged" as determined by the if statement.
    def flows(self, r):
        lpts = []
        out = array([ [a] for a in r ])
        for t in range(0, 1000):
            lpts.append(t)
            k1 = h * self.f(t, r)
            k2 = h * self.f(t + 0.5*h, r + 0.5*k1)
            k3 = h * self.f(t + 0.5*h, r + 0.5*k2)
            k4 = h * self.f(t + h, r + k3)
            r += (k1 + 2.0*k2 + 2.0*k3 + k4) / 6.0
            if isnan(amin(r)): return lpts, out
            out = insert(out, out.shape[-1], r, axis=1)
        return lpts, out
