from numpy import array, linspace, exp, append, copy, cos, pi, absolute

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
        return 4.0 * res / (200.0)**3

class Rgflows:

    def f(self, r, u):
        ufctr = 1.0 / (1.0 - u)**2
        f1 = ufctr * -(r[0]**2)
        f2 = ufctr * 0.5 * (r[2]**2 - r[0]**2)
        f3 = ufctr * r[2] * (2.0*r[1] - r[0])
        fcorr_cp = ufctr * (-2.0*r[0] + r[1] + r[2])
        fcorr_cm = ufctr * (-2.0*r[0] + r[1] - r[2])
        fcorr_sp = ufctr * (r[1] + r[2])
        fcorr_sm = ufctr * (r[1] - r[2])
        fcorr_ss = ufctr * (-r[0] - r[1])
        fcorr_st = ufctr * (r[0] - r[1])
        return array( [f1, f2, f3, fcorr_cp, fcorr_cm, fcorr_sp, fcorr_sm,
                         fcorr_ss, fcorr_st] )

    def flows(self, g):
        vals = copy(g)
        n = 1000
        h = 1.0 / n
        lpts = []
        rg_vals = array([ [] for i in range(9) ])
        for u in linspace(0.0,1.0,1000):
            lpts.append(u/(1.0-u))
            rg_vals = array([append(rg_vals[i], vals[i]) for i in range(9)])
            k1 = h * self.f(vals, u)
            k2 = h * self.f(vals + 0.5*k1, u + 0.5*h)
            k3 = h * self.f(vals + 0.5*k2, u + 0.5*h)
            k4 = h * self.f(vals + k3, u + h)
            vals += (k1 + 2.0*k2 + 2.0*k3 + k4) / 6.0

            #check g1, g2, g3 for divergent behavior and cut process
            gs = absolute(array(vals[:3]))
            old_gs = absolute(array([ a[-1] for a in rg_vals[:3] ]))
            if (u > 0.0 and (max(gs) > 60*max(old_gs))): return lpts, rg_vals
        return lpts, rg_vals
