# {Ueff} is an object contiaining information and methods for calculating the
# effective potential of a 1D-2D mixed Hubbard system with density-density
# interactions. [kpts] is the set of all momentum transfers, k4-k1, from a
# 2-particle scattering event and coinsides with the Broullin Zone. [mpts] is
# the range of 2D chemical potentials.

import numpy as np
from rg1d import *

func = Ueff()
kpts = np.linspace(-np.pi, np.pi, 40)
mpts = np.linspace(-4.0 , 4.0, 20)

for m in mpts:
    func.mu = m

    kvals = [ func(k) for k in kpts[-20:] ]
    kvals = kvals[::-1] + kvals #ueff is symmetric about the origin
    np.savetxt('kvals_'+str(m), kvals, delimiter=',') #store k4-k1 c

    xpts = range(40)
    xvals = np.fft.ifft(kvals) #inverse fourier transform to get ueff(x-x')
    np.savetxt('xvals_'+str(m), xvals.real, delimiter=',')

