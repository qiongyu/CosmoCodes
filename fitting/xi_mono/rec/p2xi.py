import math
import numpy as np
from spline import *
from scipy import interpolate

def p2xi(k, p, kr, J0, r):
#p2xi: Calculates xi(r) from P(k)
# k = k values in h/Mpc, range MUST be > range of kr/r_min
# p = P(k)
# kr = kr values where J0 is defined, MUST be equally spaced in log
# J0 = Spherical Bessel function of order 0
# r = r values at which to evaluate xi

# set up spline onto the filter, THIS IS ONLY USEABLE if power > 0
# use p2xi_lin if there are values of power < 0!!
  lk = np.log10(k)
  lp = np.log10(p)

  spform = interpolate.splrep(lk, lp, s=0)

  dlnk = math.log(kr[1]/kr[0])
  nr = np.size(r)
  m = range(0, nr)
  xi = np.zeros(nr)

  for i in m:
    kk = kr/r[i]
    lkk = np.log10(kk)
    pp = 10.0**(interpolate.splev(lkk,spform,der=0))
    
# small scale (large k) damping using exp => better numerics
    xi[i] = sum(pp * kk**3 * J0 * np.exp(-kk**2))

  xi = xi*dlnk/(2.0*math.pi**2)
 
  return xi

def p2xi_lin(k, p, kr, J0, r):
#p2xi: Calculates xi(r) from P(k)
# k = k values in h/Mpc, range MUST be > range of kr/r_min
# p = P(k)
# kr = kr values where J0 is defined, MUST be equally spaced in log
# J0 = Spherical Bessel function of order 0
# r = r values at which to evaluate xi

  lk = k
  lp = p

  spform = interpolate.splrep(lk, lp, s=0)

  dlnk = math.log(kr[1]/kr[0])
  nr = np.size(r)
  m = range(0, nr)
  xi = np.zeros(nr)

  for i in m:
    kk = kr/r[i]
    lkk = kk
    pp = interpolate.splev(lkk,spform,der=0)
    
# small scale (large k) damping using exp => better numerics
    xi[i] = sum(pp * kk**3 * J0 * np.exp(-kk**2))

  xi = xi*dlnk/(2.0*math.pi**2)
 
  return xi
