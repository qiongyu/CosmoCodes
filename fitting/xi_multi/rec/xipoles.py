import math
import numpy as np
from spline import *
from scipy import interpolate
from scipy import integrate
from matplotlib import pyplot

def P_ell(mu, weight, p2d, k, order):
# Calculates the multipoles of P(k,mu), i.e. P_\ell(k), where \ell = order
  nmu = np.size(mu)
  L = np.zeros(nmu)
  coeff = (2.*order+1)/2.
  if(order==0):
    L = L+1.
  if(order==2):
    L = L+0.5*(3.*mu**2-1.)
  if(order==4):
    L = L+0.125*(35.*mu**4-30.*mu**2+3.)
  L = L*coeff

  nk = np.size(k)
  pell = np.zeros(nk)
  for i in range(0,nk):
    pell[i] = np.sum(p2d[:,i]*L*weight)

  #pyplot.plot(k, k**3*pell, '-k')
  #pyplot.xlim(0.0,0.4)
  #pyplot.ylim(0.0,100.0)
  #pyplot.savefig('pole'+str(order)+'.ps')
  #pyplot.clf()
  return pell

def xipole(k, p, kr, Jn, r, order):
#xipole: Calculates xi_\ell from P_\ell using log splining
# k = k values in h/Mpc, range MUST be > range of kr/r_min
# p = P(k)
# kr = kr values where Jn is defined, MUST be equally spaced in log
# Jn = Spherical Bessel function of order 0
# r = r values at which to evaluate xi

# set up spline onto the filter, THIS IS ONLY USEABLE if power > 0
# use xipole_lin if there are values of power < 0!!
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
    xi[i] = sum(pp * kk**3 * Jn * np.exp(-kk**2))

  if(order==0 or order==4):
    ipow = 1.
  if(order==2):
    ipow = -1.

  xi = ipow*xi*dlnk/(2.0*math.pi**2)
 
  return xi

def xipole_lin(k, p, kr, Jn, r, order):
#xipole_lin: Calculates xi_\ell from P_\ell using linear splining
# k = k values in h/Mpc, range MUST be > range of kr/r_min
# p = P(k)
# kr = kr values where Jn is defined, MUST be equally spaced in log
# Jn = Spherical Bessel function of order 0
# r = r values at which to evaluate xi

  spform = interpolate.splrep(k, p, s=0)
#  jform = interpolate.splrep(kr, Jn, s=0)

  dlnk = math.log(kr[1]/kr[0])
  nr = np.size(r)
  m = range(0, nr)
  xi = np.zeros(nr)

  for i in m:
    kk = kr/r[i]
    pp = interpolate.splev(kk,spform,der=0)
    
# small scale (large k) damping using exp => better numerics
    xi[i] = sum(pp * kk**3 * Jn * np.exp(-kk**2))

  if(order==0 or order==4):
    ipow = 1.
  if(order==2):
    ipow = -1.

  xi = ipow*xi*dlnk/(2.0*math.pi**2)
 
  return xi

def xi_ell(x, weight, pm, ksp, Jn, kr, r, order):
  pmell = P_ell(x,weight,pm,ksp,order)
  xiell = xipole_lin(ksp,pmell,kr,Jn,r,order)

  return xiell

def P_ell_sep(mu, weight, p2d, k, order):
# Calculates the multipoles of P(k,mu), i.e. P_\ell(k), where \ell = order
  nmu = np.size(mu)
  L = np.zeros(nmu)
  coeff = (2.*order+1)/2.
  if(order==0):
    L = L+1.
  if(order==2):
    L = L+0.5*(3.*mu**2-1.)
  if(order==4):
    L = L+0.125*(35.*mu**4-30.*mu**2+3.)
  L = L*coeff

  nk = np.size(k)
  pell = np.zeros((nk,3))
  for i in range(0,nk):
    pell[i,0] = np.sum(p2d[:,i,0]*L*weight)
    pell[i,1] = np.sum(p2d[:,i,1]*L*weight)
    pell[i,2] = np.sum(p2d[:,i,2]*L*weight)

  #pyplot.plot(k, k**3*pell, '-k')
  #pyplot.xlim(0.0,0.4)
  #pyplot.ylim(0.0,100.0)
  #pyplot.savefig('pole'+str(order)+'.ps')
  #pyplot.clf()
  return pell

def xi_ell_sep(x, weight, pm, ksp, Jn, kr, r, order):
  pmell = P_ell_sep(x,weight,pm,ksp,order)
  xiell1 = xipole_lin(ksp,pmell[:,0],kr,Jn,r,order)
  xiell2 = xipole_lin(ksp,pmell[:,1],kr,Jn,r,order)
  xiell3 = xipole_lin(ksp,pmell[:,2],kr,Jn,r,order)

  return xiell1, xiell2, xiell3
