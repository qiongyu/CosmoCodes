import math
import numpy as np
from spl import *

def sbesselj0(x):
  if x <= 0.01:
    sinc = 1.0 - x**2.0/6.0 + x**4.0/120.0 - x**6.0/5040.0
  else:
    sinc = math.sin(x)/x
  return sinc

def sbesselj1(x):
  if x<= 0.01:
    sphbesj1 = x/3.0 - x**3.0/30.0 + x**5.0/840.0
  else:
    sphbesj1 = (-x*math.cos(x) + math.sin(x)) / x**2.0
  return sphbesj1


def sigmaR(k, p, R):
#sigmaR: Calculates sigma_R from P(k)
# k = k values in h/Mpc, MUST be equally spaced in log
# power = P(k)
# R = R at which sigma_R is to be computed

  dlnk = math.log(k[1]/k[0])
  kR = k*R
  nkr = np.size(kR)
  m = range(0,nkr)
  coeff = 0.5 * k**3.0 * p * 9.0 / math.pi**2 * dlnk
  bes = np.zeros(nkr)
  for i in m:
    bes[i] = (sbesselj1(kR[i])/kR[i])**2.0
  cobes = coeff*bes
  sigmaR2 = sum(cobes[1:-1])
  sigmaR2 += 0.5*(cobes[0] + cobes[-1]) 

  return math.sqrt(sigmaR2)

def sigmaR_xi(r, xi, R):
#sigmaR_xi: Calculates sigma_R from xi(r)
# r = r values in Mpc/h
# xi = xi(r)
# R = R at which sigma_R is to be computed

  dr = 0.001
  rsp = np.arange(dr, 2.0*R+1.0, dr)
  xisp = spllin(rsp, xi, r)
  dex = np.squeeze(np.where(rsp <= 2.0*R))
  nr = np.size(dex)

  ror = rsp/R
  ror3 = ror*ror*ror
  r2 = rsp*rsp
  invR = 1.0/R

  integrand = (3.0-9.0/4.0*ror+3.0/16.0*ror3)*xisp*r2
  integrand = invR*invR*invR*integrand*dr
   
  sigmaR2 = sum(integrand[1:nr-1])
  sigmaR2 += 0.5*(integrand[0] + integrand[nr-1]) 
  
  return math.sqrt(sigmaR2)
