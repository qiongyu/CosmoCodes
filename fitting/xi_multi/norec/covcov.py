import math
import numpy as np
from scipy import interpolate
from scipy import array
from scipy import special
from matplotlib import pyplot
from TwoDP import *
from scipy import integrate

def sbesselj1(x):
  nx = np.size(x)
  j1 = np.zeros(nx)
  for i in range(0,nx):
    if math.fabs(x[i]) <= 0.01:
      j1[i] = x[i]/3.0 - x[i]**3/30.0 + x[i]**5/840.0 - x[i]**7/45360.0 
    else:
      j1[i] = (math.sin(x[i])/x[i] - math.cos(x[i]))/x[i]
  return j1

def sbesselj0(x):
  nx = np.size(x)
  sinc = np.zeros(nx)
  for i in range(0,nx):
    if math.fabs(x[i]) <= 0.01:
      sinc[i] = 1.0 - x[i]**2.0/6.0 + x[i]**4.0/120.0 - x[i]**6.0/5040.0
    else:
      sinc[i] = math.sin(x[i])/x[i]
  return sinc

def SI(x):
  nx = np.size(x)
  sint = np.zeros(nx)
  for i in range(0,nx):
    if math.fabs(x[i]) < 0.01:
      sint[i] = x[i] - x[i]**2/18.0 + x[i]**5/600.0 - x[i]**7/35280.0
    else:
      sint[i], junk = special.sici(x[i])
  return sint

def ker02(r,k):
  kr = k*r
  kern0 = r**2/k*sbesselj1(kr) 
  kern2_c1 = 3.0*SI(kr)/k**3
  kern2_c2 = 3*r/k*sbesselj0(kr)
  kern2_c3 = r**2*sbesselj1(kr)
  kern2 = kern2_c1 - (kern2_c2+kern2_c3)/k
  return kern0, kern2

def volint(x,weight,params,zdist,pdw,k,spar,nbar,mesh):
#This function calculates the inverse variance as an integral over volume.
#The variance is equal to:
#[(c0*P(k) + c1/n)*stream*kaiser + c2/n]**2 (from the covariance matrix)
#       params: [c0,c1,c2]
#       pow: P(k)
#       nk: # of elements in P(k)
#       stream: streaming*kaiser term (see stream.py)
#       nbar: n(z)
#       mesh: const*r(z)**2*dz/sqrt(omega_m(1+z)**3+omega_lambda), 
#             const=omega*c/H_0, omega == solid angle in steradians
  par0 = params[0]
  par1 = params[1]
  par2 = params[2]

  nk = np.size(k)
  inttot00 = np.zeros(nk)
  inttot02 = np.zeros(nk)
  inttot22 = np.zeros(nk)

  nn = np.size(nbar)
  pell00 = np.zeros(nn)
  pell02 = np.zeros(nn)
  pell22 = np.zeros(nn)

  for cc in range(0,nk):
    for dd in range(0,nn):
      f00, f02, f22 = TwoDPn(x, pdw[:,cc], zdist[:,cc], params, nbar[dd])
      pell00[dd] = 0.5*np.sum(f00*weight)
      pell02[dd] = 0.5*np.sum(f02*weight)
      pell22[dd] = 0.5*np.sum(f22*weight)
      
    intexp = mesh/pell00
    inttot00[cc] = 1.0/np.sum(intexp)

    intexp = mesh/pell02
    inttot02[cc] = 1.0/np.sum(intexp)

    intexp = mesh/pell22
    inttot22[cc] = 1.0/np.sum(intexp)
  return inttot00, inttot02, inttot22

def cov02(k, pv00, pv22, pv02, rs):
# cov02: Calculates the full covariance matrix of monopole + quadrupole
# 	 This is for simultaneously fitting the monopole and the quadrupole
# k = k values in h/Mpc, MUST be equally spaced in log
# power = P(k)
# power = P_2(k), quadrupole spectrum
# rs = rs values at which to evaluate the covariance matrix
# volume = volume of simulation in Mpc^3/h^3
# noise = shot-noise (added equally to each mode of P(k))
# stream = streaming term

# C = 2(2l+1)(2l'+1)/V integral{ dlnk k^3/(2pi^2) W_l(krs)W_l'(krs) I_ll' }
# 2(2l+1)(2l'+1) = 2 for l=l'=0
#		   50 for l=l'=2
#		   10 for l=0, l'=2
# I_ll' = 1/2 integral(-1,1){ dmu L_l(mu)L_l'(mu) sum[P_L(k)L_L(mu) + n^-1]^2 }
# If we expand out the sum to L=2...we get some Legendre integrals.

  dlnk = math.log(k[1]/k[0])
# coeffll' = I_ll' * k^3 * 2(2l+1)(2l'+1) / 2 (<-- in 2pi**2)
  coeff00 = pv00 * k**3 * 9.0
  coeff22 = pv22 * k**3 * 225.0
  coeff02 = pv02 * k**3 * 45.0
  nrs = np.size(rs)
  nk = np.size(k)

  drs = rs[1]-rs[0] 
  covar = np.zeros((2*nrs,2*nrs))
  bins = np.zeros(nrs+1)
  bins[1:] = rs+0.5*drs
  bins[0] = bins[1]-drs

  W0 = np.zeros((nk,nrs))
  W2 = np.zeros((nk,nrs))

  m = range(0,nrs)
  for i in m:
    W0a, W2a = ker02(bins[i],k)
    W0b, W2b = ker02(bins[i+1],k)
    W0[:,i] = W0b - W0a
    W2[:,i] = W2b - W2a

# The 00 covariance matrix...(diagonal terms)
    covar[i,i] = sum(coeff00[1:-1]*W0[1:-1,i]**2)
    covar[i,i] += 0.5*(coeff00[0]*W0[0,i]**2 + coeff00[-1]*W0[-1,i]**2)
    covar[i,i] *= 1.0/(bins[i+1]**3-bins[i]**3)**2

# The 22 covariance matrix...(diagonal terms)
    covar[i+nrs,i+nrs] = sum(coeff22[1:-1]*W2[1:-1,i]**2)
    covar[i+nrs,i+nrs] += 0.5*(coeff22[0]*W2[0,i]**2 + coeff22[-1]*W2[-1,i]**2)
    covar[i+nrs,i+nrs] *= 1.0/(bins[i+1]**3-bins[i]**3)**2

    j = 0
    while j < i:
# The 00 covariance matrix...(off diagonal terms)
      covar[i,j] = sum(coeff00[1:-1]*W0[1:-1,i]*W0[1:-1,j])
      covar[i,j] += 0.5*(coeff00[0]*W0[0,i]*W0[0,j] + coeff00[-1]*W0[-1,i]*W0[-1,j])
      covar[i,j] *= 1.0/(bins[i+1]**3-bins[i]**3)/(bins[j+1]**3-bins[j]**3)
      covar[j,i] = covar[i,j]

# The 22 covariance matrix...(off diagonal terms)
      covar[i+nrs,j+nrs] = sum(coeff22[1:-1]*W2[1:-1,i]*W2[1:-1,j])
      covar[i+nrs,j+nrs] += 0.5*(coeff22[0]*W2[0,i]*W2[0,j] + coeff22[-1]*W2[-1,i]*W2[-1,j])
      covar[i+nrs,j+nrs] *= 1.0/(bins[i+1]**3-bins[i]**3)/(bins[j+1]**3-bins[j]**3)
      covar[j+nrs,i+nrs] = covar[i+nrs,j+nrs]
      j += 1

# The cross terms...between mono and quadrupoles
# The 02 covariance matrix is the transpose of the 20 matrix
  for i in m:
    for j in m:
      covar[i+nrs,j] = sum(coeff02[1:-1]*W2[1:-1,i]*W0[1:-1,j])
      covar[i+nrs,j] += 0.5*(coeff02[0]*W2[0,i]*W0[0,j] + coeff02[-1]*W2[-1,i]*W0[-1,j])
      covar[i+nrs,j] *= -1.0/(bins[i+1]**3-bins[i]**3)/(bins[j+1]**3-bins[j]**3)
      covar[j,i+nrs] = covar[i+nrs,j]

  covar = covar/(math.pi**2)*dlnk
  return covar

def redcov(cov):
#redcov: calculates the reduced covariance matrix
#cov = covariance matrix

  dim = int(np.sqrt(np.size(cov)))
  redcov = np.zeros((dim,dim))
  for i in range(0,dim):
    for j in range(0,dim):
      redcov[i,j] = cov[i,j]/math.sqrt(cov[i,i]*cov[j,j])

  return redcov

#xin = np.arange(-10.0,10.0,0.1)
#j0 = sbesselj0(xin)
#j1 = sbesselj1(xin)
#sintegral = SI(xin)

#pyplot.plot(xin,j0)
#pyplot.plot(xin,j1)
#pyplot.plot(xin,sintegral)
#pyplot.savefig('moomoo.ps')
#pyplot.close()
