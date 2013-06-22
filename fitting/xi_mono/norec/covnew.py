import math
import numpy as np
from scipy import interpolate
from scipy import array
from matplotlib import pyplot

def sbesselj1(x):
  nx = np.size(x)
  j1 = np.zeros(nx)
  for i in range(0,nx):
    if x[i] <= 0.01:
      j1[i] = x[i]/3.0 - x[i]**3/30.0 + x[i]**5/840.0 - x[i]**7/45360.0 
    else:
      j1[i] = (math.sin(x[i])/x[i] - math.cos(x[i]))/x[i]
  return j1

def volint(params,pow,nk,stream,nbar,mesh):
#This function calculates the inverse variance as an integral over volume.
#The variance is equal to:
#[(c0*P(k) + c1/n)*stream*kaiser + c2/n]**2 (from the covariance matrix)
#	params: [c0,c1,c2]
#	pow: P(k)
#	nk: # of elements in P(k)
#	stream: streaming*kaiser term (see stream.py)
#	nbar: n(z)
#	mesh: const*r(z)**2*dz/sqrt(omega_m(1+z)**3+omega_lambda), 
#	      const=omega*c/H_0, omega == solid angle in steradians
  par0 = params[0]
  par1 = params[1]
  par2 = params[2]
  inttot = np.zeros(nk)
  for cc in range(0,nk):
    intexp = 1.0/((par0*pow[cc] + par1/nbar)*stream[cc] + par2/nbar)**2.
    intexp *= mesh
    inttot[cc] = 1.0/np.sum(intexp)
  return inttot

def cov00(k, power2, krs, rs):
# cov00: Calculates the covariance matrix of the monopole omega_0, assumes
# P(k) = P_0(k) is a good approximation (i.e. no higher order terms)
# 	k = k values in h/Mpc, MUST be equally spaced in log
# 	power2 = (P(k)+1/n)**2 term...this may or may not include streaming, 
#		Kaiser, whatever. This assumes that we have done an integral 
#	    	over the volume, see the function volint, to obtain power2.
# 	krs = krs values where W0 is defined
# 	rs = rs values at which to evaluate the covariance matrix

  dlnk = math.log(k[1]/k[0])
  coeff = power2*k
  nrs = np.size(rs)
  nk = np.size(k)
  drs = rs[1]-rs[0]
 
  covar = np.zeros((nrs,nrs))
  bins = np.zeros(nrs+1)
  bins[1:] = rs+0.5*drs
  bins[0] = bins[1]-drs

  m = range(0,nrs)
  
  W0a1 = np.zeros((nk,nrs))
  W0b1 = np.zeros((nk,nrs))
 
  for i in m: 
    krsa1 = k*bins[i]
    krsb1 = k*bins[i+1]
    W0a1[:,i] = bins[i]**2*sbesselj1(krsa1)
    W0b1[:,i] = bins[i+1]**2*sbesselj1(krsb1)
    covar[i,i] = sum(coeff[1:-1]*(W0b1[1:-1,i]-W0a1[1:-1,i])**2)
    covar[i,i] += 0.5*(coeff[0]*(W0b1[0,i]-W0a1[0,i])**2 + coeff[-1]*(W0b1[-1,i]-W0a1[-1,i])**2)
    covar[i,i] *= 9.0/(bins[i+1]**3-bins[i]**3)**2

    j = 0
    while j < i: #I can nest the while because j < i always =)
      covar[i,j] = sum(coeff[1:-1]*(W0b1[1:-1,i]-W0a1[1:-1,i])*(W0b1[1:-1,j]-W0a1[1:-1,j]))
      covar[i,j] += 0.5*(coeff[0]*(W0b1[0,i]-W0a1[0,i])*(W0b1[0,j]-W0a1[0,j]) + coeff[-1]*(W0b1[-1,i]-W0a1[-1,i])*(W0b1[-1,j]-W0a1[-1,j]))
      covar[i,j] *= 9.0/(bins[i+1]**3-bins[i]**3)/(bins[j+1]**3-bins[j]**3)
      covar[j,i] = covar[i,j]
      j += 1
  
  covar = covar/(math.pi**2)*dlnk
  return covar

def cov02(k, power, power2, krs, W0, W2, rs, volume, noise):
# cov02: Calculates the full covariance matrix of monopole + quadrupole
# 	 This is for simultaneously fitting the monopole and the quadrupole
# k = k values in h/Mpc, MUST be equally spaced in log
# power = P(k)
# power = P_2(k), quadrupole spectrum
# krs = krs values where W0 is defined
# W0 = W0(krs) 
# W2 = W2(krs)
# rs = rs values at which to evaluate the covariance matrix
# volume = volume of simulation in Mpc^3/h^3
# noise = shot-noise (added equally to each mode of P(k))

# C = 2(2l+1)(2l'+1)/V integral{ dlnk k^3/(2pi^2) W_l(krs)W_l'(krs) I_ll' }
# 2(2l+1)(2l'+1) = 2 for l=l'=0
#		   50 for l=l'=2
#		   10 for l=0, l'=2
# I_ll' = 1/2 integral(-1,1){ dmu L_l(mu)L_l'(mu) sum[P_L(k)L_L(mu) + n^-1]^2 }
# If we expand out the sum to L=2...we get some Legendre integrals.

  I00 = power**2 + 0.2*power2**2 + 2.0*power*noise + noise**2
  I22 = 0.2*power**2 + 3.0*power2**2/35.0 + 4.0*power*power2/35.0 + 0.4*power*noise + 4.0*power2*noise/35.0 + 0.2*noise**2
  I02 = 2.0*power2**2/35.0 + 0.4*power*power2 + 0.4*power2*noise

  dlnk = math.log(k[1]/k[0])
# coeffll' = I_ll' * k^3 * 2(2l+1)(2l'+1) / 2 (<-- in 2pi**2)
  coeff00 = I00 * k**3 * 1.0
  coeff22 = I22 * k**3 * 25.0
  coeff02 = I02 * k**3 * 5.0
  nrs = np.size(rs)
  nk = np.size(k)
 
  spform0 = interpolate.splrep(krs,W0,s=0)
  spform2 = interpolate.splrep(krs,W2,s=0)
  covar = np.zeros((2*nrs,2*nrs))
  W0sp = np.zeros((nk,nrs))
  W2sp = np.zeros((nk,nrs))

  m = range(0,nrs)
  for i in m:
    krs1 = k*rs[i]
    temp0 = interpolate.splev(krs1,spform0,der=0)
    temp2 = interpolate.splev(krs1,spform2,der=0)
    W0sp[:,i] = np.where(krs1 < 1000.0, temp0, 0.0)
    W2sp[:,i] = np.where(krs1 < 1000.0, temp2, 0.0)
# The 00 covariance matrix...(diagonal terms)
    covar[2*i,2*i] = sum(coeff00[1:-1]*W0sp[1:-1,i]**2)
    covar[2*i,2*i] += 0.5*(coeff00[0]*W0sp[0,i]**2 + coeff00[-1]*W0sp[-1,i]**2)
# The 22 covariance matrix...(diagonal terms)
    covar[2*i+1,2*i+1] = sum(coeff22[1:-1]*W2sp[1:-1,i]**2)
    covar[2*i+1,2*i+1] += 0.5*(coeff22[0]*W2sp[0,i]**2 + coeff22[-1]*W2sp[-1,i]**2)

    j = 0
    while j < i:
# The 00 covariance matrix...(off diagonal terms)
      covar[2*i,2*j] = sum(coeff00[1:-1]*W0sp[1:-1,i]*W0sp[1:-1,j])
      covar[2*i,2*j] += 0.5*(coeff00[0]*W0sp[0,i]*W0sp[0,j] + coeff00[-1]*W0sp[-1,i]*W0sp[-1,j])
      covar[2*j,2*i] = covar[2*i,2*j]
# The 22 covariance matrix...(off diagonal terms)
      covar[2*i+1,2*j+1] = sum(coeff22[1:-1]*W2sp[1:-1,i]*W2sp[1:-1,j])
      covar[2*i+1,2*j+1] += 0.5*(coeff22[0]*W2sp[0,i]*W2sp[0,j] + coeff22[-1]*W2sp[-1,i]*W2sp[-1,j])
      covar[2*j+1,2*i+1] = covar[2*i+1,2*j+1]
      j += 1

# The cross terms...(the blanks below)
# If cov00 = ( x x x x ... ) and cov22 = ( o o o o ... ), then
# cov02 = ( x   x   x   x   ...  )
#	  (   o   o   o   o  ... ) etc.
# The 02 covariance matrix is the transpose of the 20 matrix
  c02 = np.zeros((nrs,nrs))
  for i in m:
    for j in m:
      c02[i,j] = sum(coeff02[1:-1]*W0sp[1:-1,i]*W2sp[1:-1,j])
      c02[i,j] += 0.5*(coeff02[0]*W0sp[0,i]*W2sp[0,j] + coeff02[-1]*W0sp[-1,i]*W2sp[-1,j])
  for i in m:
    for j in m:
      covar[i*2,j*2+1] = c02[i,j]
      covar[j*2+1,i*2] = c02[i,j]

  covar = covar/(math.pi**2)/volume*dlnk
  return covar

def redcov(cov, sub1, sub2):
#redcov: calculates the reduced covariance matrix
#cov = covariance matrix
#sub1
#sub2

  dim = sub2-sub1
  redcov = np.zeros((dim,dim))
  for i in range(sub1,sub2):
    for j in range(sub1,sub2):
      redcov[i-sub1,j-sub1] = cov[i,j]/math.sqrt(cov[i,i]*cov[j,j])

  return redcov
