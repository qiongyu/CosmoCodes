import math
import numpy as np
from spl import *
from sigmaR import *
from scipy import interpolate

def readcamb(tfile, n, R, sigR, omb, omc, ksp):
#readcmb: Reads in a CAMB output file, constructs P(k) and splines onto a different set of
#k's if requested (ksp)
# tfile = Name of cmbfast output file (containing the transfer fxn)
# n = spectral tilt
# R = radius to which power spectrum is to be normalized (corresponding to sigma_R)
# sigR = sigma_R value to which power spectrum is to be normalized
# omb = omega_b
# omc = omega_c
# ksp = if splining onto a different basis is desired, these are the abiscuss values, the input
# 	into sigmaR must be LOG spaced in k, so if the k's in the cmbfast output file are not
#	equally spaced in log, ksp MUST be specified, set equal to 0.0 if no splining is required

  kread, Tc, Tb, Tgamma, Tnu, moo, moo2= np.loadtxt(tfile, unpack=True)
  omt = omb + omc
# Construct the transfer fxn...
  xfer = omb/omt * Tb + omc/omt * Tc
  Pread = xfer**2 * kread**n
  k = kread
  power = Pread

# Spline if necessary
  if np.size(ksp) != 1:
    k = ksp
    power = spllog(ksp, Pread, kread)

# Normalize to sigma_R
  sigma_R = sigmaR(k, power, R)
  bias = sigR/sigma_R
  power = power*bias**2

  return power

