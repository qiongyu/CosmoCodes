import numpy as np
import math

def makeksp(kmin,kmax,kstep):
#This function returns a k-grid (ksp) onto which a CAMB transfer function 
#will be splined to produce a power spectrum P(k_sp). This should be used
#in conjunction with the function readcamb.py (ksp should be the last
#argument pased to readcamb).
#	kmin = minimum k of the grid
#	kmax = maximum k of the grid
#	kstep = number of ksp values (they will be equally spaced in log)
  lkmin = math.log(kmin) #max and min k from _ztf file
  lkmax = math.log(kmax)
  kinc = (lkmax-lkmin)/kstep #increment in log(k)
  ksp = np.exp(np.arange(lkmin,lkmax,kinc)) #k values to spline onto

  return ksp
