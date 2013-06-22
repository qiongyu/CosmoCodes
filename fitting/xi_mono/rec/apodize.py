import math
import numpy as np
from scipy import interpolate
from spl import *
from matplotlib import pyplot

def apodize(tnw, power, k, n, signl):
#apodize: Apodizes the input power spectrum using the input no wiggle power
#spectrum
#tnw = no wiggle transfer function (must be on same k-grid as power)
#power = input power spectrum
#k = k values
#n = spectral tilt
#signl = sigma_nl

  dex = np.squeeze(np.where(k > 9.0e-5))
  k0 = k[dex[0]]
  P0 = power[dex[0]]
  A = P0/k0**n
  pnw = A*k**n*tnw**2.0

  gauss = np.exp(-(k*signl)**2.0/2.0)
  p_out = gauss*power + (1.0 - gauss)*pnw

  #pyplot.plot(np.log10(k),np.log10(pnw),'r')
  #pyplot.plot(np.log10(k),np.log10(p_out),'k')
  #pyplot.savefig('pk.ps')
  #pyplot.clf()

  p02 = spllog(0.2,p_out, k)
  print 'P(k=0.2) =', p02

  return p_out
