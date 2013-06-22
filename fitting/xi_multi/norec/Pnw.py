import math
import numpy as np
from scipy import interpolate
from spl import *

def Pnw(tnw, power, k, ns):
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
  A = P0/k0**ns
  Pnw = A*k**ns*tnw**2.0

  return Pnw
