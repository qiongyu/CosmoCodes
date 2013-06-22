import numpy as np
from scipy import interpolate

def spllin(x_sp, y, x):
#spllog: spline y(x) onto y(x_sp) via log space
#x_sp = spline grid points
#y = y(x)
#x = abiscuss values

  spform = interpolate.splrep(x,y,s=0)
  y_out = interpolate.splev(x_sp,spform,der=0)

  return y_out

def spllog(x_sp, y, x):
#spllog: spline y(x) onto y(x_sp) via log space
#x_sp = spline grid points
#y = y(x)
#x = abiscuss values

  lx = np.log10(x)
  ly = np.log10(y)
  lxs = np.log10(x_sp)
  spform = interpolate.splrep(lx,ly,s=0)
  y_out = 10.0**interpolate.splev(lxs,spform,der=0)

  return y_out
