import math
import numpy as np

def mcov02(flist, avgfile, sub1, sub2):
  a = np.loadtxt(flist,dtype=({'names': ['infile'],'formats': ['S100']}))

  infile = a['infile']
  nf = np.size(infile)

  ravg, xiavg0, xiavg2 = np.loadtxt(avgfile, usecols=(0,1,2), unpack=True)
  ravg_sub = ravg[sub1:sub2]
  xiavg0_sub = xiavg0[sub1:sub2]
  xiavg2_sub = xiavg2[sub1:sub2]

  nr = np.size(ravg_sub)
  xi0 = np.zeros((nr, nf))
  xi2 = np.zeros((nr, nf))
  cov = np.zeros((2*nr,2*nr))

  for i in range(0,nf):
    xiin0, xiin2 = np.loadtxt(infile[i], usecols=(1,2), unpack=True)
    #j1, j2, rin, xiin = np.loadtxt(infile[i], unpack=True)
    xi0[:,i] = xiin0[sub1:sub2]
    xi2[:,i] = xiin2[sub1:sub2]

  for i in range(0,nr):
    for j in range(0,nr):
      diffi0 = xi0[i,:]-xiavg0_sub[i]
      diffj0 = xi0[j,:]-xiavg0_sub[j]
      diffi2 = xi2[i,:]-xiavg2_sub[i]
      diffj2 = xi2[j,:]-xiavg2_sub[j]
      cov[i,j] = sum(diffi0*diffj0)
      cov[i+nr,j+nr] = sum(diffi2*diffj2)
      cov[i,j+nr] = sum(diffi0*diffj2)
      cov[j+nr,i] = cov[i,j+nr]
  cov = cov/(nf-1.0)
#  print cov
  return cov

#moo = mcov00("real_norec_cor.txt", "real_norec_cor_avg.xi")
