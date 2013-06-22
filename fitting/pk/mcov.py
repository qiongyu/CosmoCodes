import math
import numpy as np

def mcov00(flist, avgfile, pnw):
  a = np.loadtxt(flist,dtype=({'names': ['infile'],'formats': ['S100']}))

  infile = a['infile']
  nf = np.size(infile)

  ravg, xiavg = np.loadtxt(avgfile, unpack=True)
  xiavg = xiavg/pnw

  nr = np.size(ravg)
  xi = np.zeros((nr, nf))
  cov = np.zeros((nr,nr))

  for i in range(0,nf):
    rin, xiin = np.loadtxt(infile[i], unpack=True, usecols=(0,1))
    dex = np.squeeze(np.where(rin > 0.))
    rin = rin[dex]
    xiin = xiin[dex]/pnw #Divide all mocks by BAO-less P(k)!!
    #j1, j2, rin, xiin = np.loadtxt(infile[i], unpack=True)
    xi[:,i] = xiin

  for i in range(0,nr):
    for j in range(0,nr):
      diffi = xi[i,:]-xiavg[i]
      diffj = xi[j,:]-xiavg[j]
      cov[i,j] = sum(diffi*diffj)
  cov = cov/(nf-1.0)
#  print cov
  return cov

#moo = mcov00("real_norec_cor.txt", "real_norec_cor_avg.xi")
