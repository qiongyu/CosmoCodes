import numpy as np
from scipy import integrate
#These scripts are for adding in streaming (exp or Gaussian) and the Kaiser
#factor (1+beta*mu**2)**2

def pmu_gauss(mu, k, sig, beta):
  return (1.0+beta*mu**2)**2*np.exp(-(k*sig*mu)**2)

def pmu_exp(mu, k, sig, beta):
  return (1.0+beta*mu**2)**2/(1+(k*sig*mu)**2/2.)**2

def calc_stream(k, nk, beta, sigma, flag):
  stream = k*0.0+1.0

  if(flag=='exp'):
    for cc in range(0,nk):
      tempi = integrate.fixed_quad(pmu_exp,-1,1,args=(k[cc],sigma,beta),n=8)
      stream[cc] = 0.5*tempi[0] #0.5* to normalize Legendre L_0
  elif(flag=='Gauss'):
    for cc in range(0,nk):
      tempi = integrate.fixed_quad(pmu_gauss,-1,1,args=(k[cc],sigma,beta),n=8)
      stream[cc] = 0.5*tempi[0] #0.5* to normalize Legendre L_0
  
  return stream
