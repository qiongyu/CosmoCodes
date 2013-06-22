import numpy as np

def zdistort_sep(mu, k, sig_s, type):
  if(type=='exp'):
    FoG = 1.0/(1.+(k*sig_s*mu)**2/2.)**2
  if(type=='gauss'):
    FoG = np.exp(-(k*sig_s*mu)**2)
  K1 = FoG
  K2 = 2.0*FoG*mu**2.
  K3 = FoG*mu**4.
  return K1, K2, K3

def zP_gauss_sep(x, plin, pnw, k, spar, type):
#Calculates Kaiser*FoG and P_dewiggled at the abscissa x values for 
#Gaussian quadrature at all k
#Returns 2D matrix with these values at x and k
  nk = np.size(k)
  sigperp = spar[0]
  sigpar = spar[1]
  sig_s = spar[2]

  npt = np.size(x)
  zdist = np.zeros((npt,nk,3))
  p_dw = np.zeros((npt,nk))
  for cc in range(0,nk):
    zd1, zd2, zd3 = zdistort_sep(x, k[cc], sig_s, type)
    zdist[:,cc,0] = zd1
    zdist[:,cc,1] = zd2
    zdist[:,cc,2] = zd3
    p_dw[:,cc] = Pdw(x, k[cc], plin[cc], pnw[cc], sigperp, sigpar)

  return zdist, p_dw

def TwoDP_sep(Pdw, zdist):
  P2D = zdist*0.0
  P2D[:,:,0] = zdist[:,:,0]*Pdw
  P2D[:,:,1] = zdist[:,:,1]*Pdw
  P2D[:,:,2] = zdist[:,:,2]*Pdw

  return P2D

def zdistort(mu, k, sig_s, beta, type):
  if(type=='exp'):
    FoG = 1.0/(1.+(k*sig_s*mu)**2)**2
  if(type=='gauss'):
    FoG = np.exp(-(k*sig_s*mu)**2)
  Kaiser = (1.+beta*mu**2)**2

  zdist = FoG*Kaiser
  return zdist

def Pdw(mu, k, Plin, Pnw, sigperp, sigpar):
  gau = np.exp(-0.5*(k**2*(mu**2*(sigpar**2-sigperp**2)+sigperp**2)))
  Pdw = gau*(Plin-Pnw) + Pnw

  return Pdw

def zP_gauss(x, plin, pnw, k, spar, type):
#Calculates Kaiser*FoG and P_dewiggled at the abscissa x values for 
#Gaussian quadrature at all k
#Returns 2D matrix with these values at x and k
  nk = np.size(k)
  sigperp = spar[0]
  sigpar = spar[1]
  sig_s = spar[2]
  beta = spar[3]

  npt = np.size(x)
  zdist = np.zeros((npt,nk))
  p_dw = np.zeros((npt,nk))
  for cc in range(0,nk):
    zdist[:,cc] = zdistort(x, k[cc], sig_s, beta, type)
    p_dw[:,cc] = Pdw(x, k[cc], plin[cc], pnw[cc], sigperp, sigpar)

  return zdist, p_dw

def TwoDP(Pdw, zdist):
  P2D = Pdw*zdist

  return P2D

def TwoDPn(mu, Pdw, zdist, params, nbar):
  par0 = params[0]
  par1 = params[1]
  par2 = params[2]

  Pmu = zdist*(par0*Pdw + par1/nbar)
  Ppn00 = (Pmu + par2/nbar)**2
  L2 = 0.5*(3.*mu**2-1.)
  
  Ppn02 = Ppn00*L2
  Ppn22 = Ppn02*L2

  return Ppn00, Ppn02, Ppn22
