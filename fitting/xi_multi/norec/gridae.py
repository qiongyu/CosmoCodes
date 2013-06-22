#ONE SCRIPT...TO FIT THEM ALL!!!!!!!!!!!!!!!

import math
import time
import sys
import getopt
from spl import *
from makeksp import *
from readcamb import *
from Pnw import *
from TwoDP import *
from sigmaR import *
from covcov import *
from p2xi import *
from mcov import *
from xipoles import *
import matplotlib
from matplotlib import pyplot
from matplotlib import mlab
from matplotlib.backends.backend_pdf import PdfPages
from scipy import optimize

def funcchi2_b(params):
  logb_0 = params[0]
  b_0 = np.exp(logb_0)
  beta = params[1]

  xi0_fid = xi0_1 + xi0_2*beta + xi0_3*beta**2.
  xi2_fid = xi2_1 + xi2_2*beta + xi2_3*beta**2.
  xi4_fid = xi4_1 + xi4_2*beta + xi4_3*beta**2.

  xi0_obs = np.zeros((nr_sub,nmuo2))
  xi2_obs = np.zeros((nr_sub,nmuo2))
  xi4_obs = np.zeros((nr_sub,nmuo2))

  for i in range(0,nmuo2):
    xi0_obs[:,i] = spllin(r_obs[:,i],r_fid2*xi0_fid,r_fid) / r_obs2[:,i]
    xi2_obs[:,i] = spllin(r_obs[:,i],r_fid2*xi2_fid,r_fid) / r_obs2[:,i] * L2_obs[i]
    xi4_obs[:,i] = spllin(r_obs[:,i],r_fid2*xi4_fid,r_fid) / r_obs2[:,i] * L4_obs[i]
  xi_obs = xi0_obs + xi2_obs + xi4_obs
  xi_obs0 = xi_obs*0.0
  xi_obs2 = xi_obs*0.0
  for i in range(0,nmuo2):
    xi_obs0[:,i] = xi_obs[:,i]*weighto2[i]
    xi_obs2[:,i] = xi_obs[:,i]*L2_fid[i]*weighto2[i]
  xi0 = np.sum(xi_obs0,1)
  xi2 = 5.*np.sum(xi_obs2,1)

  modelb0 = xi0*b_0
  modelb2 = xi2

  fullmodel = np.zeros(2*nr_sub)
  for i in range(0,nr_sub):
    fullmodel[i] = modelb0[i]
    fullmodel[i+nr_sub] = modelb2[i]

  model = np.zeros(2*nr_sub)
  spls = np.zeros(2*nr_sub)
  main = np.zeros(2*nr_sub)

  datanew = data - fullmodel
  basis_alpha = Ar
  alpha_ = np.dot(np.dot(basis_alpha, covinv), np.transpose(basis_alpha))
  beta_ = np.dot(np.dot(basis_alpha, covinv), datanew)
  coeffs = np.linalg.solve(alpha_,beta_)

  npar = nptm+nptq
  for i in range(0,npar):
    model = model + coeffs[i] * basis_alpha[i,:]
    spls = spls + coeffs[i] * basis_alpha[i,:]
    ap[i] = coeffs[i]

  main = fullmodel
  model = model + fullmodel
  for i in range(0,nr_sub):
    model0[i] = model[i]
    spls0[i] = spls[i]
    main0[i] = main[i]
    model2[i] = model[i+nr_sub]
    spls2[i] = spls[i+nr_sub]
    main2[i] = main[i+nr_sub]

  residual = data - model
#put in prior on b_0
  chi2 = np.dot(np.dot(residual, covinv), residual) + ((beta-beta_c)/0.2)**2

  return chi2

narg = 0
try:
  opts, args = getopt.getopt(sys.argv[1:], "m:q:c:d:u:")
except getopt.GetoptError:
  sys.exit(2)
for opt, arg in opts:
  if opt == '-m':
    nptm = int(arg)
    narg = narg+1
  if opt == '-q':
    nptq = int(arg)
    narg = narg+1
  if opt == '-c':
    if(arg == 'g'):
      gom = 'gc'
      #read in parms
      params = np.loadtxt("maxlike_result",usecols=(0,1,2,3,4))
      par0 = params[0]
      par1 = params[1]
      par2 = params[2]
      sig_s = params[3]
      par4 = params[4]
      #read in n(z) and mesh(z) information
      nbar, mesh = np.loadtxt("nz_North.txt", unpack=True)
    if(arg == 'm'):
      gom = 'mc'
    narg = narg+1
  if opt == '-d':
    dd = arg
    if(arg == '30'):
      sub1 = 6
    if(arg == '20'):
      sub1 = 4
    if(arg == '50'):
      sub1 = 11
    if(arg == '70'):
      sub1 = 16
    narg = narg+1
  if opt == '-u':
    if(arg == '150'):
      sub2 = 37
      uu = "-150"
    if(arg == '200'):
      sub2 = 50
      uu = "-200"
    narg = narg+1

if narg != 5:
  sys.exit("need -m -q -c -d -u")

print "Fitting A(r) as a polynomial of order:", nptm, "for monopole and", nptq, "for quadrupole"
#================================ PARAMETERS ==================================
rstep = 4.0
rmax = 603.0
rmin = 4.0
nr_sub = sub2-sub1

r = np.arange(rmin,rmax,rstep)
nr = np.size(r)
r_sub = r[sub1:sub2]
print r_sub

omega_m = 0.274
omega_b = 0.0457
omega_c = omega_m - omega_b
z_med = 0.57
omega_mz = omega_m*(1.+z_med)**3/(omega_m*(1.+z_med)**3+1.-omega_m)

bias = 1.6
bias2 = bias*bias
beta_c = omega_mz**0.55/bias
print "beta: ", beta_c

ns = 0.95
sigma8_lin = 0.8
## For covariance matrix...
sigperp = 6.0
sigpar = 11.0
## For model...
_sigperp = 6.0
_sigpar = 11.0
_sig_s = 3.0
#==============================================================================
#================================ READ STUFF ==================================
#read in no wiggles T(k)
knw, tnw = np.loadtxt('dr9_full.ztf_nowigs', unpack=True)
#read in kernel J_0
kr, J0, J2, J4 = np.loadtxt('Jfilt.dat', unpack=True)
#read in the mock file names
a = np.loadtxt('red_norec.txt',dtype=({'names': ['infile'],'formats': ['S100']}))
infile = a['infile']
nf = np.size(infile)
#read in weights for Gaussian quadrature
x, weight = np.loadtxt("gauss.txt", unpack=True)
#define fiducial coordinate system
nmuo2 = np.size(x)/2
mu_fid2 = x[0:nmuo2]**2.
weighto2 = weight[0:nmuo2]
r_fid = r
r_fid2 = r_fid**2.
L2_fid = 0.5*(3.*mu_fid2-1.)

#OUTPUT FILES
fout = open("cov_diag","w")
figout = PdfPages('ngrid'+dd+uu+'_poly'+str(nptm)+str(nptq)+'_'+gom+'.pdf')
#==============================================================================

ksp = makeksp(0.15e-4,200.0,1000.0)
psp = readcamb("dr9_full.ztf", ns, 8.0, sigma8_lin, omega_b, omega_c, ksp)
print 'sigmaR =', sigmaR(ksp, psp, 8.0)

#spline no wiggle transfer function onto same k-grid as P(k)
tspnw = spllog(ksp, tnw, knw)

#apodize psp
p_nw = Pnw(tspnw, psp, ksp, ns)

kcovdex = np.squeeze(np.where(ksp < 1000.0/max(r)))
ksp_c = ksp[kcovdex]
nsims = 1.0
totsim = nf*1.

if(gom == 'gc'):
  spars = [sigperp,sigpar,sig_s,beta_c]
  zdist_g, pdw_g = zP_gauss(x, psp*bias2, p_nw*bias2, ksp, spars, 'exp')

  par_in = [par0,par1,par2]
  pv00, pv02, pv22 = volint(x, weight, par_in, zdist_g[:,kcovdex], pdw_g[:,kcovdex], ksp_c, spars, nbar, mesh)
  cov = ( cov02(ksp_c, pv00, pv22, pv02, r_sub) + par4 ) / (nsims*1.0)
  print "Fitting using the Gaussian Covariance matrix..."
if(gom == 'mc'):
  cov = mcov02("red_norec.txt","red_norec_avg.xi",sub1,sub2) / (nsims*1.0)
  cov = (totsim - 1.) / (totsim - 2.*nr_sub - 2.) * cov
  print "Fitting using the mock Covariance matrix..."
covinv = np.linalg.inv(cov)

#write out the diagonals of the covariance matrix (variance values)
xierr = np.zeros(2*nr_sub)
for i in range(0,nr_sub):
  xierr[i] = math.sqrt(cov[i,i])
  xierr[i+nr_sub] = math.sqrt(cov[i+nr_sub,i+nr_sub])
  fout.write(str(r[i+sub1])+' '+str(xierr[i])+' '+str(xierr[i+nr_sub])+'\n')
fout.close()

rc = redcov(cov)

Ar = np.zeros((nptm+nptq,2*nr_sub))
for i in range(0,nptm):
  Atemp = r**(i-2.)
  Ar[i,0:nr_sub] = Atemp[sub1:sub2]
for i in range(nptm,nptm+nptq):
  Atemp = r**(i-nptm-2.)
  Ar[i,nr_sub:] = Atemp[sub1:sub2]

_spars = [_sigperp,_sigpar,_sig_s]
_zdist_g, _pdw_g = zP_gauss_sep(x, psp, p_nw, ksp, _spars, 'exp')
pm = TwoDP_sep(_pdw_g, _zdist_g)
xi0_1, xi0_2, xi0_3 = xi_ell_sep(x,weight,pm,ksp,J0,kr,r,0)
xi2_1, xi2_2, xi2_3 = xi_ell_sep(x,weight,pm,ksp,J2,kr,r,2)
xi4_1, xi4_2, xi4_3 = xi_ell_sep(x,weight,pm,ksp,J4,kr,r,4)

matplotlib.rcParams['font.size']=16
font0 = matplotlib.font_manager.FontProperties()
font0.set_size(20)
fig = pyplot.figure()

natry = 241
netry = 121
almin = 0.7
almax = 1.3
epmin = -0.3
epmax = 0.3
alpha_try = np.linspace(almin,almax,natry)
epsilon_try = np.linspace(epmin,epmax,netry)
ntot = natry*netry
chi_try = np.zeros(ntot)
ai = np.linspace(almin, almax, 100)
ei = np.linspace(epmin, epmax, 100)

bseed = math.log(bias2)
for i in range(0,2):
#for i in range(0,nf):
  print "Now fitting...", i, infile[i]
  xiin0, xiin2 = np.loadtxt(infile[i], usecols=(1,2), unpack=True)

  if(i<10):
    foutg = open('grid000'+str(i)+'.dat',"w")
  elif(i>=10 and i<100):
    foutg = open('grid00'+str(i)+'.dat',"w")
  else:
    foutg = open('grid0'+str(i)+'.dat',"w")

  data = np.zeros(2*nr_sub)
  for j in range(0,nr_sub):
    data[j] = xiin0[j+sub1]
    data[j+nr_sub] = xiin2[j+sub1]

  ap = np.zeros(nptm+nptq)

  model0 = np.zeros(nr_sub)
  spls0 = np.zeros(nr_sub)
  main0 = np.zeros(nr_sub)
  model2 = np.zeros(nr_sub)
  spls2 = np.zeros(nr_sub)
  main2 = np.zeros(nr_sub)

  ag = np.zeros(ntot)
  eg = np.zeros(ntot)
  count = 0
  for j in range(0,natry):
    alpha = alpha_try[j]
    for jj in range(0,netry):
      epsilon = epsilon_try[jj]
      print alpha, epsilon
      opep = epsilon + 1.
      mu_obs2 = 1./(1. + (1./mu_fid2-1)/opep**6.)
      L2_obs = 0.5*(3.*mu_obs2-1.)
      L4_obs = 0.125*(35.*mu_obs2**2. - 30.*mu_obs2 + 3.)
      r_obs = np.zeros((nr_sub,nmuo2))
      for jjj in range(0,nr_sub):
        r_obs[jjj,:] = r_fid[jjj+sub1]*alpha*np.sqrt(opep**4.*mu_fid2 + (1.-mu_fid2)/opep**2.)
      r_obs2 = r_obs**2.

      partry = [bseed,beta_c]
      opttry = optimize.fmin(funcchi2_b, partry, ftol=1e-6, xtol=1e-7, maxfun=20000, maxiter=20000)
      chi_try[count] = funcchi2_b(opttry)
      ag[count] = alpha_try[j]
      eg[count] = epsilon_try[jj]
      foutg.write(str(ag[count])+' '+str(eg[count])+' '+str(chi_try[count])+'\n')
      count = count + 1

  foutg.close()

  chi_try = chi_try - chi_try.min()
  chii = mlab.griddata(ag, eg, chi_try, ai, ei)
  pyplot.imshow(chii,extent=[almin,almax,epmin,epmax],origin='lower',vmin=0.0,vmax=30.0)
  cb = pyplot.colorbar()
  cb.set_label(r"$\chi^2(\alpha,\epsilon)$")
  pyplot.contour(ai,ei,chii,colors='k',origin='lower',levels=[2.3,6.18,11.8,19.5,24.])
  pyplot.xlabel(r"$\alpha$")
  pyplot.ylabel(r"$\epsilon$")

  figout.savefig()
  pyplot.clf()

figout.close()
