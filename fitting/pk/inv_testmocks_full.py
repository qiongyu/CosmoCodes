#ONE SCRIPT...TO FIT THEM ALL!!!!!!!!!!!!!!!

import math
import time
import sys
import getopt
from spl import *
from makeksp import *
from readcamb import *
from apodize import *
from sigmaR import *
from mcov import *
from stream import *
import matplotlib
from matplotlib import pyplot
from matplotlib.backends.backend_pdf import PdfPages
from scipy import optimize

def funcchi2_b(params):
  logb_0 = params
  b_0 = np.exp(logb_0)

  k_alpha = k*alpha_u
  nbasis0 = npt+1
  basis0_full = np.zeros((nbasis0-1,nk))

#SPLINING
  k3basis = k_alpha**3*xib0[0,:]
  modelb0 = spllin(k,k3basis,k_alpha)
  modelb0 = modelb0/k**3
  for i in range(0,nbasis0-1):
    basis0_full[i,:] = xib0[i+1,:]
  xibasis = basis0_full

#Normalize model to data at r=50
  modelb0 = modelb0*b_0

  model = np.zeros(nk_sub)
  spls = np.zeros(nk_sub)
  main = np.zeros(nk_sub)

  if(npt > 0):
    datanew = data[sub1:sub2] - modelb0[sub1:sub2]
    basis_alpha = xibasis[:,sub1:sub2]
    alpha_ = np.dot(np.dot(basis_alpha, covinv), np.transpose(basis_alpha))
    beta_ = np.dot(np.dot(basis_alpha, covinv), datanew)
    coeffs = np.linalg.solve(alpha_,beta_)

    npar = npt
    for i in range(0,npar):
      model = model + coeffs[i] * basis_alpha[i,:]
      spls = spls + coeffs[i] * basis_alpha[i,:]
      coeffs_out[i] = coeffs[i]

  main = modelb0
  model = model + modelb0[sub1:sub2]
  for i in range(0,nk_sub):
    model0[i] = model[i]
    spls0[i] = spls[i]
    main0[i] = main[i+sub1]

  residual = data[sub1:sub2] - model
#put in prior on b_0
  chi2 = np.dot(np.dot(residual, covinv), residual)

  return chi2

def funcchi2(params):
  print params
  alpha = params[0]
  logb_0 = params[1]
  b_0 = np.exp(logb_0)
#  snl, tilt = params

  k_alpha = k*alpha
  nbasis0 = npt+1
  print "Fitting ", nk_sub, " points"
  basis0_full = np.zeros((nbasis0-1,nk))

#SPLINING
  k3basis = k_alpha**3*xib0[0,:]
  modelb0 = spllin(k,k3basis,k_alpha)
  modelb0 = modelb0/k**3
  for i in range(0,nbasis0-1):
    basis0_full[i,:] = xib0[i+1,:]
  xibasis = basis0_full

  modelb0 = modelb0*b_0

  model = np.zeros(nk_sub)
  spls = np.zeros(nk_sub)
  main = np.zeros(nk_sub)

  if(npt > 0):
    datanew = data[sub1:sub2] - modelb0[sub1:sub2]
    
    basis_alpha = xibasis[:,sub1:sub2]
    alpha_ = np.dot(np.dot(basis_alpha, covinv), np.transpose(basis_alpha))
    beta_ = np.dot(np.dot(basis_alpha, covinv), datanew)
    coeffs = np.linalg.solve(alpha_,beta_)

    print "Coefficients: ", coeffs

    npar = npt
  
    for i in range(0,npar):
      model = model + coeffs[i] * basis_alpha[i,:]
      spls = spls + coeffs[i] * basis_alpha[i,:]
      coeffs_out[i] = coeffs[i]

  main = modelb0
  model = model + modelb0[sub1:sub2]
  for i in range(0,nk_sub):
    model0[i] = model[i]
    spls0[i] = spls[i]
    main0[i] = main[i+sub1]

  residual = data[sub1:sub2] - model
#put in prior on b_0
  chi2 = np.dot(np.dot(residual, covinv), residual)

  print "Chi^2 = ", chi2

  return chi2

narg = 0
try:
  opts, args = getopt.getopt(sys.argv[1:], "n:c:d:u:s:")
except getopt.GetoptError:
  sys.exit(2)
for opt, arg in opts:
  if opt == '-n':
    npt = int(arg)
    narg = narg+1
  if opt == '-s':
    stemp = arg
    narg = narg+1

if narg != 2:
  sys.exit("need -n -s")

sub1 = 34#34
sub2 = 69#69
nk_sub = sub2-sub1

print "Fitting A(k) as an inverse polynomial of order: ", npt

#================================ PARAMETERS ==================================
omega_m = 0.274
omega_b = 0.0457
omega_c = omega_m - omega_b
z_med = 0.55
omega_mz = omega_m*(1.+z_med)**3/(omega_m*(1.+z_med)**3+1.-omega_m)

bias = 1.6
beta = omega_mz**0.55/bias

ns = 0.95
sigma8_lin = 0.8
signl_ap = 8.0
#==============================================================================
#================================ READ STUFF ==================================
#read in no wiggles T(k)
knw, tnw = np.loadtxt('dr9_full.ztf_nowigs', unpack=True)
#read in the mock file names
a = np.loadtxt('red_norec.txt',dtype=({'names': ['infile'],'formats': ['S100']}))
infile = a['infile']
nf = np.size(infile)
if(stemp == 'f'):
  snl = signl_ap
if(stemp == 'z'):
  snl = 0.0
if(stemp == 'p'):
  snl = signl_ap+2.
if(stemp == 'm'):
  snl = signl_ap-2.
if(stemp == 'i'):
  snl = 1000.0
ss = "_snl"+str(int(snl))

#OUTPUT FILES
figout = PdfPages('ninv_poly'+str(npt)+ss+'.pdf')
fout = open("nalpha_poly"+str(npt)+ss+".txt","w")
fouts = open("nstd_poly"+str(npt)+ss+".txt","w")
#==============================================================================

ksp = makeksp(0.15e-4,200.0,1000.0)
psp = readcamb("dr9_full.ztf", ns, 8.0, sigma8_lin, omega_b, omega_c, ksp)
print 'sigmaR =', sigmaR(ksp, psp, 8.0)

k, tt = np.loadtxt(infile[0],unpack=True)
dex = np.squeeze(np.where(k > 0.))
k = k[dex]
psp = spllog(k, psp, ksp)
nk = np.size(k)

#spline no wiggle transfer function onto same k-grid as P(k), then apodize
tspnw = spllog(k, tnw, knw)
#DIVIDE data by pspnw (BAO-less P(k)), pk_mod should already be correctly
#normalized
pspnw, pk_mod = apodize(tspnw, psp, k, ns, snl)

#UNCOMMENT if want model with RSD
#rsd = calc_stream(k, nk, beta, 3.0, 'exp')
#pk_mod = pk_mod*rsd
#np.savetxt("pk_mod_rsd.txt",zip(k,pk_mod))

cov = mcov00("red_norec.txt","red_norec_avg.pk",pspnw)
print "Fitting using the mock Covariance matrix..."
covinv = np.linalg.inv(cov[sub1:sub2,sub1:sub2])

perr = np.zeros(nk)
for i in range(0,nk):
  perr[i] = math.sqrt(cov[i,i])

xib0 = np.zeros((npt+1,nk))
for i in range(1,npt+1):
    btemp = k**(1.0*i)
    xib0[i,:] = btemp
xib0[0,:] = pk_mod

matplotlib.rcParams['font.size']=7
font0 = matplotlib.font_manager.FontProperties()
font0.set_size(12)
fig = pyplot.figure()

alpha_use = np.arange(0.6, 1.4, 0.01)
dmini = 0.0001
alpha_mini = np.arange(0.6,1.4,dmini)
alpha_u = 1.0
nalpha = np.size(alpha_use)
chiopt = np.zeros(nalpha)
b0opt = np.zeros(nalpha)
b0bf = np.zeros(nf)
chibf = np.zeros(nf)
dchi = np.zeros(nf)
abf = np.zeros(nf)
da = np.zeros(nf)
dsdanp = np.zeros(nf)
dsda = np.zeros(nf)
pat1 = np.zeros(nf)
ap = np.zeros((npt,nf))
coeffs_out = np.zeros(npt+1)

bseed = math.log(bias**2.)
for i in range(0,nf):
    
  print "Now fitting...", i, infile[i]

  kin, data = np.loadtxt(infile[i], unpack=True, usecols=(0,1))
  dex = np.squeeze(np.where(kin > 0.))
  kin = kin[dex]
  data = data[dex]/pspnw #Divide data by BAO-less P(k)

  t_beg = time.clock()
  for j in range(0,nalpha):
    model0 = np.zeros(nk_sub)
    spls0 = np.zeros(nk_sub)
    main0 = np.zeros(nk_sub)
#store model in this
    alpha_u = alpha_use[j]
    aopt = optimize.fmin(funcchi2_b, [bseed], ftol=1e-4, xtol=1e-7, disp=False)
    chiopt[j] = funcchi2_b(aopt[0])
    b0opt[j] = np.exp(aopt[0])
  logb1 = np.log(b0opt[40])
 
  cd = chiopt.argmin()
  ax = fig.add_subplot(221)
  pyplot.plot(alpha_use, chiopt-chiopt[cd], marker='x')
  pyplot.xlabel(r"$\alpha$")
  pyplot.ylabel(r"$\Delta \chi^2$ from min")
  pyplot.figtext(0.15,0.87,r'min $\chi^2 =$ '+str(chiopt[cd]))

  params = [1.0, logb1]
  opt = optimize.fmin(funcchi2, params, ftol=1e-4, xtol=1e-7)
  print "On ", i
  chibf[i] = funcchi2(opt)
  dchi[i] = chibf[i] - chiopt[40]
  for pp in range(0,npt):
    ap[pp,i] = coeffs_out[pp]
  abf[i] = opt[0]
  b0bf[i] = np.exp(opt[1])
  ax = fig.add_subplot(222)
  pyplot.semilogx(kin[sub1:sub2], data[sub1:sub2], 'kx', ms=4, mew=1.2)
  pyplot.semilogx(k[sub1:sub2], model0, 'r-')
  pyplot.semilogx(k[sub1:sub2], spls0, 'r--')
  pyplot.semilogx(k[sub1:sub2], main0, 'r:')
  pyplot.errorbar(kin[sub1:sub2], data[sub1:sub2],fmt=None, yerr=perr[sub1:sub2])
  pyplot.figtext(0.57,0.65,r'$b_0 = $'+str(b0bf[i]))
  for pp in range(0,npt):
    pyplot.figtext(0.57,0.65-0.02*(pp+1),r'$a_{%d} = {%f}$' % (pp+1, ap[pp,i]))
  pyplot.figtext(0.57,0.65-0.02*(npt+1),r'$\alpha = $'+str(opt[0]))
  pyplot.figtext(0.57,0.65-0.02*(npt+2),r'$\chi^2 = $'+str(chibf[i]))
  pyplot.xlim(0.02,0.3)
  pyplot.title("Best fit")

  chi1 = funcchi2([1.0,logb1])
  b01 = b0opt[40]
  for pp in range(0,npt):
    ap[pp,i] = coeffs_out[pp]
  ax = fig.add_subplot(223)
  pyplot.semilogx(kin[sub1:sub2], data[sub1:sub2], 'kx', ms=4, mew=1.2)
  pyplot.semilogx(k[sub1:sub2], model0, 'r-')
  pyplot.semilogx(k[sub1:sub2], spls0, 'r--')
  pyplot.semilogx(k[sub1:sub2], main0, 'r:')
  pyplot.errorbar(kin[sub1:sub2], data[sub1:sub2],fmt=None, yerr=perr[sub1:sub2])
  pyplot.figtext(0.15,0.2,r'$b_0 = $'+str(b01))
  for pp in range(0,npt):
    pyplot.figtext(0.15,0.2-0.02*(pp+1),r'$a_{%d} = {%f}$' % (pp+1, ap[pp,i]))
  pyplot.figtext(0.15,0.2-0.02*(npt+1),r'$\chi^2 = $'+str(chi1))
  pyplot.xlim(0.02,0.3)
  pyplot.title(r"$\alpha = 1$")

  ax = fig.add_subplot(224)
  pnpns = np.exp(-chiopt/2) #Prob without prior
  chip = chiopt+(alpha_use-1.0)**2/0.0225 #chi^2 with prior
  chipl = chiopt+(np.log(alpha_use)/0.15)**2 #chi^2 with .15 prior on ln(alpha)
  chipi = chiopt+(1.0/alpha_use-1.0)**2/0.0225 #chi^2 with .1 prior on 1/alpha
  pwpns = np.exp(-chip/2) #Prob with prior
  pwplns = np.exp(-chipl/2) #Prob with ln(alpha) prior
  pwpins = np.exp(-chipi/2) #Prob with 1/alpha prior
  pnp = spllog(alpha_mini, pnpns, alpha_use) #spline probabilities onto tighter
  pwp = spllog(alpha_mini, pwpns, alpha_use) #alpha grid...
  pwpl = spllog(alpha_mini, pwplns, alpha_use) #alpha grid...
  pwpi = spllog(alpha_mini, pwpins, alpha_use) #alpha grid...
  pnptot = np.sum(pnp)*dmini
  pwptot = np.sum(pwp)*dmini #normalization factor 
  pwpltot = np.sum(pwpl)*dmini
  pwpitot = np.sum(pwpi)*dmini

  pnp = pnp/pnptot
  pwp = pwp/pwptot
  pwpl = pwpl/pwpltot
  pwpi = pwpi/pwpitot
  pat1[i] = np.sum(pnp[np.where(alpha_mini<1.0)])*dmini
  print pnptot, pwptot, pwpltot, pwpitot
  mnp = np.sum(alpha_mini*pnp)*dmini #mean alpha without prior
  mwp = np.sum(alpha_mini*pwp)*dmini #mean alpha with prior
  mwpl = np.sum(alpha_mini*pwpl)*dmini #mean alpha with ln prior
  mwpi = np.sum(alpha_mini*pwpi)*dmini #mean alpha with inverse prior
  sdnp = math.sqrt(np.sum((alpha_mini-mnp)**2*pnp)*dmini) #std without prior
  sdwp = math.sqrt(np.sum((alpha_mini-mwp)**2*pwp)*dmini) #std with prior
  sdwpl = math.sqrt(np.sum((alpha_mini-mwpl)**2*pwpl)*dmini) #std ln prior
  sdwpi = math.sqrt(np.sum((alpha_mini-mwpi)**2*pwpi)*dmini) #std inv prior
  da[i] = mwpl - mwp
  dsdanp[i] = sdnp
  dsda[i] = sdwpl
#  skwp = np.sum(((alpha_mini-mwp)/sdwp)**3*pwp)/pwptot
#  sknp = np.sum(((alpha_mini-mnp)/sdnp)**3*pnp)/pnptot
  pyplot.plot(alpha_mini, pnp, label="no prior")
  pyplot.plot(alpha_mini, pwp, label="with prior")
  pyplot.plot(alpha_mini, pwpl, label=r"with $\log(\alpha)$ prior")
  pyplot.plot(alpha_mini, pwpi, label=r"with $\alpha^{-1}$ prior")
  pyplot.figtext(0.77,0.43,"Without prior: ", weight='bold')
  pyplot.figtext(0.77,0.41,r"Mean $\alpha = %4.3f \pm %4.3f$" % (mnp,sdnp))
  pyplot.figtext(0.77,0.39,"With prior: ", weight='bold')
  pyplot.figtext(0.77,0.37,r"Mean $\alpha = %4.3f \pm %4.3f$" % (mwp,sdwp))
  pyplot.figtext(0.77,0.35,r"With $\log(\alpha)$ prior: ", weight='bold')
  pyplot.figtext(0.77,0.33,r"Mean $\alpha = %4.3f \pm %4.3f$" % (mwpl,sdwpl))
  pyplot.figtext(0.77,0.31,r"With $\alpha^{-1}$ prior: ", weight='bold')
  pyplot.figtext(0.77,0.29,r"Mean $\alpha = %4.3f \pm %4.3f$" % (mwpi,sdwpi))
  pyplot.xlabel(r"$\alpha$")
  pyplot.ylabel(r"P($\alpha$)")
  pyplot.legend(loc='upper left')

  figout.savefig()
  pyplot.clf()

  t_end = time.clock()
  print "It took: ", t_end-t_beg, "seconds to run through the ring."
  fout.write(str(abf[i])+' '+str(chibf[i])+' '+str(b0bf[i])+'\n')
  fouts.write(str(abf[i])+' '+str(dsda[i])+' '+str(dsdanp[i])+' '+str(pat1[i])+'\n')
  #fouts.write(str(mwpl)+' '+str(dsda[i])+' '+str(dsdanp[i])+' '+str(pat1[i])+'\n')

figout.close()
fout.close()
fouts.close()

