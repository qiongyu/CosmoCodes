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
from covnew import *
from p2xi import *
from mcov import *
from stream import *
import matplotlib
from matplotlib import pyplot
from matplotlib.backends.backend_pdf import PdfPages
from scipy import optimize

def funcchi2_b(params):
  logb_0 = params
  b_0 = np.exp(logb_0)
#  snl, tilt = params

  r_alpha = r/alpha_u
  nbasis0 = npt+1
  basis0_full = np.zeros((nbasis0-1,nr))

#SPLINING
  r2basis = r_alpha**2*xib0[0,:]
  modelb0 = spllin(r,r2basis,r_alpha)
  modelb0 = modelb0/r**2
  for i in range(0,nbasis0-1):
    basis0_full[i,:] = xib0[i+1,:]
  xibasis = basis0_full

#Normalize model to data at r=50
  m30 = spllin(50,r2basis,r_alpha)/50.0**2 #model at r=50
  deg = 2
  lfit = int(np.round((40.0-rmin)/rstep,0)) #index of r=40
  ufit = int(np.round((60.0-rmin)/rstep,0)) #index of r=60
  nco = np.polyfit(r[lfit:ufit],r[lfit:ufit]**2*data[lfit:ufit],deg)
  d30 = nco[0]*50.0**2 + nco[1]*50.0 + nco[2] #fit to data at r=50
  d30 = d30/50.0**2
  modelb0 = modelb0/m30*d30*b_0
  #modelb0 = modelb0*b_0

  model = np.zeros(nr_sub)
  spls = np.zeros(nr_sub)
  main = np.zeros(nr_sub)

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
  for i in range(0,nr_sub):
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

  r_alpha = r/alpha
  nbasis0 = npt+1
  print "Fitting ", nr_sub, " points"
  basis0_full = np.zeros((nbasis0-1,nr))

#SPLINING
  r2basis = r_alpha**2*xib0[0,:]
  modelb0 = spllin(r,r2basis,r_alpha)
  modelb0 = modelb0/r**2
  for i in range(0,nbasis0-1):
    basis0_full[i,:] = xib0[i+1,:]
  xibasis = basis0_full

#Normalize model to data at r=50
  m30 = spllin(50,r2basis,r_alpha)/50.0**2 #model at r=50
  deg = 2
  lfit = int(np.round((40.0-rmin)/rstep,0)) #index of r=40
  ufit = int(np.round((60.0-rmin)/rstep,0)) #index of r=60
  nco = np.polyfit(r[lfit:ufit],r[lfit:ufit]**2*data[lfit:ufit],deg)
  d30 = nco[0]*50.0**2 + nco[1]*50.0 + nco[2] #fit to data at r=50
  d30 = d30/50.0**2
  modelb0 = modelb0/m30*d30*b_0
#  modelb0 = modelb0*b_0

  model = np.zeros(nr_sub)
  spls = np.zeros(nr_sub)
  main = np.zeros(nr_sub)

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
  for i in range(0,nr_sub):
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
  if opt == '-c':
    if(arg == 'g'):
      gom = 'gc'
      #read in n(z) and mesh(z) information
      nbar, mesh = np.loadtxt("nz_North.txt", unpack=True)
      #read in parms
      params = np.loadtxt("maxlike_result_ff",usecols=(0,1,2,3,4))
      par0 = params[0]
      par1 = params[1]
      par2 = params[2]
      par3 = params[3]
      par4 = params[4]
    if(arg == 'm'):
      gom = 'mc'
    narg = narg+1
  if opt == '-d':
    dd = arg
    if(arg == '20'):
      sub1 = 4
    if(arg == '30'):
      sub1 = 6
    if(arg == '50'):
      sub1 = 11
    if(arg == '70'):
      sub1 = 17
    narg = narg+1
  if opt == '-u':
    if(arg == '150'):
      sub2 = 37
      uu = "-150"
    if(arg == '200'):
      sub2 = 50
      uu = "-200"
    narg = narg+1
  if opt == '-s':
    stemp = arg
    narg = narg+1

if narg != 5:
  sys.exit("need -n -c -d -u -s")

nr_sub = sub2-sub1

print "Fitting A(k) as an inverse polynomial of order: ", npt

#================================ PARAMETERS ==================================
rstep = 4.0
rmax = 253.0
rmin = 4.0

r = np.arange(rmin,rmax,rstep)
nr = np.size(r)

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
#read in kernel J_0
kr, J0 = np.loadtxt('J0filt.dat', unpack=True)
#read in the mock file names
infile = "cmass_dr10v8.dat.xi"
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
figout = PdfPages('dinv'+dd+uu+'_poly'+str(npt)+'_'+gom+ss+'.pdf')
fout = open("dalpha"+dd+uu+"_poly"+str(npt)+"_"+gom+ss+".txt","w")
fouts = open("dstd"+dd+uu+"_poly"+str(npt)+"_"+gom+ss+".txt","w")
#==============================================================================

ksp = makeksp(0.15e-4,200.0,1000.0)
psp = readcamb("dr9_full.ztf", ns, 8.0, sigma8_lin, omega_b, omega_c, ksp)
print 'sigmaR =', sigmaR(ksp, psp, 8.0)
psp_c = psp*bias**2.
nk = np.size(ksp)

#spline no wiggle transfer function onto same k-grid as P(k), then apodize
tspnw = spllog(ksp, tnw, knw)
pc = apodize(tspnw, psp_c, ksp, ns, signl_ap)

kcovdex = np.squeeze(np.where(ksp < 1000.0/max(r)))
ksp_c = ksp[kcovdex]
pc = pc[kcovdex]

if(gom == 'gc'):
  npc = np.size(kcovdex)
  stream = calc_stream(ksp_c,npc,beta,par3,'exp')
  pvol = [par0,par1,par2]
  inttot = volint(pvol,pc,npc,stream,nbar,mesh)
  cov = cov00(ksp_c, inttot, kr, r) + par4
  print "Fitting using the Gaussian Covariance matrix..."
  covinv = np.linalg.inv(cov[sub1:sub2,sub1:sub2])
if(gom == 'mc'):
  cov = mcov00("red_norec.txt","red_norec_avg.xi") 
  print "Fitting using the mock Covariance matrix..."
  covinv = np.linalg.inv(cov[sub1:sub2,sub1:sub2])

xierr = np.zeros(sub2)
for i in range(0,sub2):
  xierr[i] = math.sqrt(cov[i,i])

xib0 = np.zeros((npt+1,nr))
for i in range(1,npt+1):
    xib0[i,:] = r**(i-3.)
pm = apodize(tspnw, psp, ksp, ns, snl)

#UNCOMMENT if want model with RSD
#sig_s = 3.0
#rsd = calc_stream(ksp, nk, beta, sig_s, 'exp')
#pm = pm*rsd

xi_mod = p2xi_lin(ksp, pm, kr, J0, r)
xib0[0,:] = xi_mod

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
ap = np.zeros((npt))
coeffs_out = np.zeros(npt+1)

bseed = 0.#math.log(bias**2.)

print "Now fitting...", infile

rin, data = np.loadtxt(infile, unpack=True, usecols=(0,1))

t_beg = time.clock()
for j in range(0,nalpha):
  model0 = np.zeros(nr_sub)
  spls0 = np.zeros(nr_sub)
  main0 = np.zeros(nr_sub)
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
chibf = funcchi2(opt)
dchi = chibf - chiopt[40]
for pp in range(0,npt):
  ap[pp] = coeffs_out[pp]
abf = opt[0]
b0bf = np.exp(opt[1])
ax = fig.add_subplot(222)
pyplot.plot(rin[sub1:sub2], rin[sub1:sub2]**2*data[sub1:sub2], 'x', markersize=4)
pyplot.plot(r[sub1:sub2], r[sub1:sub2]**2*model0)
pyplot.plot(r[sub1:sub2], r[sub1:sub2]**2*spls0, '*', markersize=4)
pyplot.plot(r[sub1:sub2], r[sub1:sub2]**2*main0, 'o', markersize=4)
pyplot.errorbar(rin[sub1:sub2], rin[sub1:sub2]**2*data[sub1:sub2],fmt=None, yerr=rin[sub1:sub2]**2*xierr[sub1:sub2])
pyplot.figtext(0.57,0.63,r'$b_0^2 = $'+str(b0bf))
for pp in range(0,npt):
  pyplot.figtext(0.57,0.63-0.02*(pp+1),r'$a_{%d} = {%f}$' % (pp+1, ap[pp]))
pyplot.figtext(0.57,0.63-0.02*(npt+1),r'$\alpha = $'+str(opt[0]))
pyplot.figtext(0.57,0.63-0.02*(npt+2),r'$\chi^2 = $'+str(chibf))
pyplot.title("Best fit")

chi1 = funcchi2([1.0,logb1])
b01 = b0opt[40]
for pp in range(0,npt):
  ap[pp] = coeffs_out[pp]
ax = fig.add_subplot(223)
pyplot.plot(rin[sub1:sub2], rin[sub1:sub2]**2*data[sub1:sub2], 'x', markersize=4)
pyplot.plot(r[sub1:sub2], r[sub1:sub2]**2*model0)
pyplot.plot(r[sub1:sub2], r[sub1:sub2]**2*spls0, '*', markersize=4)
pyplot.plot(r[sub1:sub2], r[sub1:sub2]**2*main0, 'o', markersize=4)
pyplot.errorbar(rin[sub1:sub2], rin[sub1:sub2]**2*data[sub1:sub2],fmt=None, yerr=rin[sub1:sub2]**2*xierr[sub1:sub2])
pyplot.figtext(0.15,0.18,r'$b_0^2 = $'+str(b01))
for pp in range(0,npt):
  pyplot.figtext(0.15,0.18-0.02*(pp+1),r'$a_{%d} = {%f}$' % (pp+1, ap[pp]))
pyplot.figtext(0.15,0.18-0.02*(npt+1),r'$\chi^2 = $'+str(chi1))
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
pat1 = np.sum(pnp[np.where(alpha_mini<1.0)])*dmini
print pnptot, pwptot, pwpltot, pwpitot
mnp = np.sum(alpha_mini*pnp)*dmini #mean alpha without prior
mwp = np.sum(alpha_mini*pwp)*dmini #mean alpha with prior
mwpl = np.sum(alpha_mini*pwpl)*dmini #mean alpha with ln prior
mwpi = np.sum(alpha_mini*pwpi)*dmini #mean alpha with inverse prior
sdnp = math.sqrt(np.sum((alpha_mini-mnp)**2*pnp)*dmini) #std without prior
sdwp = math.sqrt(np.sum((alpha_mini-mwp)**2*pwp)*dmini) #std with prior
sdwpl = math.sqrt(np.sum((alpha_mini-mwpl)**2*pwpl)*dmini) #std ln prior
sdwpi = math.sqrt(np.sum((alpha_mini-mwpi)**2*pwpi)*dmini) #std inv prior
da = mwpl - mwp
dsdanp = sdnp
dsda = sdwpl
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
fout.write(str(abf)+' '+str(chibf)+' '+str(b0bf)+'\n')
fouts.write(str(abf)+' '+str(dsda)+' '+str(dsdanp)+' '+str(pat1)+'\n')

figout.close()
fout.close()
fouts.close()
