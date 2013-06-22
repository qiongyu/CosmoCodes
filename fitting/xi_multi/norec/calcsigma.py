import math
import numpy as np
import matplotlib
from matplotlib import pyplot
from matplotlib.backends.backend_pdf import PdfPages

a = np.loadtxt('grids.txt',dtype=({'names': ['infile'],'formats': ['S100']}))
infile = a['infile']
nf = np.size(infile)
malpha = np.zeros(nf)
mepsilon = np.zeros(nf)
calpha = np.zeros(nf)
cepsilon = np.zeros(nf)
cae = np.zeros(nf)

na = 241
dex = np.arange(0,na)
pa = np.zeros(na)
ca = np.zeros(na)
alpha = np.zeros(na)
ca_ez = np.zeros(na)

fout = open("stdevs.dat","w")
figout = PdfPages("chi2.pdf")

matplotlib.rcParams['font.size']=7
font0 = matplotlib.font_manager.FontProperties()
font0.set_size(10)
fig = pyplot.figure()

fig = pyplot.figure()
for i in range(0,nf):
  ag, eg, cg = np.loadtxt(infile[i],unpack=True)
  tophat = np.squeeze(np.where( (eg < 0.15) & (eg > -0.15) ))
  ag = ag[tophat]
  eg = eg[tophat]
  cg = cg[tophat]

  ne = np.size(tophat)/na
  ce = np.zeros(ne)
  pe = np.zeros(ne)
  epsilon = np.zeros(ne)
  pae = np.zeros(na*ne)

  cg = cg+ (np.log(ag)/0.15)**2.
  pg = np.exp(-0.5*cg)

  for j in range(0,na): #at each alpha...
    ld = j*ne
    hd = (j+1)*ne
    temp_pa = pg[ld:hd] #p(alpha) = sum[over epsilon] p(alpha,epsilon)
    pa[j] = np.sum(temp_pa)
    ca[j] = -2.*np.log(pa[j])
    alpha[j] = ag[ld]
    ca_ez[j] = cg[ld+ne/2]
  for j in range(0,ne): #at each epsilon...
    d = dex*ne+j
    temp_pe = pg[d]
    pe[j] = np.sum(temp_pe)
    ce[j] = -2.*np.log(pe[j])
    epsilon[j] = eg[d[0]]

  norma = np.sum(pa)
  pa = pa/norma

  malpha[i] = np.sum(pa*alpha)
  calpha[i] = np.sum(pa*(alpha-malpha[i])**2)

  norme = np.sum(pe)
  pe = pe/norme

  mepsilon[i] = np.sum(pe*epsilon)
  cepsilon[i] = np.sum(pe*(epsilon-mepsilon[i])**2)

  for j in range(0,na*ne):
    pae[j] = pg[j]*(ag[j]-malpha[i])*(eg[j]-mepsilon[i])
  normp = np.sum(pg)
  cae[i] = np.sum(pae)/normp

  rho = cae[i]/(math.sqrt(calpha[i]*cepsilon[i]))
  
  for j in range(0,na*ne):
    pae[j] = pg[j]*(ag[j]-malpha[i])*(ag[j]-malpha[i]) 
  testa = np.sum(pae)/normp
  for j in range(0,na*ne):
    pae[j] = pg[j]*(eg[j]-mepsilon[i])*(eg[j]-mepsilon[i]) 
  teste = np.sum(pae)/normp
  print testa, teste

  print malpha[i], calpha[i], mepsilon[i], cepsilon[i], cae[i], rho

  fout.write(str(malpha[i])+' '+str(calpha[i])+' '+str(mepsilon[i])+' '+str(cepsilon[i])+' '+str(cae[i])+' '+str(rho)+'\n')

  ax = fig.add_subplot(221)
  pyplot.plot(alpha, ca-ca.min(), 'k-')
  pyplot.plot(alpha, ca_ez-ca_ez.min(), 'k--')
  pyplot.xlabel(r"$\alpha$", fontproperties=font0)
  pyplot.ylabel(r"$\Delta\chi^2(\alpha)$", fontproperties=font0)
  pyplot.figtext(0.34,0.86,r'Min $\chi^2=%4.2f$'%ca.min(), fontproperties=font0)

  ax = fig.add_subplot(222)
  pyplot.plot(epsilon, ce-ce.min(), 'k-')
  pyplot.xlabel(r"$\epsilon$", fontproperties=font0)
  pyplot.ylabel(r"$\Delta\chi^2(\epsilon)$", fontproperties=font0)
  pyplot.xlim(-0.3,0.3)
  pyplot.figtext(0.76,0.86,r'Min $\chi^2=%4.2f$'%ce.min(), fontproperties=font0)

  ax = fig.add_subplot(223)
  pyplot.plot(alpha, pa, 'k-')
  pyplot.xlabel(r"$\alpha$", fontproperties=font0)
  pyplot.ylabel(r"$p(\alpha)$", fontproperties=font0)
  pyplot.figtext(0.32,0.42,r'$\langle\alpha\rangle=%4.3f\pm%4.3f$'%(malpha[i],math.sqrt(calpha[i])), fontproperties=font0)

  ax = fig.add_subplot(224)
  pyplot.plot(epsilon, pe, 'k-')
  pyplot.xlabel(r"$\epsilon$", fontproperties=font0)
  pyplot.ylabel(r"$p(\epsilon)$", fontproperties=font0)
  pyplot.xlim(-0.3,0.3)
  pyplot.figtext(0.745,0.42,r'$\langle\epsilon\rangle=%4.3f\pm%4.3f$'%(mepsilon[i],math.sqrt(cepsilon[i])), fontproperties=font0)
  pyplot.figtext(0.78,0.39,r'$C_{\alpha\epsilon}=%5.4f$'%(cae[i]), fontproperties=font0)
  pyplot.figtext(0.784,0.36,r'$\rho_{\alpha\epsilon}=%4.2f$'%(rho), fontproperties=font0)
  figout.savefig()
  pyplot.clf()

fout.close()
figout.close()
