import matplotlib as mpl
mpl.use('Agg')
#
import os
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import airy
import resread
import qplot

pname = ''
punit = ''
pp_operation = ''
ppo_parameter = ''
cwd = './'

mode = 'w' # w - write, a - append

p0name = ''
p0unit = ''

def ppo(currentpvalue,conf,t,currentp0value=''):
    'Post-processing operations'
    global mode
    resread.data_folder = 'results/'
    resread.read_parameters()
    for i in np.arange(len(pp_operation)):
	#s = pname + '-' + pp_operation[i] + '--' + conf[1:] + '--' + t
	s1 = ''
	s2 = ''
	if currentp0value!='':
	    s1 = p0name + '-'
	    s2 = '-' + punit + '-'
	s3 = '-' + pp_operation[i] + '--' + conf[1:] + '--' + t
	s = cwd+'/'+s1+str(currentp0value)+s2+pname+'-'+str(currentpvalue)+'-'+punit+s3+'.png'
	write2file = False
	if pp_operation[i]=='energy':
	    tmp = resread.t_data('energy')
	    tmp = tmp[-1,1:]
	    if pname=='a0' and ppo_parameter[i]!=0:
		for j in np.arange(len(tmp)):
		    tmp[j] = ppo_parameter[i]*tmp[j]/currentpvalue**2
	    write2file = True
	elif pp_operation[i]=='density':
	    plt.clf()
	    qplot.density(ppo_parameter[i],'xy',data_folder='results/',save2=s)
	elif pp_operation[i]=='mollweide':
	    plt.clf()
	    qplot.mollweide(data_folder='results/',save2=s)
	    tmp = [qplot.directivity, qplot.directivity_lat,\
	    qplot.directivity_lng, qplot.directivity_lat1,\
	    qplot.directivity_lng1, qplot.directivity_lat2,\
	    qplot.directivity_lng2]
	    write2file = True
	elif pp_operation[i]=='ph_spectrum':
	    plt.clf()
	    qplot.spectrum(particles='g',data_folder='results/',save2=None)
	    sp = qplot.resread.t_data('spectrum_ph'+qplot.resread.t,qplot.resread.deps_ph)
	    for j in np.arange(len(sp[:,0])):
		if sp[j,1]>0:
		    sp[j,1] = np.log(sp[j,1])/np.log(10.)
		else:
		    sp[j,1] = 0
	    j = 0
	    epsmax = sp[1,0]
	    while True:
		if j>=len(sp[:,0]):
		    break
		if sp[j,1]==0:
		    break
		epsmax = sp[j,0]
		j+=1
	    plt.xlim([0,1.2*epsmax])
	    epsmin = epsmax/5.
	    n = 25
	    eps = np.linspace(epsmin,epsmax,n)
	    alpha = 1.
	    for k in np.arange(11):
		spapprox = np.zeros(n)
		Fperpapprox = np.sqrt(qplot.resread.a0y**2+qplot.resread.a0z**2)
		gammaapprox = np.sqrt(qplot.resread.a0y**2+qplot.resread.a0z**2)
		chi = 2*np.pi*3.862e-11/qplot.resread.lmbda*Fperpapprox*gammaapprox
		for j in np.arange(n):
		    kappa = (eps[j]/(alpha*chi*gammaapprox*0.511))**(2./3)
		    # airy(x) = [Ai(x), Ai'(x), Bi(x), Bi'(x)]
		    spapprox[j] = np.log( -np.exp(-2./3*kappa**1.5)/np.sqrt(2*np.pi)/(kappa+0.80049)**0.75 - 2./kappa*airy(kappa)[1] )/np.log(10.)
		chi = sp[int(eps[0]/qplot.resread.deps_ph),1] - spapprox[0]
		for j in np.arange(n):
		    spapprox[j] += chi
		diff = 0
		for j in np.arange(n):
		    diff += spapprox[j] - sp[int(eps[j]/qplot.resread.deps_ph),1]
		if diff>0:
		    alpha -= 1/2.**(k+1)
		elif diff<0:
		    alpha += 1/2.**(k+1)
	    tmp = [alpha,(eps[0]-eps[1])/(spapprox[1]-spapprox[0])]
	    write2file = True
	    plt.plot(eps,spapprox,'co')
	    plt.savefig(s)
	elif pp_operation[i]=='ions':
	    resread.t = '%g' % (resread.t_end - resread.dt) # str(resread.t_end-resread.dt)
	    iphsp = resread.particles('phasespace_'+str(resread.icmr[0])+'_',['q','g'])
	    gammaarray = [1.01,1.5,2.5,10]
	    Ni = 0
	    Niarray = np.zeros(len(gammaarray))
	    for j in np.arange(len(iphsp[0,:])):
		Ni += iphsp[0,j]
		for k in np.arange(len(gammaarray)):
		    if iphsp[1,j]>gammaarray[k]:
			Niarray[k] += iphsp[0,j]
	    for j in np.arange(len(gammaarray)):
		Niarray[j] = Niarray[j]/Ni
	    tmp = Niarray
	    write2file = True
	if write2file==True:
	    f = open(cwd+'/'+s1+pname+s3,mode)
	    tmp = np.insert(tmp,0,currentpvalue)
	    if currentp0value!='':
		tmp = np.insert(tmp,0,currentp0value)
	    s = ''
	    for j in np.arange(len(tmp)):
		s += str(tmp[j]) + '\t'
	    f.write(s[:-1]+'\n')
	    f.close()
    mode = 'a'

def prepare_conf(currentpvalue,use_p0=False):
    tmpf = open('conf')
    tmp = tmpf.readlines()
    tmpf.close()
    tmpf = open('conf','w')
    if use_p0==False:
	pn = pname
	pu = punit
    else:
	pn = p0name
	pu = p0unit
    for i in np.arange(0,len(tmp)-1,4):
	if tmp[i].strip()==pn:
	    tmpf.write(tmp[i])
	    tmpf.write(str(currentpvalue)+'\n')
	    if pu!='':
		tmpf.write(pu+'\n')
	    else:
		tmpf.write('#\n')
	    tmpf.write(tmp[i+3])
	elif pn=='a0' and tmp[i].strip()=='deps':
	    tmpf.write(tmp[i])
	    b = False
	    for j in np.arange(len(pp_operation)):
		if pp_operation[j]=='ph_spectrum':
		    #tmpf.write('%e'%(ppo_parameter[j]*currentpvalue**3)+'\n')
		    tmpf.write('%e'%(ppo_parameter[j]*currentpvalue**3*200/(200+currentpvalue))+'\n')
		    b = True
		    break
	    if b!=True:
		tmpf.write(tmp[i+1])
	    tmpf.write(tmp[i+2])
	    tmpf.write(tmp[i+3])
	elif pn=='a0' and tmp[i].strip()=='enthp_ph':
	    tmpf.write(tmp[i])
	    b = False
	    for j in np.arange(len(pp_operation)):
		if pp_operation[j]=='ph_spectrum':
		    tmpf.write('%d'%(5-4.*320/(currentpvalue+320))+'\n')
		    b = True
		    break
	    if b!=True:
		tmpf.write(tmp[i+1])
	    tmpf.write(tmp[i+2])
	    tmpf.write(tmp[i+3])
	else:
	    tmpf.write(tmp[i])
	    tmpf.write(tmp[i+1])
	    tmpf.write(tmp[i+2])
	    tmpf.write(tmp[i+3])
    tmpf.write('$')
    tmpf.close()

def geometric_progression(first,last,N):
    'Returns sequence of N values a[i] from first to last\n\
    (including) so that a[i+1]/a[i] = a[i]/a[i-1] for any value of i'
    a = np.empty(N)
    a[0] = first
    b = (float(last)/first)**(1./(N-1))
    for i in np.arange(1,N):
	a[i] = a[i-1]*b
    return a
