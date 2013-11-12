import matplotlib as mpl
mpl.use('Agg')
#
import os
import matplotlib.pyplot as plt
import numpy as np
#from scipy.special import airy
import resread
import qplot

p0_name = ''
p0_unit = ''
p1_name = ''
p1_unit = ''
pp_operation = ''
ppo_parameter = ''
cwd = './'
#mode = 'w' # w - write, a - append

def geometric_progression(first,last,N):
    'Returns sequence of N values a[i] from first to last\n\
    (including) so that a[i+1]/a[i] = a[i]/a[i-1] for any value of i'
    a = np.empty(N)
    a[0] = first
    if N!=1:
	b = (float(last)/first)**(1./(N-1))
    else:
	b = 1
    for i in np.arange(1,N):
	a[i] = a[i-1]*b
    return a

def p_substitution(pn,pv,pu=''):
    tmpf = open('conf')
    tmp = tmpf.readlines()
    tmpf.close()
    tmpf = open('conf','w')
    for i in np.arange(0,len(tmp)-1,4):
	if tmp[i].strip()==pn:
	    tmpf.write(tmp[i])
	    tmpf.write(str(pv)+'\n')
	    if pu!='':
		tmpf.write(pu+'\n')
	    else:
		tmpf.write('#\n')
	    tmpf.write(tmp[i+3])
	else:
	    tmpf.write(tmp[i])
	    tmpf.write(tmp[i+1])
	    tmpf.write(tmp[i+2])
	    tmpf.write(tmp[i+3])
    tmpf.write('$')
    tmpf.close()

def prepare_conf(pv0,pv1,config):
    p_substitution(p0_name,pv0,p0_unit)
    p_substitution(p1_name,pv1,p1_unit)
    # Exceptions
    if config=='.laser-piston':
	a0 = pv0*(p0_name=='a0') + pv1*(p1_name=='a0')
	ne = pv0*(p0_name=='ne') + pv1*(p1_name=='ne')
	if a0!=0:
	    p_substitution('t_end',np.ceil(22+1.e-7*a0**3))
	    if ne!=0:
		p_substitution('deps',2+1.7e-5*a0**2)

def ppo(pv0,pv1,config,t):
    'Post-processing operations'
    s1 = cwd+'/results-'+t+'/'
    s2 = p0_name+'-'+str(pv0)+p0_unit+'-'+p1_name+'-'+str(pv1)+p1_unit+'-'
    os.system('mkdir '+s1)
    a0 = pv0*(p0_name=='a0') + pv1*(p1_name=='a0')
    ne = pv0*(p0_name=='ne') + pv1*(p1_name=='ne')
    for operation in pp_operation:
	s = s1 + s2 + operation + '.png'
	resread.data_folder = 'results/'
	resread.read_parameters()
	if operation=='density':
	    plt.clf()
	    qplot.density(0.01*np.floor(resread.output_period*100*np.floor(resread.t_end/resread.output_period)),'xy',data_folder='results/',save2=s)
	    plt.clf()
	    qplot.density(0.01*np.floor(0.25*resread.output_period*100*np.floor(resread.t_end/resread.output_period)),'xy',data_folder='results/',save2=s[:-4]+'025.png')
	if operation=='spectrum':
	    plt.clf()
	    qplot.spectrum(0.01*np.floor(resread.output_period*100*np.floor(resread.t_end/resread.output_period)),data_folder='results/',save2=s)
	if operation=='tracks':
	    plt.clf()
	    qplot.tracks(['x','t'],'ie',colors='mg',data_folder='results/',save2=s)
	if operation=='i:x-ux':
	    plt.clf()
	    qplot.particles(0.01*np.floor(0.25*resread.output_period*100*np.floor(resread.t_end/resread.output_period)),['x','ux'],'i',colors='c',alpha=0.02,data_folder='results/',save2=None)
	    qplot.particles(0.01*np.floor(0.5*resread.output_period*100*np.floor(resread.t_end/resread.output_period)),['x','ux'],'i',colors='r',alpha=0.02,data_folder='results/',save2=None)
	    qplot.particles(0.01*np.floor(0.75*resread.output_period*100*np.floor(resread.t_end/resread.output_period)),['x','ux'],'i',colors='g',alpha=0.02,data_folder='results/',save2=None)
	    qplot.particles(0.01*np.floor(resread.output_period*100*np.floor(resread.t_end/resread.output_period)),['x','ux'],'i',colors='b',alpha=0.02,data_folder='results/',save2=s)
	if operation=='energy':
	    tmp = resread.t_data('energy')
	    plt.clf()
	    plt.plot(tmp[:,0],tmp[:,1],'k',tmp[:,0],tmp[:,2],'g',tmp[:,0],tmp[:,3],'r',tmp[:,0],tmp[:,4],'b',tmp[:,0],tmp[:,5],'m',tmp[:,0],tmp[:,1]+tmp[:,2]+tmp[:,3]+tmp[:,4]+tmp[:,5],'--k')
	    plt.xlabel('ct/lambda')
	    plt.ylabel('Energy, J')
	    plt.savefig(s)
	    energy0 = tmp[0,1]
	    tmp = tmp[-1,1:]
	    for i in np.arange(len(tmp)):
		tmp[i] = tmp[i]/energy0
	    f = open(s1+operation,ppo.energy2file)
	    ppo.energy2file = 'a' # function attribute, similar to 'static' in C++; 'a' - append, 'w' - write
	    tmp = np.insert(tmp,0,pv1)
	    tmp = np.insert(tmp,0,pv0)
	    s = ''
	    for a in tmp:
		s += str(a) + '\t'
	    f.write(s[:-1]+'\n')
	    f.close()

ppo.energy2file = 'w' # attribute must be initialized
