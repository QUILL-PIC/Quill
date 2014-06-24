#!/usr/bin/python

import matplotlib as mpl
mpl.use('Agg')
#
import matplotlib.pyplot as plt
import numpy as np
import os
import resread
import tcmap

__name__ =\
'qplot - provides functions for visualization of quill results'
__doc__ = 'see source'

def density(t=0,plane='xy',max_w=0,max_e_density=0,max_p_density=0,max_g_density=0,max_i_density=0,data_folder='../results/',axis=[],save2=''):
    'Plots w and particle densities.'
    resread.t = '%g' % t
    resread.data_folder = data_folder
    resread.read_parameters()
    #
    tcmap.red()
    tcmap.green()
    tcmap.blue()
    tcmap.orange()
    tcmap.purple()
    #
    plt.title('Particle densities')
    plt.xlabel(plane[0]+' (wavelengths)')
    plt.ylabel(plane[1]+' (wavelengths)')
    if plane=='xz':
	xlength = resread.nx*resread.dx
	ylength = resread.nz*resread.dz
    elif plane=='yz':
	xlength = resread.ny*resread.dy
	ylength = resread.nz*resread.dz
    else:
	xlength = resread.nx*resread.dx
	ylength = resread.ny*resread.dy
    if axis!=[]:
	plt.axis(axis)
    #
    edensity = np.multiply(resread.density('rho',plane),-1)
    w = resread.density('w',plane)
    if max_e_density==0:
	max_e_density = edensity.max()
	print 'qplot.density: max_e_density = ', max_e_density
    if max_w==0:
	max_w = w.max()
	print 'qplot.density: max_w = ', max_w
    if max_p_density==0:
	max_p_density = max_e_density
    if max_g_density==0:
	max_g_density = max_e_density
    if max_i_density==0:
	max_i_density = max_e_density
    #
    plt.imshow(resread.density('rho_ph',plane),'tcmap_blue',interpolation='none',vmin=0,vmax=max_g_density,origin='lower',extent=(0,xlength,0,ylength))
    plt.imshow(w,'tcmap_orange',interpolation='none',vmin=0,vmax=max_w,origin='lower',extent=(0,xlength,0,ylength))
    if resread.icmr!=[]:
	plt.imshow(resread.density('irho_'+str(resread.icmr[0])+'_',plane),'tcmap_purple',interpolation='none',vmin=0,vmax=max_i_density,origin='lower',extent=(0,xlength,0,ylength))
    plt.imshow(edensity,'tcmap_green',interpolation='none',vmin=0,vmax=max_e_density,origin='lower',extent=(0,xlength,0,ylength))
    plt.imshow(resread.density('rho_p',plane),'tcmap_red',interpolation='none',vmin=0,vmax=max_p_density,origin='lower',extent=(0,xlength,0,ylength))
    #
    if save2=='':
	plt.show()
    else:
	plt.savefig(save2)

def particles(t=0,space=['x','y'],particles='geip',colors='bgmrcyk',r=5,alpha=0.01,cmap='jet',gamma=0,data_folder='../results/',axis=[],save2=''):
    'Plots particles as dots in (phase)*space*.\n\
    \n\
    Examples:\n\
    qplot.particles(10,[\'x\',\'y\']),\n\
    qplot.particles(15,[\'x\',\'y\',\'g\'],\'i\',gamma=1).'
    resread.t = '%g' % t
    resread.data_folder = data_folder
    resread.read_parameters()
    plt.xlabel(space[0])
    plt.ylabel(space[1])
    if axis!=[]:
	plt.axis(axis)
    lp = list(particles)
    s = []
    for p in lp:
	if p=='e':
	    s.append('')
	elif p=='p':
	    s.append('_p')
	elif p=='g':
	    s.append('_ph')
	elif p=='i':
	    for cmr in resread.icmr:
		s.append('_'+str(cmr)+'_')
    if len(space)==2:
	c = list(colors)
	i=0
	for suffix in s:
	    phspace = resread.particles('phasespace'+suffix,space)
	    plt.scatter(phspace[0,:],phspace[1,:],s=r*r,color=c[i],alpha=alpha,edgecolors='None')
	    i+=1
    elif len(space)==3:
	cmap = tcmap.get(cmap,gamma=gamma)
	plt.title(space[2])
	phspace = resread.particles('phasespace'+s[0],space)
	plt.scatter(phspace[0,:],phspace[1,:],s=r*r,c=phspace[2,:],cmap=cmap,edgecolors='None')
	plt.colorbar()
    if save2=='':
	plt.show()
    elif save2!=None:
	plt.savefig(save2)

def tracks(space=['x','y'],particles='geip',t0=0,t1=0,colors='bgmrcyk',cmaps=['jet'],clims2all=0,axis=[],save2='',r=2,data_folder='../results/'):
    'Plots particle tracks as lines in 2D or dots in 3D (phase)*space*\n\
    at [tr_start+*t0*,tr_start+*t1*]\n\
    Examples:\n\
    tracks(space=[\'t\',\'x\',\'g\'],particles=\'e\',clims2all=1)\n\
    tracks(colors=\'\',axis=[9,11,8,12],save2=\'/tmp/tmp.png\') # plots\
    each track with its own color and saves picture to a file\n\
    '
    # clims2all=1 -> color limits for cmaps[i] (vmin,vmax) the same
    # for all particles (trajectories) of this type
    resread.data_folder = data_folder
    resread.read_parameters()
    track_names = []
    for filename in os.listdir(data_folder):
	if filename.find('track')==0:
	    track_names.append(filename)
    tracks = []
    cmr = []
    for trackname in track_names:
	resread.t = ''
	tmp = resread.particles(trackname,space)
	if t1!=0 and np.floor(t1/resread.dt)<len(tmp[0,:]):
	    tracks.append(tmp[:,np.floor(t0/resread.dt):np.floor(t1/resread.dt)])
	else:
	    tracks.append(tmp[:,9*np.floor(t0/resread.dt):])
	cmr.append(float( trackname[ trackname.find('_')+1 : trackname.find('_',trackname.find('_')+1) ] ))
    def cmr2int(a):
	'Converts cmr *a* to the int index of *colors* or *cmaps*'
	if len(space)==2:
	    if colors=='':
		return 0
	    else:
		maxv = len(colors)
	elif len(space)==3:
	    maxv = len(cmaps)
	if particles.find('i')!=-1:
	    s = particles[:particles.find('i')] + (len(resread.icmr)-1)*'i' + particles[particles.find('i'):]
	else:
	    s = particles
	if a==-1 and s.find('e')!=-1:
	    return s.find('e')%maxv
	elif a==1 and s.find('p')!=-1:
	    return s.find('p')%maxv
	elif a==0 and s.find('g')!=-1:
	    return s.find('g')%maxv
	elif s.find('i')!=-1:
	    for j in np.arange(len(resread.icmr)):
		if a==resread.icmr[j]:
		    return (particles.find('i')+j)%maxv
	else:
	    return 0
    if len(space)==2:
	if colors!='':
	    for i in np.arange(len(track_names)):
		if cmr[i]==-1 and particles.find('e')!=-1:
		    plt.plot(tracks[i][0,:],tracks[i][1,:],colors[cmr2int(cmr[i])])
		elif cmr[i]==1 and particles.find('p')!=-1:
		    plt.plot(tracks[i][0,:],tracks[i][1,:],colors[cmr2int(cmr[i])])
		elif cmr[i]==0 and particles.find('g')!=-1:
		    plt.plot(tracks[i][0,:],tracks[i][1,:],colors[cmr2int(cmr[i])])
		elif particles.find('i')!=-1:
		    for j in np.arange(len(resread.icmr)):
			if cmr[i]==resread.icmr[j]:
			    plt.plot(tracks[i][0,:],tracks[i][1,:],colors[cmr2int(cmr[i])])
	else:
	    for i in np.arange(len(track_names)):
		if cmr[i]==-1 and particles.find('e')!=-1:
		    plt.plot(tracks[i][0,:],tracks[i][1,:])
		elif cmr[i]==1 and particles.find('p')!=-1:
		    plt.plot(tracks[i][0,:],tracks[i][1,:])
		elif cmr[i]==0 and particles.find('g')!=-1:
		    plt.plot(tracks[i][0,:],tracks[i][1,:])
		elif particles.find('i')!=-1:
		    for j in np.arange(len(resread.icmr)):
			if cmr[i]==resread.icmr[j]:
			    plt.plot(tracks[i][0,:],tracks[i][1,:])
    elif len(space)==3:
	plt.title(space[2])
	if clims2all==1:
	    def pname(a):
		if a==-1:
		    return 'e'
		elif a==1:
		    return 'p'
		elif a==0:
		    return 'g'
		else:
		    return 'i'
	    for m in [-1,1,0]+resread.icmr:
		if particles.find(pname(m))!=-1:
		    tmp0 = []
		    tmp1 = []
		    tmp2 = []
		    for i in np.arange(len(track_names)):
			if cmr[i]==m:
			    for j in np.arange(len(tracks[i][0,:])):
				tmp0.append(tracks[i][0,j])
				tmp1.append(tracks[i][1,j])
				tmp2.append(tracks[i][2,j])
		    plt.scatter(tmp0,tmp1,c=tmp2,cmap=cmaps[cmr2int(m)],s=r*r,faceted=False)
	else:
	    for i in np.arange(len(track_names)):
		tmp0 = []
		tmp1 = []
		tmp2 = []
		for j in np.arange(len(tracks[i][0,:])):
		    tmp0.append(tracks[i][0,j])
		    tmp1.append(tracks[i][1,j])
		    tmp2.append(tracks[i][2,j])
		if particles.find('e')!=-1 and cmr[i]==-1:
		    plt.scatter(tmp0,tmp1,c=tmp2,cmap=cmaps[cmr2int(cmr[i])],s=r*r,faceted=False)
		if particles.find('p')!=-1 and cmr[i]==1:
		    plt.scatter(tmp0,tmp1,c=tmp2,cmap=cmaps[cmr2int(cmr[i])],s=r*r,faceted=False)
		if particles.find('g')!=-1 and cmr[i]==0:
		    plt.scatter(tmp0,tmp1,c=tmp2,cmap=cmaps[cmr2int(cmr[i])],s=r*r,faceted=False)
		if particles.find('i')!=-1:
		    for k in np.arange(len(resread.icmr)):
			if cmr[i]==resread.icmr[k]:
			    plt.scatter(tmp0,tmp1,c=tmp2,cmap=cmaps[cmr2int(cmr[i])],s=r*r,faceted=False)
    else:
	print 'qplot.tracks: warning: ambiguous value for *space*'
    if axis!=[]:
	plt.axis(axis)
    plt.xlabel(space[0])
    plt.ylabel(space[1])
    if save2=='':
	plt.show()
    else:
	plt.savefig(save2)

def rpattern(t=None,particles='geip',colors='bgmrcyk',dphi=0.1,save2='',data_folder='../results/'):
    'Plots radiation pattern of the emitted energy\n\
    Examples:\n\
    rpattern() # plots radiation patterns for all particles\n\
    at t_end'
    resread.data_folder = data_folder
    resread.read_parameters()
    if t==None:
	t = resread.t_end - resread.dt
    resread.t = '%g' % t
    lp = list(particles)
    s = []
    for p in lp:
	if p=='e':
	    s.append('')
	elif p=='p':
	    s.append('_p')
	elif p=='g':
	    s.append('_ph')
	elif p=='i':
	    for cmr in resread.icmr:
		s.append('_'+str(cmr)+'_')
    c = list(colors)
    ci = 0
    f = plt.figure()
    ax = f.add_subplot(111,polar=True)
    for suffix in s:
	qgphi = resread.particles('phasespace'+suffix,['q','g','phi'])
	phi = np.arange(-np.pi,np.pi+dphi,dphi)
	n = len(phi)
	rp = np.zeros(n)
	for i in np.arange(len(qgphi[0,:])):
	    j = np.floor((qgphi[2,i]+np.pi)/dphi)
	    rp[j] += np.fabs(qgphi[0,i])*qgphi[1,i]
	rp[0] += rp[n-1]
	rp[n-1] = rp[0]
	phi[n-1] = phi[0]
	rpint = 0
	for i in np.arange(n):
	    rpint += rp[i]
	if rpint!=0:
	    for i in np.arange(n):
		rp[i] = rp[i]*2*np.pi/(dphi*rpint)
	ax.plot(phi,rp,color=c[ci])
	ci+=1
    #ax.set_theta_direction(-1)
    #ax.set_theta_offset(np.pi/2)
    if save2=='':
	plt.show()
    else:
	plt.savefig(save2)

def spectrum(t=None,particles='geip',colors='bgmrcyk',sptype='simple',axis=[],save2='',data_folder='../results/'):
    'spectrum() # plots spectrum for all particles\n\
    at t_end\n\
    Examples:\n\
    spectrum(10,sptype=\'loglog\') # energy distr.\n\
    in log-log axes'
    plt.xlabel('kinetic energy, MeV')
    plt.ylabel('dN/deps, a.u.')
    if axis!=[]:
	plt.axis(axis)
    resread.data_folder = data_folder
    resread.read_parameters()
    if t==None:
	t = resread.t_end - resread.dt
    resread.t = '%g' % t
    s = []
    deps = []
    ci = []
    i = -1
    for p in list(particles):
	i+=1
	if p=='e':
	    s.append('')
	    deps.append(resread.deps)
	    ci.append(i%len(colors))
	elif p=='p':
	    s.append('_p')
	    deps.append(resread.deps_p)
	    ci.append(i%len(colors))
	elif p=='g':
	    s.append('_ph')
	    deps.append(resread.deps_ph)
	    ci.append(i%len(colors))
	elif p=='i':
	    for cmr in resread.icmr:
		s.append('_'+str(cmr)+'_')
		deps.append(resread.deps_i)
		ci.append(i%len(colors))
    for i in np.arange(len(s)):
	sp = resread.t_data('spectrum'+s[i]+resread.t,deps[i])
	if sptype=='energy':
	    for j in np.arange(len(sp[:,0])):
		sp[j,1] = sp[j,0]*sp[j,1]
	elif sptype=='loglog':
	    for j in np.arange(len(sp[:,0])):
		sp[j,1] = sp[j,0]*sp[j,1]
		if sp[j,1]>0:
		    sp[j,1] = np.log10( sp[j,1] )
		else:
		    sp[j,1] = 0
		sp[j,0] = np.log10( sp[j,0] )
	for j in np.arange(len(sp[:,0])):
	    if sp[j,1]>0:
		sp[j,1] = np.log(sp[j,1])/np.log(10.)
	    else:
		sp[j,1] = 0
	plt.plot(sp[:,0],sp[:,1],colors[ci[i]])
    if save2!=None:
	if save2=='':
	    plt.show()
	else:
	    plt.savefig(save2)

directivity = 0
directivity_lat = 0
directivity_lng = 0
directivity_lat1 = 0
directivity_lng1 = 0
directivity_lat2 = 0
directivity_lng2 = 0
def mollweide(t=None,nlongitude=80,nlatitude=40,Nlevels=15,save2='',data_folder='../results/'):
    'Plots photon radiation pattern in Mollweide projection (z-axis\n\
    sticks out of the north pole and y-axis sticks out of the\n\
    projection center) and computes (antenna-like) directivity.\n\
    Examples:\n\
    ...'
    global directivity, directivity_lat, directivity_lng,\
    directivity_lat1, directivity_lng1, directivity_lat2,\
    directivity_lng2
    resread.data_folder = data_folder
    resread.read_parameters()
    if t==None:
	t = resread.t_end - resread.dt
    resread.t = '%g' % t
    #
    lat = np.linspace(-np.pi/2,np.pi/2,nlatitude)
    lng = np.linspace(-np.pi,np.pi,nlongitude)
    rp = np.zeros((nlatitude,nlongitude))
    #
    qgthphi = resread.particles('phasespace_ph',['q','g','theta','phi'])
    for i in np.arange(len(qgthphi[0,:])):
	# theta - latitude, from xy-plane, [-pi/2,pi/2]
	# phi - longitude, in xy-plane, (-pi,pi]
	y1 = (qgthphi[2,i]+np.pi/2)*(nlatitude-1)/np.pi
	j = np.floor(y1)
	y1 = y1 - j
	y2 = 1 - y1
	x1 = (qgthphi[3,i]+np.pi)*(nlongitude-1)/(2*np.pi)
	k = np.ceil(x1)
	x1 = k - x1
	x2 = 1 - x1
	if j!=nlatitude:
	    rp[j,k] += y1*x1*qgthphi[0,i]*qgthphi[1,i]
	    rp[j+1,k] += y2*x1*qgthphi[0,i]*qgthphi[1,i]
	    rp[j+1,k-1] += y2*x2*qgthphi[0,i]*qgthphi[1,i]
	    rp[j,k-1] += y1*x2*qgthphi[0,i]*qgthphi[1,i]
	else:
	    rp[j,k] += y1*x1*qgthphi[0,i]*qgthphi[1,i]
	    rp[j,k-1] += y1*x2*qgthphi[0,i]*qgthphi[1,i]
    tmp1 = 0
    tmp2 = 0
    for i in np.arange(nlongitude):
	tmp1 += rp[0,i]
	tmp2 += rp[nlatitude-1,i]
    for i in  np.arange(nlongitude):
	rp[0,i] = tmp1
	rp[nlatitude-1,i] = tmp2
    for i in np.arange(1,nlatitude-1):
	rp[i,0] += rp[i,nlongitude-1]
	rp[i,nlongitude-1] = rp[i,0]
    for j in np.arange(nlongitude):
	for i in np.arange(1,nlatitude-1):
	    rp[i,j] = rp[i,j]/np.cos(lat[i])
	rp[0,j] = rp[0,j]*8/np.cos(lat[1])/(nlongitude-1)
	rp[nlatitude-1,j] = rp[nlatitude-1,j]*8/np.cos(lat[1])/(nlongitude-1)
    integral = 0
    integral1 = 0
    for j in np.arange(nlongitude):
	for i in np.arange(1,nlatitude-1):
	    integral += rp[i,j]*np.cos(lat[i])
	    integral1 += np.cos(lat[i])
    integral += rp[0,0]*np.cos(lat[1])/8
    integral += rp[nlatitude-1,0]*np.cos(lat[1])/8
    integral1 += np.cos(lat[1])/8
    integral1 += np.cos(lat[1])/8
    for j in np.arange(nlongitude):
	for i in np.arange(nlatitude):
	    rp[i,j] = rp[i,j]*integral1/integral
    directivity = rp.max()
    directivity_lat = lat[np.floor(rp.argmax()/nlongitude)]
    directivity_lng = lng[rp.argmax()%nlongitude]
    directivity_lat1 = lat[np.floor(rp[:,:nlongitude/2].argmax()/(nlongitude/2))]
    directivity_lng1 = lng[rp[:,:nlongitude/2].argmax()%(nlongitude/2)]
    directivity_lat2 = lat[np.floor(rp[:,nlongitude/2:].argmax()/(nlongitude-nlongitude/2))]
    directivity_lng2 = lng[nlongitude/2 + rp[:,nlongitude/2:].argmax()%(nlongitude/2)]
    f = plt.figure()
    #ax = f.add_subplot(111) # qwe
    #ax.imshow(rp,interpolation='none') # qwe
    ax = f.add_subplot(111,projection='mollweide')
    ax.contour(lng,lat,rp,Nlevels,cmap='jet',origin=None)
    #
    if save2!=None:
	if save2=='':
	    plt.show()
	else:
	    plt.savefig(save2)

def field(t=0,field='ex',plane='xy',fmax=None,data_folder='../results/',axis=[],save2=''):
    'Plots fields.'
    resread.t = '%g' % t
    resread.data_folder = data_folder
    resread.read_parameters()
    #
    plt.title(field)
    plt.xlabel(plane[0]+' (wavelengths)')
    plt.ylabel(plane[1]+' (wavelengths)')
    if plane=='xz':
	xlength = resread.nx*resread.dx
	ylength = resread.nz*resread.dz
    elif plane=='yz':
	xlength = resread.ny*resread.dy
	ylength = resread.nz*resread.dz
    else:
	xlength = resread.nx*resread.dx
	ylength = resread.ny*resread.dy
    if axis!=[]:
	plt.axis(axis)
    #
    f = resread.density(field,plane)
    if fmax==None:
	fmax = np.max( [-np.min(f), np.max(f)] )
	print 'fmax = ', fmax
    #
    plt.imshow(f,'bwr',interpolation='none',vmin=-fmax,vmax=fmax,origin='lower',extent=(0,xlength,0,ylength))
    #
    if save2=='':
	plt.show()
    else:
	plt.savefig(save2)

