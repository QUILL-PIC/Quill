#!/usr/bin/python

import numpy as np
import math
import os

__name__ = 'resread - provides functions for extractind data arrays\n\
from quill output files'
__doc__ = 'see source'

dx = 0
dy = 0
dz = 0
dt = 0
nx = 0
ny = 0
nz = 0
output_period = 0
n_ion_populations = 0
icmr = []
t_end = 0
tr_start = 0
deps = 0
deps_p = 0
deps_ph = 0
deps_i = 0
a0y = 0
a0z = 0
lmbda = 0
ne = 0
xsigma = 0
filmwidth = 0
nerflow = 0
Tlflow = 0
mcrlflow = 0
vlflow = 0
Trflow = 0
vrflow = 0

data_folder = '../results/'
t = '0'

v = 1

def xi( x, t ):
    return x - v * t

def read_parameters(log=None):
    'Reads nx, ny, etc. from *log*.'
    global dx,dy,dz,dt,nx,ny,nz,output_period,n_ion_populations,icmr,t_end,tr_start,\
    deps,deps_p,deps_ph,deps_i,a0y,a0z,lmbda,ne,xsigma,filmwidth,nerflow,\
    Tlflow, mcrlflow, vlflow, Trflow, vrflow
    if log is None:
	log = data_folder+'log'
    icmr = []
    f = open(log)
    for line in f:
	if line=='dx\n':
	    dx = float(f.next())
	elif line=='dy\n':
	    dy = float(f.next())
	elif line=='dz\n':
	    dz = float(f.next())
	elif line=='dt\n':
	    dt = float(f.next())
	elif line=='nx\n':
	    nx = int(f.next())
	elif line=='ny\n':
	    ny = int(f.next())
	elif line=='nz\n':
	    nz = int(f.next())
	elif line=='output_period\n':
	    output_period = float(f.next())
	elif line=='n_ion_populations\n':
	    n_ion_populations = int(f.next())
	elif line=='icmr\n':
	    icmr.append(float(f.next()))
	elif line=='t_end\n':
	    t_end = float(f.next())
	elif line=='tr_start\n':
	    tr_start = float(f.next())
	elif line=='deps\n':
	    deps = float(f.next())
	elif line=='deps_p\n':
	    deps_p = float(f.next())
	elif line=='deps_ph\n':
	    deps_ph = float(f.next())
	elif line=='deps_i\n':
	    deps_i = float(f.next())
	elif line=='a0y\n':
	    a0y = float(f.next())
	elif line=='a0z\n':
	    a0z = float(f.next())
	elif line=='lambda\n':
	    lmbda = float(f.next())
	elif line=='ne\n':
	    ne = float(f.next())
	elif line=='xsigma\n':
	    xsigma = float(f.next())
	elif line=='filmwidth\n':
	    filmwidth = float(f.next())
	elif line=='nerflow\n':
	    nerflow = float(f.next())
	elif line=='Tlflow\n':
	    Tlflow = float(f.next())
	elif line=='mcrlflow\n':
	    mcrlflow = float(f.next())
	elif line=='vlflow\n':
	    vlflow = float(f.next())
	elif line=='Trflow\n':
	    Trflow = float(f.next())
	elif line=='vrflow\n':
	    vrflow = float(f.next())
    f.close()

def density(name='rho',plane='xy'):
    'Returns 2d data for plane *plane* from file\n\
    data_folder+*name*+t.'
    f = open(data_folder+name+t)
    data = f.readlines()
    f.close()
    n = nx*ny + nx*nz + ny*nz
    density = np.empty(n)
    for i in np.arange(0,n,1):
	density[i] = float(data[i])
    if (plane!='xy') & (plane!='xz') & (plane!='yz'):
	print 'resread.density: warning: ambiguous value for *plane* -\n\
	'+plane+', value \'xy\' used instead'
	plane = 'xy'
    if plane=='xy':
	density = np.reshape(density[:-ny*nz],(nx,ny+nz))[:,:-nz]
    elif plane=='xz':
	density = np.reshape(density[:-ny*nz],(nx,ny+nz))[:,ny:]
    else:
	density = np.reshape(density[nx*(ny+nz):],(ny,nz))
    density = density.transpose()
    return density

def particles(name='phasespace',s=['x','y','g']):
    'Returns characteristics *s* for particles from the file\n\
    data_folder+*name*+t.'
    f = open(data_folder+name+t)
    data = f.readlines()
    f.close()
    n = len(data)/9
    m = len(s)
    a = np.empty((m,n))
    for i in np.arange(0,m,1):
	if s[i]=='q':
	    for j in np.arange(0,n,1):
		a[i][j] = float(data[9*j])
	elif s[i]=='x':
	    for j in np.arange(0,n,1):
		a[i][j] = float(data[9*j+1])
	elif s[i]=='y':
	    for j in np.arange(0,n,1):
		a[i][j] = float(data[9*j+2])
	elif s[i]=='z':
	    for j in np.arange(0,n,1):
		a[i][j] = float(data[9*j+3])
	elif s[i]=='ux':
	    for j in np.arange(0,n,1):
		a[i][j] = float(data[9*j+4])
	elif s[i]=='uy':
	    for j in np.arange(0,n,1):
		a[i][j] = float(data[9*j+5])
	elif s[i]=='uz':
	    for j in np.arange(0,n,1):
		a[i][j] = float(data[9*j+6])
	elif s[i]=='g':
	    for j in np.arange(0,n,1):
		a[i][j] = float(data[9*j+7])
	elif s[i]=='chi':
	    for j in np.arange(0,n,1):
		a[i][j] = float(data[9*j+8])
	elif s[i]=='t': # for qplot.tracks()
		for j in np.arange(n):
		    a[i][j] = tr_start + j*dt
	elif s[i]=='xi': # for qplot.tracks()
		for j in np.arange(n):
		    a[i][j] = xi( float(data[9*j+1]), ( tr_start + j*dt) )
	elif s[i]=='vx':
	    for j in np.arange(0,n,1):
		a[i][j] = float(data[9*j+4])/float(data[9*j+7])
	elif s[i]=='vy':
	    for j in np.arange(0,n,1):
		a[i][j] = float(data[9*j+5])/float(data[9*j+7])
	elif s[i]=='vz':
	    for j in np.arange(0,n,1):
		a[i][j] = float(data[9*j+6])/float(data[9*j+7])
	elif s[i]=='phi': # measured in xy plane countercloclwise from x-axis, lies in (-pi,pi]
	    for j in np.arange(0,n,1):
		a[i][j] = math.atan2(float(data[9*j+5]),float(data[9*j+4]))
	elif s[i]=='theta': # measured from xy-plane, lies in [-pi/2,pi/2]
	    for j in np.arange(0,n,1):
		a[i][j] = math.atan2(float(data[9*j+6]),np.sqrt(float(data[9*j+4])**2+float(data[9*j+5])**2))
	else:
	    print 'resread.particles: warning: ambiguous value for\n\
	    *s*['+str(i)+'] - '+s[i]+', value \'x\' used instead'
	    for j in np.arange(0,n,1):
		a[i][j] = float(data[9*j+1])
    return a

def t_data(name='energy',step=None):
    'Returns array of rows containing value of t and data from file\n\
    data_folder+*name*.'
    if step==None:
	step = dt
    f = open(data_folder+name)
    i = 0
    data = []
    for line in f:
	a = line.split('\t')
	data.append(i*step)
	i+=1
	for b in a:
	    data.append(float(b))
    data = np.reshape(data,(i,len(data)/i))
    return data

def tracks():
    'Returns a list of tracks from the data_folder. Each track is a dictionary with keys: t, x, y, z, ux, uy, uz, q, g, file'
    read_parameters()
    track_names = [x for x in os.listdir(data_folder) if x.startswith('track')]
    tracks = [read_track(x) for x in track_names]
    return tracks

def read_track(track_name):
    'Reads track from the specified track file. The returned track is a dictionary with keys: t, x, y, z, ux, uy, uz, q, g, file'
    raw_data = np.loadtxt(data_folder + track_name)
    raw_track = raw_data.reshape(9, -1, order='F')
    track_size = raw_track[0].size
    track = {'x' : raw_track[1],
             'y' : raw_track[2],
	     'z' : raw_track[3],
	     'ux' : raw_track[4],
	     'uy' : raw_track[5],
	     'uz' : raw_track[6],
	     'file' : data_folder + track_name,
	     'q' : raw_track[0],
	     'g' : raw_track[7],
	     't' : np.linspace(0, dt * track_size, track_size)}
    track['vx'] = track['ux'] / track['g']
    track['vy'] = track['uy'] / track['g']
    track['vz'] = track['uz'] / track['g']
    return track
