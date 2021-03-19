#!/usr/bin/python

import numpy as np
import sys
import os
import expression_parser

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
nfilm = 0
filmwidth = 0
nerflow = 0
Tlflow = 0
mcrlflow = 0
vlflow = 0
Trflow = 0
vrflow = 0
catching = False
dump_photons = False
particles_for_output = 'e'
output_mode = 0

data_folder = '../results/'
t = '0'

v = 1

def xi( x, t ):
    return x - v * t

def read_parameters(log=None):
    'Reads nx, ny, etc. from *log*.'
    global dx,dy,dz,dt,nx,ny,nz,output_period,n_ion_populations,icmr,t_end,tr_start,\
    deps,deps_p,deps_ph,deps_i,a0y,a0z,lmbda,ne,xsigma,nfilm,filmwidth,nerflow,\
    Tlflow, mcrlflow, vlflow, Trflow, vrflow, catching, dump_photons, particles_for_output, output_mode
    if log is None:
        log = os.path.join(data_folder,'log')
    icmr = []
    reset_globals()
    f = open(log)
    for line in f:
        if line=='dx\n':
            dx = float(next(f))
        elif line=='dy\n':
            dy = float(next(f))
        elif line=='dz\n':
            dz = float(next(f))
        elif line=='dt\n':
            dt = float(next(f))
        elif line=='nx\n':
            nx = int(next(f))
        elif line=='ny\n':
            ny = int(next(f))
        elif line=='nz\n':
            nz = int(next(f))
        elif line=='output_period\n':
            output_period = float(next(f))
        elif line=='n_ion_populations\n':
            n_ion_populations = int(next(f))
        elif line=='icmr\n':
            icmr.append(float(next(f)))
        elif line=='t_end\n':
            t_end = float(next(f))
        elif line=='tr_start\n':
            tr_start = float(next(f))
        elif line=='deps\n':
            deps = float(next(f))
        elif line=='deps_p\n':
            deps_p = float(next(f))
        elif line=='deps_ph\n':
            deps_ph = float(next(f))
        elif line=='deps_i\n':
            deps_i = float(next(f))
        elif line=='a0y\n':
            a0y = float(next(f))
        elif line=='a0z\n':
            a0z = float(next(f))
        elif line=='lambda\n':
            lmbda = float(next(f))
        elif line=='ne\n':
            ne = float(next(f))
        elif line=='xsigma\n':
            xsigma = float(next(f))
        elif line=='nfilm\n':
            nfilm = float(next(f))
        elif line=='filmwidth\n':
            filmwidth = float(next(f))
        elif line=='nerflow\n':
            nerflow = float(next(f))
        elif line=='Tlflow\n':
            Tlflow = float(next(f))
        elif line=='mcrlflow\n':
            mcrlflow = float(next(f))
        elif line=='vlflow\n':
            vlflow = float(next(f))
        elif line=='Trflow\n':
            Trflow = float(next(f))
        elif line=='vrflow\n':
            vrflow = float(next(f))
        elif line.strip() == 'catching':
            ss = next(f).strip()
            if ss == 'on':
                catching = True
        elif line.strip() == 'dump_photons':
            ss = next(f).strip()
            if ss == 'on':
                dump_photons = True
        elif line.strip() == 'particles_for_output':
            particles_for_output = next(f).strip().replace('ph','g')
        elif line.strip() == 'output_mode':
            output_mode = int(next(f))
    f.close()

def density(name='rho',plane='xy', log=None):
    'Returns 2d data for plane *plane* from file\n\
    data_folder+*name*+t.'
    filename = os.path.join(data_folder, name + t)
    if output_mode == 1 and name != 'w' and name != 'inv':
        sys.path.append(os.path.join('..', 'chameleon'))
        try:
            import chameleon
        except Exception as e:
            raise ImportError('Chameleon cannot be imported. Check if it is compiled. Error: ' + str(e))
        chameleon.configure(log if log is not None else os.path.join(data_folder, 'log'))
        return chameleon.read2d(filename, plane)
    else:
        with open(filename) as f:
            density = np.array([float(x) for x in f])
        n = nx*ny + nx*nz + ny*nz
        if density.size != n:
            raise Exception("The size of data in [%s] equal to %d doesn't match n=%d" % (name + t, density.size, n))
        if (plane!='xy') & (plane!='xz') & (plane!='yz'):
            print('resread.density: warning: ambiguous value for *plane* - {0}, value \'xy\' used instead'.format(plane))
            plane = 'xy'
        if plane=='xy':
            density = np.reshape(density[:-ny*nz],(nx,ny+nz))[:,:-nz]
        elif plane=='xz':
            density = np.reshape(density[:-ny*nz],(nx,ny+nz))[:,ny:]
        else:
            density = np.reshape(density[nx*(ny+nz):],(ny,nz))
        density = density.transpose()
        return density

def particles(name='phasespace', s=['x','y','g'], every=1):
    'Returns characteristics *s* for particles from the file\n\
    data_folder+*name*+t.'
    f = open(os.path.join(data_folder, name+t))
    data = f.readlines() if every == 1 else [line for i, line in enumerate(f.readlines()) if (i//9) % every == 0]
    f.close()
    n = len(data)//9
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
        # phi = 0, theta = 0 - direction of x axis
        # theta = pi / 2 - direction of z axis
        # phi = pi / 2, theta = 0 - direction of y axis
        elif s[i]=='phi': 
            a[i,:] = np.arctan2(np.asfarray(data[5::9]), np.asfarray(data[4::9]))
        elif s[i]=='theta':
            x = np.asfarray(data[4::9])
            y = np.asfarray(data[5::9])
            a[i,:] = np.arctan2(np.asfarray(data[6::9]), np.sqrt(x * x + y * y))
        else:
            print('resread.particles: warning: ambiguous value for s[{0}] - {1}, value \'x\' used instead'.format(i, s[i]))
            for j in np.arange(0,n,1):
                a[i][j] = float(data[9*j+1])
    return a

def t_data(name='energy', step=None, silent=False):
    'Returns array of rows containing value of t and data from file\n\
    data_folder+*name*.'
    if not silent:
        print ('Fetching t_data from file: {0}; data_folder = {1}'.format(name, data_folder))
    if step==None:
        step = dt
    f = open(os.path.join(data_folder, name))
    i = 0
    data = []
    for line in f:
        a = line.split('\t')
        data.append(i*step)
        i+=1
        for b in a:
            data.append(float(b))
    if i == 0:
        raise ValueError('The file is empty!')
    data = np.reshape(data,(i,len(data)//i))
    return data



def tracks(particles='e', filter=None):
    'Returns a list of tracks from the data_folder. Each track is a dictionary with keys: t, x, y, z, ux, uy, uz, q, g, file'
    read_parameters()
    suffix_dict = {'e': '-1', 'p': '1', 'g': '0'}
    track_names = [x for x in os.listdir(data_folder) if x.startswith('track_'+suffix_dict[particles[0]])]
    tracks = [read_track(x) for x in track_names]
    if filter is not None:
        filter_expr = expression_parser.to_polish(filter)
        tracks = [t for t in tracks if expression_parser.evaluate(filter_expr, lambda var_name: t[var_name])]
    return tracks

def read_track(track_name):
    'Reads track from the specified track file. The returned track is a dictionary with keys: t, x, y, z, ux, uy, uz, q, g, file'
    filename = os.path.join(data_folder, track_name)
    raw_data = np.loadtxt(filename)
    raw_track = raw_data.reshape(9, -1, order='F')
    track_size = raw_track[0].size
    track = {'x' : raw_track[1],
             'y' : raw_track[2],
             'z' : raw_track[3],
             'ux' : raw_track[4],
             'uy' : raw_track[5],
             'uz' : raw_track[6],
             'file' : filename,
             'q' : raw_track[0],
             'g' : raw_track[7],
             'chi' : raw_track[8],
             't' : np.linspace(0, dt * (track_size - 1), track_size)}
    track['vx'] = track['ux'] / track['g']
    track['vy'] = track['uy'] / track['g']
    track['vz'] = track['uz'] / track['g']
    return track

def smooth(xs, lr, squeeze = True):
    'Returns smoothed array; *lr* >= len(xs) results no smoothing;\n\
    if *squeze* == True, the length of resulting array is equal to *lr*,\n\
    otherwise the length of the resulting array is equal to len(xs);\n\
    for example, try\n\
    squeeze(range(15), 5)'
    n = len(xs)
    if lr >= n:
        return xs
    else:
        a = np.zeros(n)
        sigma = 1.0 * n / lr
        ns = int(np.ceil(sigma))
        for i in np.arange(ns, n - ns):
            for j in np.arange(-ns, ns + 1):
                a[i] += xs[i+j] * np.cos(0.5 * np.pi * j / sigma)
        tmp = 0
        for j in np.arange(-ns, ns + 1):
            tmp += np.cos(0.5 * np.pi * j / sigma)
        a = a / tmp
        for i in np.arange(ns):
            tmp = 0
            for j in np.arange(-i, i + 1):
                a[i] += xs[i+j] * np.cos(0.5 * np.pi * j / sigma)
                a[-(i + 1)] += xs[-(i+1)+j] * np.cos(0.5 * np.pi * j / sigma)
                tmp += np.cos(0.5 * np.pi * j / sigma)
            a[i] = a[i] / tmp
            a[-(i+1)] = a[-(i+1)] / tmp
        if squeeze == True:
            b = np.zeros(lr)
            b[0] = a[0]
            b[-1] = a[-1]
            for i in np.arange(1, lr - 1):
                x = 1.0 * (n - 1) * i / (lr - 1)
                j = int(np.floor(x))
                x1 = x - j
                x2 = 1 - x1
                b[i] = a[j] * x2 + a[j+1] * x1
            return b
        else:
            return a

def onaxis(filename, sx = 1, sy = 1, sz = 1, av = 'None'):
    'sx, sy, sz are the number of neighboring points using for averaging and smoothing.\n\
    av == \'y\' or av == \'z\' results density integrated along y or z axis, respectively.\n\
    WARNING: if sx != 1, the x-distance between points is alternating.'
    lr = int(nx / sx)
    if filename == 'x':
        a = np.arange(nx) * dx
    else:
        if av == 'y':
            a = np.sum(density(filename), 0) / ny
        elif av == 'z':
            a = np.sum(density(filename, 'xz'), 0) / nz
        else:
            a = np.zeros(nx)
            for i in np.linspace(1 - sy, sy - 1, 2 * sy - 1):
                a += density(filename)[int(ny/2+i),:]
            for i in np.linspace(1 - sz, sz - 1, 2 * sz - 1):
                a += density(filename,'xz')[int(nz/2+i),:]
            a = a / (2 * (sy + sz) - 2)
    return smooth(a, lr)

def reset_globals():
    global dx,dy,dz,dt,nx,ny,nz,output_period,n_ion_populations,icmr,t_end,tr_start,\
    deps,deps_p,deps_ph,deps_i,a0y,a0z,lmbda,ne,xsigma,nfilm,filmwidth,nerflow,\
    Tlflow, mcrlflow, vlflow, Trflow, vrflow, catching, particles_for_output
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
    nfilm = 0
    filmwidth = 0
    nerflow = 0
    Tlflow = 0
    mcrlflow = 0
    vlflow = 0
    Trflow = 0
    vrflow = 0
    catching = False
    dump_photons = False
    particles_for_output = 'e'
    output_mode = 0
