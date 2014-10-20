#!/usr/bin/python

import matplotlib as mpl
mpl.use('Agg')
#
import numpy as np
import matplotlib.pyplot as plt
import os
import qplot

__name__ = 'movie'
__doc__ = 'Makes gif movie from quill results'

qplot.resread.data_folder = '../results/'
qplot.resread.read_parameters()

t0 = 0
t1 = 2 #qplot.resread.t_end + 0.5*qplot.resread.dt
dt = 0.015 # qplot.resread.dt

#for i in np.arange( (t1-t0)/qplot.resread.output_period ):
for i in np.arange( (t1-t0)/dt ):
    #t = 1e-2*np.floor(1e2*(t0+i*qplot.resread.output_period))
    t = i*dt
    s = 'frame%05d' % i + '.png'
    print t, s
    #qplot.density(t,max_w=1e9,max_e_density=1e9,max_i_density=1e4,save2=s)
    #qplot.density(t,max_w=1e6,max_e_density=1e9,save2=s)
    #qplot.density(t,max_w=1e9,max_e_density=1e9,max_i_density=20,save2=s)
    #qplot.field(t,'by',fmax=10,save2=s)
    #qplot.tracks(['y','z','t'],'p',t0=t,t1=t+0.6,cmaps=['Reds'],r=0.5,axis=[9,16,10,15],save2=s)
    #qplot.tracks(['x','y','t'],'p',t0=t,t1=t+0.4,cmaps=['Reds'],r=0.5,axis=[13,19,10,15],save2=s)
    #qplot.tracks(['x','y'],'p',t0=0,t1=t,colors=['#dd0000'],axis=[13.6,14.6,13,13.7],save2=s)
    #qplot.tracks(['x','y'],'p',t0=0,t1=t,colors='',axis=[13.6,14.6,13,13.7],save2=s)
    #
    qplot.tracks(['z','y','x'],'p',t0=t,t1=t+0.1,r=1.5,cmaps=['brg_r'],axis=[10,12,10,12])
    cb = plt.colorbar()
    cb.set_label('x')
    plt.clim(13.7,14.7)
    plt.savefig(s)
    #
    plt.clf()

os.system('convert -density 100 -delay 15 -loop 0 frame*.png movie.gif')
os.system('rm frame*.png')
