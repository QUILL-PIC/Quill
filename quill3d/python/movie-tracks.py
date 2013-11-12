#!/usr/bin/python

import matplotlib as mpl
mpl.use('Agg')
#
import numpy as np
import matplotlib.pyplot as plt
import os
import resread
import qplot

__name__ = 'movie'
__doc__ = 'Makes gif movie from quill results'

t0 = 1
t1 = 7

resread.data_folder = '../results/'
resread.read_parameters()

dt = 5*resread.dt

for t in np.arange(t0,t1,dt):
    s = 'frame%05d' % int(t/dt) + '.png'
    print t, s
    #qplot.tracks(space=['x','y','t'],t0=t,t1=t+0.3,particles='ep',cmaps=['Blues','Reds'],axis=[7,11,7,11],r=2,save2=s)
    #qplot.tracks(space=['x','y','t'],t0=t,t1=t+0.3,particles='p',cmaps=['Reds'],axis=[20,30,5,15],r=2,save2=s)
    #qplot.tracks(space=['x','y','t'],t0=t,t1=t+0.3,particles='e',cmaps=['Blues'],axis=[8,13,8,12],save2=s)
    qplot.tracks(['x','y','g'],t0=t,t1=t+0.3,particles='e',cmaps=['gist_earth_r'],axis=[9,12,9,11],save2=s)
    #plt.savefig(s)
    plt.clf()

os.system('convert -density 100 -delay 15 -loop 0 frame*.png movie.gif')
os.system('rm frame*.png')
