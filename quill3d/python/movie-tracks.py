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
t1 = 2
dt = 0.015

for i in np.arange( (t1-t0)/dt ):
    t = i*dt
    s = 'frame%05d' % i + '.png'
    print t, s
    qplot.tracks(['z','y','x'],'p',t0=t,t1=t+0.1,r=1.5,cmaps=['brg_r'],axis=[10,12,10,12])
    cb = plt.colorbar()
    cb.set_label('x')
    plt.clim(13.7,14.7)
    plt.savefig(s)
    #
    plt.clf()

os.system('convert -density 100 -delay 15 -loop 0 frame*.png movie.gif')
os.system('rm frame*.png')
