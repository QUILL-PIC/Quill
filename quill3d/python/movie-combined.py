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

qplot.resread.data_folder = '/home/eugn/2014.11.25-positron-trapping/quill3d/results_ne50-a0212/'
qplot.resread.read_parameters()

dt = qplot.resread.output_period
t0 = qplot.resread.tr_start + 1e-2 * np.floor( 1e2 * dt )
t1 = qplot.resread.t_end + 0.5*dt

qplot.resread.v = 0.3

cm = qplot.tcmap.copper()

for i in np.arange( (t1-t0)/qplot.resread.output_period ):
    t = 1e-2*np.floor(1e2*(t0+i*qplot.resread.output_period))
    s = 'frame%05d' % i + '.png'
    print t, s
    #
    plt.clf()
    #
    qplot.density(t,max_w=1e9,max_e_density=3e2,max_g_density=1e9,extent='xi')
    #qplot.density(t,max_w=3e5,max_e_density=4e2,extent='xi')
    #qplot.tracks(['xi','y'],'e',t0=t-qplot.resread.tr_start-qplot.resread.output_period,t1=t-qplot.resread.tr_start)
    qplot.tracks(['xi','y','t'],'e',t0=0,t1=t-qplot.resread.tr_start,cmaps=[cm],r=0.7)
    plt.clim(t-dt,t)
    plt.xlim(30,48)
    plt.savefig(s)


os.system('convert -density 100 -delay 15 -loop 0 frame*.png movie-combined.gif')
os.system('rm frame*.png')
