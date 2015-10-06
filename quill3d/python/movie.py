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

#qplot.resread.data_folder = '../results/'
qplot.resread.data_folder = '/share/dmserebr/Quill/quill3d/results_normal_ions_off_2/results_ne210-a0210/'
#qplot.resread.data_folder = '/share/dmserebr/Quill/quill3d/results_normal_ions_off_2/results_ne105-a0210/'
#qplot.resread.data_folder = '/share/dmserebr/Quill/quill3d/results_normal_ions_off_2/results_ne210-a0105/'
qplot.resread.read_parameters()

t0 = 0
t1 = qplot.resread.t_end + 0.5*qplot.resread.dt
dt = qplot.resread.output_period

for i in np.arange( (t1-t0)/qplot.resread.output_period ):
    t = 1e-2*np.floor(1e2*(t0+i*qplot.resread.output_period))
    s = 'frame%05d' % i + '.png'
    print t, s
    #qplot.density(t,max_w=3e5,max_e_density=800,max_i_density=1e9,max_g_density=1e3,save2=s)
    #qplot.density(t,max_w=3e5,max_e_density=400,max_i_density=1e9,max_g_density=1e3,save2=s)
    #qplot.density(t,max_w=3e5/4,max_e_density=800,max_i_density=1e9,max_g_density=3e2,save2=s)
    #
    plt.subplot(111,axisbg='black')
    qplot.particles(t, ['x','y','chi'], 'e', r = 1, cmap = 'YlGnBu_r', data_folder = qplot.resread.data_folder)
    plt.xlim(11, 13.5)
    plt.ylim(0, 27)
    plt.clim(0, 0.2)
    plt.savefig(s)
    #
    plt.clf()


#os.system('convert -density 100 -delay 15 -loop 0 frame*.png movie.gif')
#os.system('rm frame*.png')
