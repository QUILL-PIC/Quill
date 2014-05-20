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
t1 = qplot.resread.t_end + 0.5*qplot.resread.dt
dt = qplot.resread.dt

for i in np.arange( (t1-t0)/qplot.resread.output_period ):
    t = 1e-2*np.floor(1e2*(t0+i*qplot.resread.output_period))
    s = 'frame%05d' % i + '.png'
    print t, s
    qplot.density(t,max_w=1e9,max_e_density=1e9,max_i_density=5e2,save2=s)
    plt.clf()

os.system('convert -density 100 -delay 15 -loop 0 frame*.png movie.gif')
os.system('rm frame*.png')
