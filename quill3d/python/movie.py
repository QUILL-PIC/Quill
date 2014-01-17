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

resread.data_folder = '../results/'
resread.read_parameters()

t0 = resread.output_period
t1 = resread.t_end

for i in np.arange(1+np.floor(0.5+(t1-t0)/resread.output_period)):
    t = 1e-2*np.floor(1e2*(t0+i*resread.output_period))
    s = 'frame%05d' % int(t/float(resread.output_period)) + '.png'
    print t, s
    #qplot.density(t,save2=s)
    #qplot.particles(t,['x','g'],'e',axis=[9,12,0,100],save2=s)
    qplot.mollweide(t,save2=s)
    plt.savefig(s)
    plt.clf()

os.system('convert -density 100 -delay 15 -loop 0 frame*.png movie.gif')
os.system('rm frame*.png')
