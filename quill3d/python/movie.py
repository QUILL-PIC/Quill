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

t0 = 0
t1 = 12

resread.data_folder = '../results/'
resread.read_parameters()

for t in np.arange(t0,t1+0.5*resread.output_period,resread.output_period):
    s = 'frame%05d' % int(t/float(resread.output_period)) + '.png'
    print t, s
    #qplot.density(t,'xy',axis=[8,13,8,12],save2=s)
    qplot.particles(t,['x','g'],'e',axis=[9,12,0,100],save2=s)
    plt.savefig(s)
    plt.clf()

os.system('convert -density 100 -delay 15 -loop 0 frame*.png movie.gif')
os.system('rm frame*.png')
