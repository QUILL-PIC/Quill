#!/usr/bin/python

import matplotlib.cm as cm
import matplotlib.colors as clrs
import numpy as np

__name__ = 'tcmap - provides functions for \
partially-transparent-colormap generation\n'
__doc__ =\
'Provided functions:\n\
\tget(name=\'jet\',N=256,gamma=1)\n\
\tvary(cmap,N=256,gamma=1)\n\
\tred(N=256,gamma=0.4)\n\
\tgreen(N=256,gamma=0.4)\n\
\tblue(N=256,gamma=0.4)\n\
\torange(N=256,gamma=0.4)\n\
\tpurple(N=256,gamma=0.4)\n\
\n\
Function get(name=\'jet\',N=256,gamma=1) returns\n\
:class:`matplotlib.colors.Colormap` instance with alpha channel\n\
changing between fully transparent to fully opaque. *name* corresponds\n\
to one of the predefined names of :class:`matplotlib.colors.Colormap`,\n\
*N* is the number of colormap levels, and *gamma* is the gamma\n\
correction for alpha channel that nevertheless can be set less than\n\
zero: if *gamma*>0 then alpha is proportional to value^{gamma},\n\
otherwise 1-alpha is proportional to value^{-gamma}, hence, in the\n\
latter case low values are plotted as opaque, and high values are\n\
plotted as transparent. Furthermore, if *gamma* is set to 0, opaque\n\
colormap is generated.\n\
\n\
Function vary(cmap,N=256,gamma=1) is the analogue of the function\n\
get(name=\'jet\',N=256,gamma=1) except the first mandatory parameter\n\
cmap that should be set to some :class:`matplotlib.colors.Colormap`\n\
instance.\n\
\n\
The rest functions (red, etc.) return some partially transparent\n\
colormaps that look fine if used togeather. When awoked, these\n\
functions also register corresponding colormaps that since can be\n\
easily used by means of their names (\'tcmap_red\', etc.).  *N* and\n\
*gamma* has the same meaning as for get(...) function.\n\
\n\
Below is the example of tcmap.py usage:\n\
$ ipython --pylab\n\
>> import tcmap\n\
>> a = random.rand(10,10)\n\
>> cmap = tcmap.get(\'hsv\')\n\
>> imshow(a,cmap)\n\
>> colorbar()\n\
>> cmap = tcmap.vary(cmap,gamma=2)\n\
>> imshow(a,cmap)\n\
>> colorbar()\n\
>> hold(True)\n\
>> cmap = tcmap.red()\n\
>> imshow(a,cmap)\n\
>> colorbar()\n\
>> tcmap.blue()\n\
>> imshow(a,\'tcmap_blue\')\n\
>> colorbar()'

def vary(cmap,N=256,gamma=1.0):
    'Varies transparency of a colormap *cmap*.'
    m = cmap
    m._init() # now 2D m._lut array is available
    # m._lut is the lookup table - array that contains [[r,g,b,a],...]
    # data of a colormap as well as _rgba_bad (for masked values),
    # _rgba_over (for high out-of-range values) and _rgba_under (for
    # low out-of-range values) data
    v = np.linspace(0,1,m.N)
    if gamma==0:
        m._lut[:-3,3] = 1
    elif gamma>0:
        m._lut[:-3,3] = v**gamma
    else:
        m._lut[:-3,3] = 1-v**(-gamma)
    return m

def get(name='jet',N=256,gamma=1.0):
    'Partially-transparent-colormap generator.'
    m = cm.get_cmap(name,N)
    m = vary(m,N,gamma)
    return m

def red(N=256,gamma=0.4):
    'Returns and registers \'tcmap_red\' colormap.'
    cdict = { 'red': ((0,1,1),(0.5,1,1),(1,0.7,0.7)), 'green': ((0,0.237,0.237),(0.5,0.237,0.237),(1,0.166,0.166)), 'blue': ((0,0.3,0.3),(0.5,0.3,0.3),(1,0.21,0.21)) }
    m = clrs.LinearSegmentedColormap('tcmap_red',cdict,N=N)
    m = vary(m,gamma=gamma)
    cm.register_cmap(cmap=m)
    return m

def green(N=256,gamma=0.4):
    'Returns and registers \'tcmap_green\' colormap.'
    cdict = { 'red': ((0,0.42,0.42),(0.5,0.42,0.42),(1,0.294,0.294)), 'green': ((0,0.59,0.59),(0.5,0.59,0.59),(1,0.413,0.413)), 'blue': ((0,0,0),(1,0,0)) }
    m = clrs.LinearSegmentedColormap('tcmap_green',cdict,N=N)
    m = vary(m,gamma=gamma)
    cm.register_cmap(cmap=m)
    return m

def blue(N=256,gamma=0.4):
    'Returns and registers \'tcmap_blue\' colormap.'
    cdict = { 'red': ((0,0,0),(1,0,0)), 'green': ((0,0.61,0.61),(0.5,0.61,0.61),(1,0.427,0.427)), 'blue': ((0,1,1),(0.5,1,1),(1,0.7,0.7)) }
    m = clrs.LinearSegmentedColormap('tcmap_blue',cdict,N=N)
    m = vary(m,gamma=gamma)
    cm.register_cmap(cmap=m)
    return m

def orange(N=256,gamma=0.4):
    'Returns and registers \'tcmap_orange\' colormap.'
    cdict = { 'red': ((0,1,1),(0.5,1,1),(1,0.7,0.7)), 'green': ((0,0.28,0.28),(0.5,0.28,0.28),(1,0.196,0.196)), 'blue': ((0,0.07,0.07),(0.5,0.07,0.07),(1,0.049,0.049)) }
    m = clrs.LinearSegmentedColormap('tcmap_orange',cdict,N=N)
    m = vary(m,gamma=gamma)
    cm.register_cmap(cmap=m)
    return m

def purple(N=256,gamma=0.4):
    'Returns and registers \'tcmap_purple\' colormap.'
    cdict = { 'red': ((0,0.454,0.454),(0.5,0.454,0.454),(1,0.318,0.318)), 'green': ((0,0.454,0.454),(0.5,0.454,0.454),(1,0.318,0.318)), 'blue': ((0,0.61,0.61),(0.5,0.61,0.61),(1,0.427,0.427)) }
    m = clrs.LinearSegmentedColormap('tcmap_purple',cdict,N=N)
    m = vary(m,gamma=gamma)
    cm.register_cmap(cmap=m)
    return m

def copper():
    #'Returns \'tcmap_copper\' colormap.'
    m = cm.get_cmap('jet',256)
    m._init()
    c0 = np.array([221,0,0]) #dd0000
    c0 = c0/255.
    c1 = np.array([0,0,0])
    c1 = c1/255.
    alpha = 0.1
    m._lut[:-3,0] = np.linspace(c0[0],c1[0],256)
    m._lut[:-3,1] = np.linspace(c0[1],c1[1],256)
    m._lut[:-3,2] = np.linspace(c0[2],c1[2],256)
    m._lut[:-3,3] = np.linspace(alpha,1,256)
    m._lut[-3,3] = alpha
    m._lut[-3,0] = c0[0]
    m._lut[-3,1] = c0[1]
    m._lut[-3,2] = c0[2]
    m._lut[-2,3] = 1
    m._lut[-2,0] = c1[0]
    m._lut[-2,1] = c1[1]
    m._lut[-2,2] = c1[2]
    m._lut[-1,3] = 1
    m._lut[-1,0] = c1[0]
    m._lut[-1,1] = c1[1]
    m._lut[-1,2] = c1[2]
    return m
