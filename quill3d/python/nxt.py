#!/usr/bin/python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
import resread

df = '../results/'
yav = 1
inp = 'bicubic' # interpolation
cmap = 'bwr'
iaa = 2 # /int for antialiasing
f2draw = 'by'

def onaxis(data_folder,filename,t,yav=1):
    resread.data_folder = data_folder
    resread.read_parameters()
    resread.t = t
    a = np.zeros(resread.nx)
    if filename=='x':
	for i in np.arange(resread.nx):
	    a[i] = i*resread.dx
    else:
	for j in np.arange(yav):
	    a += resread.density(filename)[resread.ny/2+j,:]/yav
    return a

def smooth(arr,sigma=1):
    n0 = len(arr)
    n = np.floor(n0/sigma)
    delta = float(n0-1)/(n-1)
    a = np.zeros(n)
    for i in np.arange(1,n-1):
	tmp = 0
	for j in np.arange( np.floor((i-1)*delta), np.floor((i+1)*delta) ):
	    b = np.cos( -0.5*np.pi*(i-j/delta)**2 )**2
	    tmp += b
	    a[i] += arr[j]*b
	a[i] = a[i]/tmp
    return a

resread.data_folder = df
resread.read_parameters()

sigma =  resread.output_period/resread.dx/iaa
print 'sigma = ', sigma

tmp = np.array([])
#for t in np.arange(0, 39, resread.output_period):
for t in np.arange(0, resread.t_end+0.5*resread.dt, resread.output_period):
    print str('%g'%t)
    if f2draw=='rho':
	n = -onaxis(df,f2draw,str('%g'%t),yav)
    else:
	n = onaxis(df,f2draw,str('%g'%t),yav)
    n = smooth(n,sigma)
    tmp = np.append(tmp,n)
tmp = tmp.reshape(len(tmp)/len(n),len(n))
tmp = tmp.transpose() # n in (i: x-t, j: t)

nxt = tmp
#norm = resread.ne/(np.pi/(2.818e-13*resread.lmbda*resread.lmbda))
#nxt = nxt/norm

f = plt.figure()
ax = f.add_subplot(111)

extent = (0,resread.t_end,0,resread.nx*resread.dx)
plt.imshow(nxt,cmap=cmap,interpolation=inp,origin='lower',extent=extent)
plt.clim(-10,10)
plt.colorbar()

plt.xlabel(r'$ct/\lambda$')
plt.ylabel(r'$x/\lambda$')

plt.xlim(0,int(resread.t_end))

plt.show()
