#!/usr/bin/python

# This script plots results written by rmt2d.py

import matplotlib.pyplot as plt
import numpy as np

data_folder = ''

p0name = 'ne'
pname = 'a0'
conf = '.film-and-lp-model'
date = '13-06-11--17-15-36'

zlog = False #True
linewidth = 1.5

lmbda = 0.91e-4
l = 2*np.pi # in lmbda/2\pi
tau = 3.3*2*np.pi # FWHM in lmbda/2\pi
mp = 1836
imcr = 2.0

ncr = 9.11e-28*(3e10*2*np.pi/lmbda)**2/(4*np.pi*4.8e-10**2)

extent = (0,1,0,1)

# a0 -> intensity normalized on 1 W/cm^2
a = 2/np.log(10.)
b = np.log(1e-7*np.pi*9.11e-28**2*3e10**5/4.8e-10**2/lmbda**2)/np.log(10.)
# ne -> density normalized on 1 cm^{-3}
c = 1/np.log(10.)
d = 0

def lgI(a0):
    return a*np.log(a0) + b

def lgne(ne):
    return c*np.log(ne) + d

def get(pp_operation,column2plot,get_extent=False):
    global extent
    f = open(data_folder+p0name+'-'+pname+'-'+pp_operation+'--'+conf[1:]+'--'+date)
    i = 0
    data = []
    for line in f:
	a1 = line.split('\t')
	i+=1
	for b1 in a1:
	    data.append(float(b1))
    f.close()
    data = np.reshape(data,(i,len(data)/i))
    i = 0
    while True:
	i += 1
	if data[i,1]==data[0,1]:
	    break
    j = int(len(data[:,0])/i)
    data2plot = np.reshape(data[:j*i,column2plot],(j,i))
    if get_extent==True:
	extent = (lgI(data[0,1]),lgI(data[i-1,1]),lgne(data[0,0]),lgne(data[j*i-1,0]))
	aa = extent[0]
	bb = extent[1]
	cc = extent[2]
	dd = extent[3]
	extent = (aa-0.5*(bb-aa)/(i-1),bb+0.5*(bb-aa)/(i-1),cc-0.5*(dd-cc)/(j-1),dd+0.5*(dd-cc)/(j-1))
    if zlog==True:
	for i in np.arange(len(data2plot[:,0])):
	    for j in np.arange(len(data2plot[0,:])):
		data2plot[i,j] = np.log(data2plot[i,j])/np.log(10.)
    return data2plot

lng1 = get('mollweide',6,True)
lng2 = get('mollweide',8)

lng = lng2 - lng1

f = plt.figure()
ax = f.add_subplot(111)
plt.title('tmp')
plt.xlabel('lg I[W/cm^2]')
plt.ylabel('lg n_e[cm^{-3}]')
plt.xlim(extent[0],extent[1])
plt.ylim(extent[2],extent[3])
aspect = (extent[1]-extent[0])/(extent[3]-extent[2])
im = ax.imshow(lng,'Spectral_r',aspect=aspect,interpolation='nearest',origin='lower',extent=extent) # ,vmin=-3,vmax=1)
f.colorbar(im)

n = 10
a0arr = np.linspace(25,2500,n)
x = np.empty(n)
y = np.empty(n)
for i in np.arange(n):
    x[i] = lgI(a0arr[i])
    y[i] = lgne(2*a0arr[i]*ncr/l)
ax.plot(x,y,'--k',linewidth=linewidth)
ax.text(22.75,22.5,'S=2/l')

a0arr = np.linspace(150,2500,n)
x = np.empty(n)
y = np.empty(n)
for i in np.arange(n):
    x[i] = lgI(a0arr[i])
    y[i] = lgne(8*np.pi*ncr*a0arr[i]**2/mp/imcr)
ax.plot(x,y,'--k',linewidth=linewidth)
ax.text(23,23.7,'4')

for i in np.arange(n):
    y[i] = lgne(ncr/(2*np.pi*l*2)*mp*imcr)
ax.plot(x,y,'--k',linewidth=linewidth)
ax.text(24.5,22.7,'6')

plt.show()
