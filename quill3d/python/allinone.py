#!/usr/bin/python

# This script plots results written by rmt2d.py

import matplotlib.pyplot as plt
import numpy as np

p0name = 'ne'
pname = 'a0'
conf = '.film-and-lp-model'
date = '13-06-11--17-15-36'

lmbda = 0.91e-4
zlog = False

Slines = []

interp = 'nearest'
cmap = 'GnBu'

extent = (0,1,0,1)

def get(pp_operation,column2plot,get_extent=False):
    global extent
    f = open(p0name+'-'+pname+'-'+pp_operation+'--'+conf[1:]+'--'+date)
    i = 0
    data = []
    for line in f:
	a = line.split('\t')
	i+=1
	for b in a:
	    data.append(float(b))
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
	# a0 -> intensity normalized on 1 W/cm^2
	a = 2/np.log(10.)
	b = np.log(1e-7*np.pi*9.11e-28**2*3e10**5/4.8e-10**2/lmbda**2)/np.log(10.)
	# ne -> density normalized on 1 cm^{-3}
	c = 1/np.log(10.)
	d = 0
	extent = (a*np.log(data[0,1])+b,a*np.log(data[i-1,1])+b,c*np.log(data[0,0])+d,c*np.log(data[j*i-1,0])+d)
	a = extent[0]
	b = extent[1]
	c = extent[2]
	d = extent[3]
	extent = (a-0.5*(b-a)/(i-1),b+0.5*(b-a)/(i-1),c-0.5*(d-c)/(j-1),d+0.5*(d-c)/(j-1))
    if zlog==True:
	for i in np.arange(len(data2plot[:,0])):
	    for j in np.arange(len(data2plot[0,:])):
		data2plot[i,j] = np.log(data2plot[i,j])/np.log(10.)
    return data2plot

def Sdraw(axis):
    def f(lgne,S):
	a = 2/np.log(10.)
	b = np.log(1e-7*np.pi*9.11e-28**2*3e10**5/4.8e-10**2/lmbda**2)/np.log(10.)
	return a*(lgne*np.log(10.) - np.log(np.pi/(lmbda**2*2.818e-13)) - np.log(S))+b
    for S in Slines:
	axis.plot([f(extent[2],S),f(extent[3],S)],[extent[2],extent[3]],'--k',linewidth=2)

l_energy = get('energy',2,True)
e_energy = get('energy',3)
p_energy = get('energy',4)
ph_energy = get('energy',5)
i_energy = get('energy',6)
sum_energy = l_energy + e_energy + p_energy + ph_energy + i_energy
print 'Error:', sum_energy.min(), sum_energy.max()

alpha = get('ph_spectrum',2)
a0 = get('energy',1)
chi = np.multiply(np.power(a0,2),2*np.pi*3.862e-11/lmbda)
# it is supposed that Fperp=a0 and gamma=sqrt(alpha)*a0
chi = np.multiply(np.power(alpha,0.5),chi)
print 'Chi:', chi.min(), chi.max()

f = plt.figure()
ax = f.add_subplot(231)
plt.title('Photon energy')
plt.xlabel('lg I[W/cm^2]')
plt.ylabel('lg n_e[cm^{-3}]')
plt.xlim(extent[0],extent[1])
plt.ylim(extent[2],extent[3])
aspect = (extent[1]-extent[0])/(extent[3]-extent[2])
ax.imshow(ph_energy,cmap,aspect=aspect,interpolation=interp,origin='lower',extent=extent)
Sdraw(ax)

ax = f.add_subplot(232)
plt.title('Chi')
plt.xlabel('lg I[W/cm^2]')
#plt.ylabel('lg n_e[cm^{-3}]')
plt.xlim(extent[0],extent[1])
plt.ylim(extent[2],extent[3])
aspect = (extent[1]-extent[0])/(extent[3]-extent[2])
ax.imshow(l_energy,cmap,aspect=aspect,interpolation=interp,origin='lower',extent=extent)
Sdraw(ax)

ax = f.add_subplot(233)
plt.title('Summ energy')
plt.xlabel('lg I[W/cm^2]')
#plt.ylabel('lg n_e[cm^{-3}]')
plt.xlim(extent[0],extent[1])
plt.ylim(extent[2],extent[3])
aspect = (extent[1]-extent[0])/(extent[3]-extent[2])
ax.imshow(sum_energy,cmap,aspect=aspect,interpolation=interp,origin='lower',extent=extent)
Sdraw(ax)

ax = f.add_subplot(234)
plt.title('Electron energy')
plt.xlabel('lg I[W/cm^2]')
plt.ylabel('lg n_e[cm^{-3}]')
plt.xlim(extent[0],extent[1])
plt.ylim(extent[2],extent[3])
aspect = (extent[1]-extent[0])/(extent[3]-extent[2])
ax.imshow(e_energy,cmap,aspect=aspect,interpolation=interp,origin='lower',extent=extent)
Sdraw(ax)

ax = f.add_subplot(235)
plt.title('Positron energy')
plt.xlabel('lg I[W/cm^2]')
#plt.ylabel('lg n_e[cm^{-3}]')
plt.xlim(extent[0],extent[1])
plt.ylim(extent[2],extent[3])
aspect = (extent[1]-extent[0])/(extent[3]-extent[2])
ax.imshow(p_energy,cmap,aspect=aspect,interpolation=interp,origin='lower',extent=extent)
Sdraw(ax)

ax = f.add_subplot(236)
plt.title('Ion energy')
plt.xlabel('lg I[W/cm^2]')
#plt.ylabel('lg n_e[cm^{-3}]')
plt.xlim(extent[0],extent[1])
plt.ylim(extent[2],extent[3])
aspect = (extent[1]-extent[0])/(extent[3]-extent[2])
ax.imshow(i_energy,cmap,aspect=aspect,interpolation=interp,origin='lower',extent=extent)
Sdraw(ax)

plt.show()
