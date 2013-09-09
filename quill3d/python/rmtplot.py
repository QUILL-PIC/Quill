#!/usr/bin/python

# This script plots results written by rmt.py

import matplotlib.pyplot as plt
import numpy as np

conf = '.film-and-lp-model'
date = '13-05-05--20-07-27'
pp_operation = ['energy','ph_spectrum']
pname = 'a0'
linestyle = [['-.^y','--sg',':vr','-*b','-.Dm','--k.'],['-co'],['-mo','-mo','-mo']]

xlog = True
ylog = False

n = 0
for name in pp_operation:
    f = open(pname+'-'+name+'--'+conf[1:]+'--'+date)
    data = []
    nlines = 0
    for line in f:
	nlines += 1
	a = line.split('\t')
	for b in a:
	    data.append(float(b))
	if name=='energy':
	    data.append(data[-1]+data[-2]+data[-3]+data[-4]+data[-5]) # it's assumed that only one type of ions is present
    data = np.reshape(data,(nlines,len(data)/nlines))
    if ylog==True:
	for i in np.arange(len(data[:,0])):
	    for j in np.arange(1,len(data[0,:])):
		if data[i,j]>0:
		    data[i,j] = np.log(data[i,j])/np.log(10.)
		else:
		    data[i,j] = -10
    x = data[:,0]
    if pname=='a0':
	for i in np.arange(len(x)):
	    x[i] = x[i]*x[i]*3.32e18 # qwe - WARNING - lambda!
    if xlog==True:
	for i in np.arange(len(x)):
	    x[i] = np.log( x[i] )/np.log(10.)
    for i in np.arange(1,len(data[0,:])):
	plt.plot(x,data[:,i],linestyle[n][i-1])
    n += 1
plt.show()
