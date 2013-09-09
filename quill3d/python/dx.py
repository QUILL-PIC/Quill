#!/usr/bin/python

from numpy import pi
from numpy import sqrt

def dx(ne=0,dt=0):
    'returns dx [lambda] for given *ne* [ncr] and *dt* [lambda]'
    if dt==0:
	return 1/(pi*sqrt(ne))
    else:
	return dt/(1-ne*dt*dt*4*pi*pi/4.04)

ne = 75.7659
dt = 0.014
print 'dt_crit = ',dx(ne)
print 'dx = ',dx(ne,dt)
print '1/dx = ',1/dx(ne,dt)
