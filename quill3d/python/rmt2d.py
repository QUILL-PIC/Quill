#!/usr/bin/python

# This script runs quill3d many times for various parameters

import os
import datetime as dt
import rmtlib

conf = '.film-and-lp-model'

rmtlib.p0name = 'ne'
p0value = rmtlib.geometric_progression(2e22,1e24,10)
rmtlib.p0unit = ''

rmtlib.pname = 'a0'
pvalue = rmtlib.geometric_progression(50,2500,10)
rmtlib.punit = ''

rmtlib.pp_operation = ['density','energy','ph_spectrum','mollweide','ions'] # post-processing operations
rmtlib.ppo_parameter = [6,159,6e-9,0,0] # t for output, energy normalization, etc.

rmtlib.cwd = os.getcwd() # current directory
t = dt.datetime.now()

rmtlib.mode = 'w'

print rmtlib.p0name+', '+rmtlib.p0unit+',',p0value
print rmtlib.pname+', '+rmtlib.punit+',',pvalue
print rmtlib.pp_operation, rmtlib.ppo_parameter
os.chdir('../')
os.system('./parse.sh '+conf+' > conf')

for p0 in p0value:
    rmtlib.prepare_conf(p0,True)
    for p in pvalue:
	rmtlib.prepare_conf(p)
	print rmtlib.p0name+' = '+str(p0)+' '+rmtlib.p0unit
	print rmtlib.pname+' = '+str(p)+' '+rmtlib.punit
	os.system('cat conf | ./quill')
	rmtlib.ppo(p,conf,t.strftime('%y-%m-%d--%H-%M-%S'),p0)
