#!/usr/bin/python

# This script runs quill3d many times for various parameters

import os
import datetime as dt
import rmtlib

import numpy as np # qwe

config = '.laser-piston'

# parameter 0
rmtlib.p0_name = 'ne'
p0_value = rmtlib.geometric_progression(7e22,3e23,2)
rmtlib.p0_unit = ''

# parameter 1
rmtlib.p1_name = 'a0'
p1_value = rmtlib.geometric_progression(400,2500,4)
rmtlib.p1_unit = ''

# see list of pp_operations in rmtlib.py
rmtlib.pp_operation = ['density','spectrum','tracks','i:x-ux','energy'] # post-processing operations

rmtlib.cwd = os.getcwd() # current directory
t = dt.datetime.now()

print rmtlib.p0_name+', '+rmtlib.p0_unit+',',p0_value
print rmtlib.p1_name+', '+rmtlib.p1_unit+',',p1_value
print rmtlib.pp_operation
os.chdir('../')
os.system('./parse.sh '+config+' > conf')

for p0 in p0_value:
    for p1 in p1_value:
	rmtlib.prepare_conf(p0,p1,config)
	print rmtlib.p0_name+' = '+str(p0)+' '+rmtlib.p0_unit
	print rmtlib.p1_name+' = '+str(p1)+' '+rmtlib.p1_unit
	os.system('rm results/*')
	os.system('cat conf | ./quill')
	rmtlib.ppo(p0,p1,config,t.strftime('%y-%m-%d--%H-%M-%S'))
