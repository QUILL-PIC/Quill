#!/usr/bin/python

# This script runs quill3d many times for various parameters

import os
import datetime as dt
import rmtlib

conf = '.film-and-lp-model'

rmtlib.pname = 'a0'
pvalue = [80,160,320,640,1280,2560]
rmtlib.punit = ''

rmtlib.pp_operation = ['density','energy','ph_spectrum','mollweide'] # post-processing operations
rmtlib.ppo_parameter = [6,175,2e-8,0] # t for output, energy normalization, etc.

rmtlib.cwd = os.getcwd() # current directory
t = dt.datetime.now()

rmtlib.mode = 'w'

print rmtlib.pname+', '+rmtlib.punit+',',pvalue
print rmtlib.pp_operation, rmtlib.ppo_parameter
os.chdir('../')
os.system('./parse.sh '+conf+' > conf')

for p in pvalue:
    rmtlib.prepare_conf(p)
    print rmtlib.pname+' = '+str(p)+' '+rmtlib.punit
    os.system('cat conf | ./quill')
    rmtlib.ppo(p,conf,t.strftime('%y-%m-%d--%H-%M-%S'))
