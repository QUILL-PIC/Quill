# this script plots results obtained by parallel.cpp of run_config.cpp

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import resread

def parallel(t, resread_function, param0 = 'ne', param1 = 'a0', search_in =
    '../'):
  'returns tuple tau such that tau[0] is the unsorted (!) parameter array\
  [[param0, param1], ...], and tau[1] is the array of results returned by\
  resread_function (energy, spectrum, etc.)'
  nfiles = 0
  for filename in os.listdir(search_in):
    if filename.find('results_') == 0:
      nfiles += 1
  tau = []
  tau1 = []
  tau.append(np.empty((nfiles,2)))
  i = 0
  for filename in os.listdir(search_in):
    if filename.find('results_') == 0:
      print filename
      ip0 = filename.find( param0 )
      ip1 = filename.find( param1 )
      im = filename.find( "-" )
      # value of parameters
      tau[0][i,0] = float( filename[ip0+len( param0 ):im] )
      tau[0][i,1] = float( filename[ip1+len( param1 ):] )
      i += 1
      #
      if resread_function == 'spectrum':
        resread.data_folder = search_in + filename + '/'
        resread.read_parameters()
        # [[epsilon[0], dN/depsilon[0]],...]
        tau1.append(resread.t_data('spectrum' + '%g' % t, resread.deps))
  tau.append(np.array(tau1))
  return tau

tau = parallel(40, 'spectrum')

# slow sort; us element, us[i], defines index of (sorted) i-th param0 value in
# initially unsorted array tau[0][:,0]
arrts = tau[0][:,0] # array to sort
us = np.arange(len(arrts))
print arrts
for i in np.arange(len(arrts)):
  for j in np.arange(len(arrts)-1):
    if arrts[j+1] < arrts[j]: # resulting arr values increase with index
      tmp = arrts[j+1]
      arrts[j+1] = arrts[j]
      arrts[j] = tmp
      tmp = us[j+1]
      us[j+1] = us[j]
      us[j] = tmp
print arrts
print us

n1 = len(tau[1])
n2 = len(tau[1][0,:])
print 'img shape:', n1, n2
sp = np.empty((n1, n2))
for i in np.arange(n1):
  for j in np.arange(n2):
    sp[i,j] = tau[1][us[i],j,1]
# exclusion of cold electrons (nthr starting points)
nthr = 20
for i in np.arange(n1):
  for j in np.arange(0, nthr):
    sp[i,j] = 0
# normalization
#intsp = np.zeros(n2)
#for i in np.arange(n1):
#  for j in np.arange(n2):
#    intsp[i] += sp[i,j]
#for i in np.arange(n1):
#  if intsp[i] != 0:
#    sp[i] /= intsp[i]
#for i in np.arange(n1):
#  if sp[i].max() != 0:
#    sp[i] /= sp[i].max()

lw = 0.7 # linewidth
plt.figure(figsize=(3.5,3),dpi=300)
font = {'family' : 'Liberation Serif',
        'weight' : 'normal',
        'size'   : 10}
matplotlib.rc('font', **font)
matplotlib.rc('lines', linewidth=lw)
matplotlib.rc('axes', linewidth=lw)

epsmax = tau[1][0,-1,0]
nemin = tau[0][:,0].min()
nemax = tau[0][:,0].max()
#print epsmax, nemin, nemax
extent = (0, epsmax, 0, n1 - 1)
aspect = float(epsmax / (n1 - 1))
plt.imshow(sp,cmap='jet',interpolation='none',origin='lower',extent=extent,aspect=aspect)
cb = plt.colorbar()
cb.set_label(r'$dN_e/d\epsilon$, MeV$^{-1}$')
plt.clim(0,1e9)
plt.xlabel(r'$(\gamma - 1) mc^2$, MeV')
#plt.yticks([2.6, 4, 5.9, 8.9, 13, 20,
#  30],['0.021','0.024','0.028','0.037','0.052','0.095','0.22'])
plt.yticks([2, 10, 15, 18, 21, 23, 26.5, 30], ['0.02', '0.04', '0.06', '0.08',
  '0.1', '0.12', '...', '0.22'])
plt.ylabel(r'$n_e/n_{cr}$')
plt.savefig('plot_parallel.pdf',bbox_inches='tight')
