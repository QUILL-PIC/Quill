import time
import matplotlib.pyplot as plt
#import resread
import chameleon

t0 = time.clock()
## resread together with text mode
## bubble: 5.1 M, 0.45 s
#resread.read_parameters()
#resread.t = '0'
#ne = -resread.density('rho', 'xy')
#
# chameleon
# text mode; 5.1 M, 0.21 s, binary: 4.7 M, 0.02 s
chameleon.configure('../results3/log', True)
ne = -chameleon.read2d('../results3/rho27', 'yz')
t1 = time.clock()

print('reading time:', t1 - t0, 's')

plt.imshow(ne, 'jet', interpolation='none', origin='lower')
plt.colorbar()
#plt.savefig('example-py.png')
plt.show()
