import numpy as np
import numpy.ctypeslib as nc

ncc = nc.ctypes

lib = nc.load_library('chameleon-core', '.')

# fills up conf (internal chameleon global variable) with values from log file,
# prints variable name and value if flag is True
def configure(log, flag = False):
    f = lib.configure
    f.argtypes = [ncc.c_char_p, ncc.c_bool]
    f(log, flag)

def get(name):
    f = lib.get
    f.argtypes = [ncc.c_char_p]
    f.restype = ncc.c_double
    return np.double(f(name))

def read(filename, plane):
    f = lib.read
    f.argtypes = [ncc.c_char_p, ncc.c_char_p]
    nx = int(get('nx'))
    ny = int(get('ny'))
    nz = int(get('nz'))
    if plane == "xy":
        f.restype = nc.ndpointer(ncc.c_double, shape = (nx * ny,))
        v = f(filename, plane)
        v = np.reshape(v, (nx, ny))
    elif plane == "xz":
        f.restype = nc.ndpointer(ncc.c_double, shape = (nx * nz,))
        v = f(filename, plane)
        v = np.reshape(v, (nx, nz))
    elif plane == "yz":
        f.restype = nc.ndpointer(ncc.c_double, shape = (ny * nz,))
        v = f(filename, plane)
        v = np.reshape(v, (ny, nz))
    return np.transpose(v)

