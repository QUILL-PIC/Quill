#!/usr/bin/python
import numpy as np
import resread

def check(data_folder = '../results/', t=None):
    'Verifies energy conservation in a Quill run'
    resread.data_folder = data_folder
    resread.read_parameters()
    data = resread.t_data('energy', silent = True)
    data_deleted = None
    if resread.catching:
        data_deleted = resread.t_data('energy_deleted', silent=True)
    
    index = -1
    index_del = -1
    if t==None:
        t = data[-1][0]
    else:
        index = next((i for i in range(len(data)) if data[i][0] > t), None)
        if resread.catching:
            index_del = next((i for i in range(len(data_deleted)) if data_deleted[i][0] > t), None)
        
        if index == None or index_del == None:
            raise ValueError("T = {0} is too big; T max = {1}".format(t, data[-1][0]))
    print ("T = {0}, T max = {1}, T end = {2}, data folder = {3}".format(t, data[-1][0], resread.t_end, data_folder))
    print ("dx = {0}, dy = {1}, dz = {2}, dt = {3}, ne = {4}, nfilm = {5}".format(resread.dx, resread.dy, resread.dz, resread.dt, resread.ne, resread.nfilm))
    print ("W in box (t={0}) / W_0 = {1}".format(t, np.sum(data[index][1:]) / np.sum(data[0][1:])))
    if resread.catching:
        print ("W total  (t={0}) / W_0 = {1}".format(t, (np.sum(data[index][1:]) + np.sum(data_deleted[index_del][1:])) / np.sum(data[0][1:])))


