#! /usr/bin/env python


import sys, os
import numpy as np
import matplotlib
matplotlib.use("Qt4Agg")
import matplotlib.pyplot as plt
from matplotlib.pylab import *
from matplotlib import rc, font_manager
import ctypes
import numpy.ctypeslib as npct

from obspy.core import read, Stream, Trace, UTCDateTime


def rotate(u, v, w, M, inverse):
# rotates from UVW system to XYZ system by multiplying with the rotation matrix M
    # compute the inverse of M
    M_inv = np.linalg.inv(M)

    # do the transformation
    if inverse:
        x, y, z = np.dot(M_inv, [u, v, w])
    if not inverse:
        x, y, z = np.dot(M, [u, v, w])


    # Replace all negative zeros. These might confuse some further
    # processing programs.
    y = np.array(y).ravel()
    y[y == -0.0] = 0
    z = np.array(z).ravel()
    z[z == -0.0] = 0
    x = np.array(x).ravel()
    x[x == -0.0] = 0

 
    return x, y, z

def plot_deramp(fog, ramp, fog_corr, x, fog_sort, ramp_sort, _mean):
    params = {'text.usetex': True, 
          'text.latex.preamble': [r'\usepackage{cmbright}', r'\usepackage{amsmath}']}
    plt.rcParams.update(params)
    sizeOfFont = 15
    fontProperties = {'family':'sans-serif','sans-serif':['Helvetica'],
        'weight' : 'bold', 'size' : sizeOfFont}
    rc('font',**fontProperties)

    x_axis = []
    mean = []
    for i in range(3):
        x_axis.append(np.arange(len(fog[i])))
        mean.append(np.ones(len(ramp_sort[i]))*_mean[i])

    fig = plt.figure(figsize=(16, 16))
    ax0 = []
    ax1 = []
    for i in range(3):
        ax0.append(plt.subplot2grid((2, 3), (0, i)))
        ax1.append(plt.subplot2grid((2, 3), (1, i)))

    for i in range(3):
        line_raw, = ax0[i].plot(x_axis[i], fog[i], color='k')
        line_corr, = ax0[i].plot(x_axis[i], fog_corr[i], color='g')

        line_fr, = ax1[i].plot(ramp[i], fog[i], 'o', markersize=1, color='k')
        line_sort, = ax1[i].plot(ramp_sort[i], fog_sort[i], 'o', markersize=1.5, color='r')
        line_fr_corr, = ax1[i].plot(ramp[i], fog_corr[i], 'o', markersize=1, color='y')
        line_mod, = ax1[i].plot(ramp_sort[i], x[i], 'o', markersize=1, color='g')
        line_m, = ax1[i].plot(ramp_sort[i], mean[i], color='b')
    
    label_raw = 'raw'
    label_corr = 'corrected'
    label_fr = 'raw'
    label_sort = 'sorted'
    label_mod = 'model'
    label_fr_corr = 'corrected'
    label_m = 'average'
    
    lines0 = (line_raw, line_corr)
    labels0 = (label_raw, label_corr)
    lines1 = (line_fr, line_sort, line_mod, line_fr_corr, line_m)
    labels1 = (label_fr, label_sort, label_mod, label_fr_corr, label_m)

    ax0[0].set_title('HJ1')
    ax0[1].set_title('HJ2')
    ax0[2].set_title('HJ3')

    ax0[0].set_ylabel('rate')
    ax1[0].set_ylabel('fog')

    ax0[0].set_xlabel('samples')
    ax0[1].set_xlabel('samples')
    ax0[2].set_xlabel('samples')

    ax1[0].set_xlabel('ramp')
    ax1[1].set_xlabel('ramp')
    ax1[2].set_xlabel('ramp')

    bta0 = (0.5, 0.5, 0.5, 0.5)
    bta1 = (0., 0.75, 1.0, 0.25)

    ax0[2].legend(lines0, labels0, loc='best', bbox_to_anchor=bta0)
    ax1[2].legend(lines1, labels1, loc='best', bbox_to_anchor=bta1, ncol=2)

    plt.show()



def deramp(fog, ramp, navg, plot_flagg):
# bind in the deramp.so library
    deramplib = ctypes.CDLL('/home/fbernauer/deRamp_test/blueseis_sandbox-master/deramp.so')

# define the argument variables   
    array_2d_double = npct.ndpointer(dtype=np.double, ndim=2, flags='CONTIGUOUS')
    array_1d_double = npct.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS')
    array_2d_int = npct.ndpointer(dtype=np.int64, ndim=2, flags='CONTIGUOUS')

    deramplib.deramp.argtypes = [array_2d_double, array_2d_double, array_2d_double, array_1d_double, ctypes.c_int, ctypes.c_int, array_1d_double]

    fog_corr = []
    x = []
    ramp_sort = []
    fog_sort = []
    mean = []

    for i in range(len(fog)):
        data = np.ascontiguousarray(np.transpose(np.asarray([ramp[i], fog[i]], dtype=np.double))) 

# define length of moving average window
        n_avg0 = navg
        n_avg = navg
        if len(data) < n_avg0:
            n_avg = int(len(data) / 2)

        data_orig = np.ascontiguousarray(data.copy(), dtype=np.double)
        data_sort = np.ascontiguousarray(np.zeros((len(data), 2), dtype=np.double))
        
        _mean = np.ascontiguousarray(0.0, dtype=np.double)
        _x = np.ascontiguousarray(np.zeros(len(data)), dtype=np.double)

#######################################################
        n = deramplib.deramp(data, data_orig, data_sort, _x, len(data), n_avg, _mean)
#######################################################

        if plot_flagg:
            x.append(_x[:n])
            ramp_sort.append(np.transpose(data_sort)[1][:n])
            fog_sort.append(np.transpose(data_sort)[0][:n])
            mean.append(_mean)

            fog_corr.append(np.transpose(data_orig)[1])
            
        if not plot_flagg:
            fog_corr.append(np.transpose(data_orig)[1])

    return fog_corr, x, ramp_sort, fog_sort, mean
    

def deramp_median(fog, ramp, navg, plot_flagg):
# bind in the deramp.so library
    deramplib = ctypes.CDLL('/home/fbernauer/deRamp_test/blueseis_sandbox-master/deramp_median.so')

# define the argument variables   
    array_2d_double = npct.ndpointer(dtype=np.double, ndim=2, flags='CONTIGUOUS')
    array_1d_double = npct.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS')
    array_2d_int = npct.ndpointer(dtype=np.int64, ndim=2, flags='CONTIGUOUS')

    deramplib.deramp_median.argtypes = [array_2d_double, array_2d_double, array_2d_double, array_1d_double, ctypes.c_int, ctypes.c_int, array_1d_double]

    fog_corr = []
    x = []
    ramp_sort = []
    fog_sort = []
    mean = []

    for i in range(len(fog)):
        data = np.ascontiguousarray(np.transpose(np.asarray([ramp[i], fog[i]], dtype=np.double))) 

# define length of moving average window
        n_avg0 = navg
        n_avg = navg
        if len(data) < n_avg0:
            n_avg = int(len(data) / 2)

        data_orig = np.ascontiguousarray(data.copy(), dtype=np.double)
        data_sort = np.ascontiguousarray(np.empty((len(data), 2), dtype=np.double))

        _mean = np.ascontiguousarray(0.0, dtype=np.double)
        _x = np.ascontiguousarray(np.zeros(len(data)), dtype=np.double)

#######################################################
        n = deramplib.deramp_median(data, data_orig, data_sort, _x, len(data), n_avg, _mean)
#######################################################
           
        if plot_flagg:
            x.append(_x[:n])
            ramp_sort.append(np.transpose(data_sort)[1][:n])
            fog_sort.append(np.transpose(data_sort)[0][:n])
            mean.append(_mean)

            fog_corr.append(np.transpose(data_orig)[1])
            
        if not plot_flagg:
            fog_corr.append(np.transpose(data_orig)[1])

    return fog_corr, x, ramp_sort, fog_sort, mean

def deramp_mode(fog, ramp, navg, plot_flagg):
# bind in the deramp.so library
    deramplib = ctypes.CDLL('/home/fbernauer/deRamp_test/blueseis_sandbox-master/deramp_mode.so')

# define the argument variables   
    array_2d_double = npct.ndpointer(dtype=np.double, ndim=2, flags='CONTIGUOUS')
    array_1d_double = npct.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS')
    array_2d_int = npct.ndpointer(dtype=np.int64, ndim=2, flags='CONTIGUOUS')

    deramplib.deramp_mode.argtypes = [array_2d_double, array_2d_double, array_2d_double, array_1d_double, ctypes.c_int, ctypes.c_int, array_1d_double]

    fog_corr = []
    x = []
    ramp_sort = []
    fog_sort = []
    mean = []

    for i in range(len(fog)):
        data = np.ascontiguousarray(np.transpose(np.asarray([ramp[i], fog[i]], dtype=np.double))) 

# define length of moving average window
        n_avg0 = navg
        n_avg = navg
        if len(data) < n_avg0:
            n_avg = int(len(data) / 2)

        data_orig = np.ascontiguousarray(data.copy(), dtype=np.double)
        data_sort = np.ascontiguousarray(np.empty((len(data), 2), dtype=np.double))

        _mean = np.ascontiguousarray(0.0, dtype=np.double)
        _x = np.ascontiguousarray(np.zeros(len(data)), dtype=np.double)

#######################################################
        n = deramplib.deramp_mode(data, data_orig, data_sort, _x, len(data), n_avg, _mean)
#######################################################
           
        if plot_flagg:
            x.append(_x[:n])
            ramp_sort.append(np.transpose(data_sort)[1][:n])
            fog_sort.append(np.transpose(data_sort)[0][:n])
            mean.append(_mean)

            fog_corr.append(np.transpose(data_orig)[1])
            
        if not plot_flagg:
            fog_corr.append(np.transpose(data_orig)[1])

    return fog_corr, x, ramp_sort, fog_sort, mean

