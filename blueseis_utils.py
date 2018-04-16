#! /usr/bin/env python


import sys, os
import numpy as np
import ctypes
import numpy.ctypeslib as npct

from obspy.core import read, Stream, Trace, UTCDateTime

def rotate(u, v, w, M, inverse):

    # compute the inverse of M
    M_inv = np.linalg.inv(M)

    # do the transformation
    if not inverse:
        x, y, z = np.dot(M_inv, [u, v, w])
    else:
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


def deramp(fog, ramp, navg):
# bind in the deramp.so library
    deramplib = ctypes.CDLL('./deramp.so')

# define the argument variables   
    array_2d_double = npct.ndpointer(dtype=np.double, ndim=2, flags='CONTIGUOUS')
    array_2d_int = npct.ndpointer(dtype=np.int64, ndim=2, flags='CONTIGUOUS')

    deramplib.deramp.argtypes = [array_2d_int, array_2d_double, array_2d_double, ctypes.c_int, ctypes.c_int]

    fog_corr = []

    for i in range(len(fog)):
        data = np.ascontiguousarray(np.transpose(np.asarray([ramp[i], fog[i]], dtype=np.int64))) 

# define length of moving average window
        n_avg0 = navg
        n_avg = navg
        if len(data) < n_avg0:
            n_avg = int(len(data) / 2)

        data_orig = np.ascontiguousarray(data.copy(), dtype=np.double)
        data_sort = np.ascontiguousarray(np.empty((len(data), 2), dtype=np.double))

#######################################################
        deramplib.deramp(data, data_orig, data_sort, len(data), n_avg)
#######################################################
           
        fog_corr.append(np.transpose(data_orig)[1])

    return fog_corr

