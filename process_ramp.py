#! /usr/bin/env python

"""


USAGE: ./process_ramp.py <rotfile> <rampfile> <out_dir> <npts> <overwrite> <n_avg>


        rotfile :  input file containing rotation rate data (in miniseed format)
        rampfile:  input file containing ramp data (in miniseed format)
        out_dir :  path to output directory
        npts    :  number of samples used for the correction
        overwrite: 0 or 1. if 1 existing file will be overwritten!
        n_avg   :  window length for moving average

"""

import sys, os
import numpy as np
import matplotlib.pyplot as plt
import ctypes
import numpy.ctypeslib as npct

from scipy.signal import resample
from obspy.core import read, Trace, UTCDateTime


def process_ramp(infname_rot, infname_ramp, out_folder, n_samp, overwrite, navg):
# bind in the deramp.so library
    deramplib = ctypes.CDLL('./deramp.so')

# define the argument variables   
    array_2d_double = npct.ndpointer(dtype=np.double, ndim=2, flags='CONTIGUOUS')
    array_2d_int = npct.ndpointer(dtype=np.int64, ndim=2, flags='CONTIGUOUS')

    deramplib.deramp.argtypes = [array_2d_int, array_2d_double, array_2d_double, ctypes.c_int, ctypes.c_int]

# check if output directory exists. If not, create it.
    if not out_folder.endswith('/'):
        out_folder = out_folder+'/'
    if not os.path.exists(out_folder):
        os.makedirs(out_folder)    

# check if out file exists
    outfname = out_folder+infname_rot.split('/')[-1]+".pro"
    if os.path.isfile(outfname) and overwrite == 0:
        print 'ERROR: file '+outfname+' already exists!'
        sys.exit(1)
    if os.path.isfile(outfname) and overwrite == 1:
        print 'Warning: file '+outfname+' will be overwritten!'

    out_file = open(outfname, 'wb')

# read input files
    st_rot = read(infname_rot)
    st_ramp = read(infname_ramp)

# loop over traces
    for i in range(len(st_rot)):
        GoOn = True
        tr = Trace(header=st_rot[i].stats)

# check start times and channels
        if st_rot[i].stats.starttime != st_ramp[i].stats.starttime:
            print "ERROR: Rotation rate and ramp data for trace "+str(i)+" do not have same start time!"
            sys.exit(1)
        if st_rot[i].stats.channel[-1] != st_ramp[i].stats.channel[-1]:
            print "ERROR: Rotation rate and ramp data for trace "+str(i)+" do not have same channel id"
            sys.exit(1)

# handle traces with different length
        if len(st_rot[i].data) != len(st_ramp[i].data):
            n_max = min(len(st_rot[i].data), len(st_ramp[i].data))
        else:
            n_max = len(st_rot[i].data)

        start = 0
        stop = n_samp
        if stop > n_max:
            stop = n_max
        k = 0
        print ''
        print 'processing trace '+str(i+1)+' of '+str(len(st_rot))
        while GoOn:
# get data
            sys.stdout.write('start: %10i | stop: %10i | total %15i\r' %(start, stop, n_max))                
            sys.stdout.flush()
            data, GoOn = get_samples(st_rot[i], st_ramp[i], start, stop)

# define length of moving average window
            n_avg0 = navg
            n_avg = navg
            if len(data) < n_avg0:
                n_avg = int(len(data) / 2)

            start = stop
            stop += n_samp
            if stop > n_max:
                stop = n_max
            data_orig = np.ascontiguousarray(data.copy(), dtype=np.double)
            data_sort = np.ascontiguousarray(np.empty((len(data), 2), dtype=np.double))

#######################################################
            deramplib.deramp(data, data_orig, data_sort, len(data), n_avg)
#######################################################
           
            data_corr = np.transpose(data_orig)[1]

            tr.data = np.append(tr.data, data_corr)
# write processed data
        out_file.seek(0, 2)
        tr.write(out_file, format='mseed', reclen=512, encoding='FLOAT64')
        out_file.flush()
    
    print ''
    out_file.flush()
    out_file.close()







def get_samples(tr1, tr2, start, stop):
    if len(tr1) <= stop or len(tr2) <= stop:
        GoOn = False
    else:
        GoOn = True
    data1 = tr1.data[start:stop]
    data2 = tr2.data[start:stop]
    return np.ascontiguousarray(np.transpose(np.asarray([data2, data1], dtype=np.int64))), GoOn


def main():
    if len(sys.argv) < 6:
        print __doc__
        sys.exit(1)

    rotfname = sys.argv[1]
    rampfname = sys.argv[2]
    out_folder = sys.argv[3]
    n_samp = int(sys.argv[4])
    ow = int(sys.argv[5])
    n_avg = int(sys.argv[6])

    process_ramp(rotfname, rampfname, out_folder, n_samp, ow, n_avg)

if __name__ == "__main__":
    main()
