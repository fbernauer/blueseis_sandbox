#! /usr/bin/env python
"""
This script takes exactly three rotation rate files and the corresponding three ramp files as input. 
It rotates the three rotation rate traces into the original blueseis system (by using the required rotation matrix),
deramps the traces useing the ramp traces,
and rotates the three rotation rate traces back.

The script can handle day-files (or shorter) in miniseed format.

The output files will be named after miniseed convention: NN.SSSSS.LL.CCC.D.YYYY.DDD.pro (note the suffix .pro!)

The input rotation rate and ramp files have to cantain the same time spans.

"""

import sys, os
import numpy as np
import ctypes
import numpy.ctypeslib as npct
import argparse

from obspy.core import read, Stream, Trace, UTCDateTime
from blueseis_utils import rotate, deramp


def check_time_uvw(st_rot_u, st_rot_v, st_rot_w, st_ramp_u, st_ramp_v, st_ramp_w):
    for i in range(len(st_rot_u)):
        if st_rot_u[i].stats.starttime != st_rot_v[i].stats.starttime\
        or st_rot_u[i].stats.starttime != st_rot_w[i].stats.starttime\
        or st_rot_w[i].stats.starttime != st_rot_v[i].stats.starttime\
        or st_rot_u[i].stats.endtime != st_rot_v[i].stats.endtime\
        or st_rot_u[i].stats.endtime != st_rot_w[i].stats.endtime\
        or st_rot_w[i].stats.endtime != st_rot_v[i].stats.endtime:

            print "Warning: Traces do not have the same start or end time. Traces will be trimmed!"
            t0 = max(st_rot_u[i].stats.starttime, st_rot_v[i].stats.starttime, st_rot_w[i].stats.starttime)
            t1 = min(st_rot_u[i].stats.endtime, st_rot_v[i].stats.endtime, st_rot_w[i].stats.endtime)
            st_rot_u[i].trim(t0, t1)
            st_rot_v[i].trim(t0, t1)
            st_rot_w[i].trim(t0, t1)

def check_traces(tr_rot, tr_ramp, i, n_samp):
    tr = Trace(header=tr_rot.stats)

# check start times and channels
    if tr_rot.stats.starttime != tr_ramp.stats.starttime or tr_rot.stats.endtime != tr_ramp.stats.endtime:
        t1 = max(tr_rot.stats.starttime, tr_ramp.stats.starttime)
        t2 = min(tr_rot.stats.endtime, tr_ramp.stats.endtime)
        print "Warning: Rotation rate and ramp data for trace "+str(i)+" of channel "+tr_rot.stats.channel+" do not have same time stamps!"
        tr_rot.trim(t1, t2)
        tr_ramp.trim(t1, t2)
    if tr_rot.stats.channel[-1] != tr_ramp.stats.channel[-1]:
        print "ERROR: Rotation rate and ramp data for trace "+str(i)+" do not have same channel id"
        sys.exit(1)

# handle traces with different length
    if len(tr_rot.data) != len(tr_ramp.data):
        n_max = min(len(tr_rot.data), len(tr_ramp.data))
    else:
        n_max = len(tr_rot.data)
    start = 0
    stop = n_samp
    if stop > n_max:
        stop = n_max
    return tr, start, stop, n_max

def create_outfile(st_in, out_folder):
    net = st_in[0].stats.network
    sta = st_in[0].stats.station
    loc = st_in[0].stats.location
    cha = st_in[0].stats.channel
    dq = 'D'
    yyyy = str(st_in[0].stats.starttime.year)
    ddd = str(st_in[0].stats.starttime.julday).zfill(3)
    fname = net+'.'+sta+'.'+loc+'.'+cha+'.'+dq+'.'+yyyy+'.'+ddd+'.pro'
    outfname = out_folder+fname
    return open(outfname, 'wb')

def get_samples(tr1, tr2, start, stop):
    if len(tr1) <= stop or len(tr2) <= stop:
        GoOn = False
    else:
        GoOn = True
    data1 = tr1.data[start:stop]
    data2 = tr2.data[start:stop]
    return data1, data2, GoOn


def process_data(infnames_rot, infnames_ramp, out_folder, n_samp, overwrite, navg, matrix_file):
# check if output directory exists. If not, create it.
    if not out_folder.endswith('/'):
        out_folder = out_folder+'/'
    if not os.path.exists(out_folder):
        os.makedirs(out_folder)    

# check if out file exists
    outfnames = []
    for infname_rot in infnames_rot:
        outfname = out_folder+infname_rot.split('/')[-1]+".pro"
        if os.path.isfile(outfname) and overwrite == 0:
            print 'ERROR: file '+outfname+' already exists!'
            sys.exit(1)
        if os.path.isfile(outfname) and overwrite == 1:
            print 'Warning: file '+outfname+' will be overwritten!'
        outfnames.append(outfname)

# check list of input files
    if len(infnames_rot) != 3 and len(infnames_ramp) != 3:
        print "ERROR: You must process at least three input channels!"
        sys.exit(1)

# read input files
    st_rot = Stream()
    st_ramp = Stream()
    for infname_rot in infnames_rot:
        st_rot += read(infname_rot)
    for infname_ramp in infnames_ramp:
        st_ramp += read(infname_ramp)

    st_rot_u = st_rot.select(channel='HJ1')
    out_file_u = create_outfile(st_rot_u, out_folder)
    st_rot_v = st_rot.select(channel='HJ2')
    out_file_v = create_outfile(st_rot_v, out_folder)
    st_rot_w = st_rot.select(channel='HJ3')
    out_file_w = create_outfile(st_rot_w, out_folder)
    st_ramp_u = st_ramp.select(channel='YR1')
    st_ramp_v = st_ramp.select(channel='YR2')
    st_ramp_w = st_ramp.select(channel='YR3')

# read rotation matrix file
    #M = np.transpose(np.loadtxt(matrix_file, skiprows=1))
    M = np.loadtxt(matrix_file, skiprows=1)

# check start and end times of u, v, w:
    check_time_uvw(st_rot_u, st_rot_v, st_rot_w, st_ramp_u, st_ramp_v, st_ramp_w)


# loop over traces
    
    for i in range(len(st_rot_u)):
        GoOn = True

        tr_u, start_u, stop_u, n_max_u = check_traces(st_rot_u[i], st_ramp_u[i], i, n_samp)
        tr_v, start_v, stop_v, n_max_v = check_traces(st_rot_v[i], st_ramp_v[i], i, n_samp)
        tr_w, start_w, stop_w, n_max_w = check_traces(st_rot_w[i], st_ramp_w[i], i, n_samp)
        print st_rot_u[i]
        print st_ramp_u[i]
        print st_rot_v[i]
        print st_ramp_v[i]
        print st_rot_w[i]
        print st_ramp_w[i]

        print ''
        print 'processing trace '+str(i+1)+' of '+str(len(st_rot_u))
        while GoOn:
# get data
            sys.stdout.write('U: start: %10i | stop: %10i | total %15i\r\n' %(start_u, stop_u, n_max_u))                
            sys.stdout.write('V: start: %10i | stop: %10i | total %15i\r\n' %(start_v, stop_v, n_max_v))                
            sys.stdout.write('W: start: %10i | stop: %10i | total %15i\r\n' %(start_w, stop_u, n_max_w))                
            sys.stdout.flush()
            fu, ru, GoOn = get_samples(st_rot_u[i], st_ramp_u[i], start_u, stop_u)
            fv, rv, GoOn = get_samples(st_rot_v[i], st_ramp_v[i], start_v, stop_v)
            fw, rw, GoOn = get_samples(st_rot_w[i], st_ramp_w[i], start_w, stop_w)

            start_u = stop_u
            stop_u += n_samp
            if stop_u > n_max_u:
                stop_u = n_max_u

            start_v = stop_v
            stop_v += n_samp
            if stop_v > n_max_v:
                stop_v = n_max_v

            start_w = stop_w
            stop_w += n_samp
            if stop_w > n_max_w:
                stop_w = n_max_w

#######################################################
            sys.stdout.write('rotating ...\r\n')
            fx, fy, fz = rotate(fu, fv, fw, M, True)

            fog = [fx, fy, fz]
            ramp = [ru, rv, rw]
            
            sys.stdout.write('deramping ...\r\n')
            f_corr = deramp(fog, ramp, navg)

            sys.stdout.write('rotating back again ...\r\n')
            fu_corr, fv_corr, fw_corr = rotate(f_corr[0], f_corr[1], f_corr[2], M, False)
#######################################################

            tr_u.data = np.append(tr_u.data, fu_corr)
            tr_v.data = np.append(tr_v.data, fv_corr)
            tr_w.data = np.append(tr_w.data, fw_corr)
# write processed data
        out_file_u.seek(0, 2)
        tr_u.write(out_file_u, format='mseed', reclen=512, encoding='FLOAT64')
        out_file_u.flush()

        out_file_v.seek(0, 2)
        tr_v.write(out_file_v, format='mseed', reclen=512, encoding='FLOAT64')
        out_file_v.flush()

        out_file_w.seek(0, 2)
        tr_w.write(out_file_w, format='mseed', reclen=512, encoding='FLOAT64')
        out_file_w.flush()
    
    print ''
    out_file_u.flush()
    out_file_u.close()
    out_file_v.flush()
    out_file_v.flush()
    out_file_w.close()
    out_file_w.close()







def main():
    parser = argparse.ArgumentParser(description=__doc__, prog='process_blueseis.py')
    
    parser.add_argument('-F', metavar='', type=str, nargs=3, required=True,
                        default='', dest='fogfiles', help='three files containing raw rotation rates, in minseed')
    parser.add_argument('-R', metavar='', type=str, nargs=3, required=True,
                        default=5, dest='rampfiles', help='three files containing ramp data, in minseed')
    parser.add_argument('-O', metavar='', type=str,
                        default='./', dest='outfolder', help='output directory')
    parser.add_argument('-l', metavar='', type=int,
                        default=75000, dest='nsamp', help='[int] length of window to be processed, in samples')
    parser.add_argument('-o', metavar='', type=int,
                        default=0, dest='ow', help='[0 or 1] 0: existing output file will not be over written, 1: existing output fill will be over written')
    parser.add_argument('-m', metavar='', type=int,
                        default=100, dest='navg', help='[int] length of moving average window, in samples')
    parser.add_argument('-M', metavar='', type=str, required=True,
                        default='', dest='matrixfile', help='file containing BlueSeis rotation matrix')
    
    args = parser.parse_args()

    rotfnames = args.fogfiles
    rampfnames = args.rampfiles
    out_folder = args.outfolder
    n_samp = args.nsamp
    ow = args.ow
    n_avg = args.navg
    matrix_file = args.matrixfile

    process_data(rotfnames, rampfnames, out_folder, n_samp, ow, n_avg, matrix_file)

if __name__ == "__main__":
    main()
