#! /usr/bin/env python
"""

"""

import sys, os
import numpy as np
import ctypes
import numpy.ctypeslib as npct
import argparse

from obspy.core import read, Stream, Trace, UTCDateTime
from blueseis_utils import rotate, deramp, deramp_median, deramp_mode, plot_deramp

#"""
#The script can handle day-files (or shorter) in miniseed format.
#
#The output files will be named after miniseed convention: NN.SSSSS.LL.CCC.D.YYYY.DDD.pro (note the suffix .pro!)
#
#The input rotation rate and ramp files have to contain the same time spans.
#
#"""

######################################################################################################################
# This script takes exactly one rotation rate file and the corresponding ramp file as input. 
# It deramps the traces using the ramp traces.
#
# The script can handle day-files (or shorter) in miniseed format.
#
# The output files will be named after miniseed convention: NN.SSSSS.LL.CCC.D.YYYY.DDD.pro (note the suffix .pro!)
#
# The input rotation rate and ramp files have to contain the same time spans.
#######################################################################################################################


def check_time_uvw(st_rot_u, st_ramp_u):
    for i in range(len(st_rot_u)):
        if st_rot_u[i].stats.starttime != st_rot_v[i].stats.starttime\
        or st_rot_u[i].stats.starttime != st_rot_w[i].stats.starttime\
        or st_rot_w[i].stats.starttime != st_rot_v[i].stats.starttime\
        or st_rot_u[i].stats.endtime != st_rot_v[i].stats.endtime\
        or st_rot_u[i].stats.endtime != st_rot_w[i].stats.endtime\
        or st_rot_w[i].stats.endtime != st_rot_v[i].stats.endtime:

            print("Warning: Traces do not have the same start or end time. Traces will be trimmed!")
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
        print("Warning: Rotation rate and ramp data for trace "+str(i)+" of channel "+tr_rot.stats.channel+" do not have same time stamps!")
        tr_rot.trim(t1, t2)
        tr_ramp.trim(t1, t2)
    if tr_rot.stats.channel[-1] != tr_ramp.stats.channel[-1]:
        print("ERROR: Rotation rate and ramp data for trace "+str(i)+" do not have same channel id")
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

def create_outfile(st_in, out_folder, method, suff, overwrite):
    net = st_in[0].stats.network
    sta = st_in[0].stats.station
    loc = st_in[0].stats.location
    cha = st_in[0].stats.channel
    dq = 'D'
    yyyy = str(st_in[0].stats.starttime.year)
    ddd = str(st_in[0].stats.starttime.julday).zfill(3)
    if not suff == '':
        fname = net+'.'+sta+'.'+loc+'.'+cha+'.'+dq+'.'+yyyy+'.'+ddd+'.pro_'+method+'.'+suff
    else:
        fname = net+'.'+sta+'.'+loc+'.'+cha+'.'+dq+'.'+yyyy+'.'+ddd+'.pro_'+method
        
    outfname = out_folder+fname
    # check if out file exists
    if os.path.isfile(outfname) and overwrite == 0:
        print('ERROR: file '+outfname+' already exists!')
        sys.exit(1)
    if os.path.isfile(outfname) and overwrite == 1:
        print('Warning: file '+outfname+' will be overwritten!')
    #    outfnames.append(outfname)
    return open(outfname, 'wb')

def get_samples(tr1, tr2, start, stop):
    if len(tr1) <= stop or len(tr2) <= stop:
        GoOn = False
    else:
        GoOn = True
    data1 = tr1.data[start:stop]
    data2 = tr2.data[start:stop]
    return data1, data2, GoOn


def process_data(infnames_rot, infnames_ramp, out_folder, n_samp, overwrite, navg, method, suff, plot_flagg):
# check if output directory exists. If not, create it.
    if not out_folder.endswith('/'):
        out_folder = out_folder+'/'
    if not os.path.exists(out_folder):
        os.makedirs(out_folder)    


# check list of input files
    #if len(infnames_rot) != 3 and len(infnames_ramp) != 3:
    #    print("ERROR: You must process at least three input channels!")
    #    sys.exit(1)

# read input files
    st_rot = Stream()
    st_ramp = Stream()
    for infname_rot in infnames_rot:
        st_rot += read(infname_rot)
    for infname_ramp in infnames_ramp:
        st_ramp += read(infname_ramp)

    st_rot_u = st_rot.select(channel='HJ1')
    out_file_u = create_outfile(st_rot_u, out_folder, method, suff, overwrite)
    st_ramp_u = st_ramp.select(channel='YR1')


# loop over traces
    
    for i in range(len(st_rot_u)):
        GoOn = True

        tr_u, start_u, stop_u, n_max_u = check_traces(st_rot_u[i], st_ramp_u[i], i, n_samp)
        print(st_rot_u[i])
        print(st_ramp_u[i])

        print('')
        print('processing trace '+str(i+1)+' of '+str(len(st_rot_u)))
        while GoOn:
# get data
            sys.stdout.write('U: start: %10i | stop: %10i | total %15i\r\n' %(start_u, stop_u, n_max_u))                
            sys.stdout.flush()
            fu, ru, GoOn = get_samples(st_rot_u[i], st_ramp_u[i], start_u, stop_u)

            start_u = stop_u
            stop_u += n_samp
            if stop_u > n_max_u:
                stop_u = n_max_u

#######################################################

            fog = [fu]
            ramp = [ru]
            
            sys.stdout.write('deramping ...\r\n')
            if method == 'mean':
                f_corr, x, ramp_sort, fog_sort, mean = deramp(fog, ramp, navg, plot_flagg)
            if method == 'median':
                f_corr, x, ramp_sort, fog_sort, mean = deramp_median(fog, ramp, navg, plot_flagg)
            if method == 'mode':
                f_corr, x, ramp_sort, fog_sort, mean = deramp_mode(fog, ramp, navg, plot_flagg)

            sys.stdout.write('rotating back again ...\r\n')
#######################################################
            if plot_flagg:
                plot_deramp(fog, ramp, f_corr, x, ramp_sort, fog_sort, mean)
    
            tr_u.data = np.append(tr_u.data, f_corr)
# write processed data
        out_file_u.seek(0, 2)
        tr_u.write(out_file_u, format='mseed', reclen=512, encoding='FLOAT64')
        out_file_u.flush()

    print('')
    out_file_u.flush()
    out_file_u.close()




def main():
    parser = argparse.ArgumentParser(description=__doc__, prog='process_bs.py')
    
    parser.add_argument('-F', metavar='', type=str, nargs=1, required=True,
                        default='', dest='fogfiles', help='three files containing raw rotation rates, in minseed')
    parser.add_argument('-R', metavar='', type=str, nargs=1, required=True,
                        default=5, dest='rampfiles', help='three files containing ramp data, in minseed')
    parser.add_argument('-O', metavar='', type=str,
                        default='./', dest='outfolder', help='output directory')
    parser.add_argument('-l', metavar='', type=int,
                        default=75000, dest='nsamp', help='[int] length of window to be processed, in samples')
    parser.add_argument('-o', default=False,
                        action='store_true', dest='ow', help='if "-o" option is set, existing output file will be over written')
    parser.add_argument('-m', metavar='', type=int,
                        default=100, dest='navg', help='[int] length of moving average window, in samples')
    parser.add_argument('-a', metavar='', type=str, required=True,
                        default='', dest='method', help='specify method: "mean", "median" or "mode"')
    parser.add_argument('-s', metavar='', type=str, required=False,
                        default='', dest='suff', help='specify any suffix to the out put file name.')
    parser.add_argument('-p', metavar='', type=int, required=False,
                        default=0, dest='plot', help='[0 or 1]: results will be plotted or not.')
    
    args = parser.parse_args()

    rotfnames = args.fogfiles
    rampfnames = args.rampfiles
    out_folder = args.outfolder
    n_samp = args.nsamp
    if args.ow:
        ow = 1
    if not args.ow:
        ow = 0
    n_avg = args.navg
    method = args.method
    suff = args.suff
    plot_flagg = args.plot

    process_data(rotfnames, rampfnames, out_folder, n_samp, ow, n_avg, method, suff, plot_flagg)

if __name__ == "__main__":
    main()
