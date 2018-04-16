#! /usr/bin/env python
"""
This is a simple plot script, that reproduces the plots showing the results
"""


import sys, os
import numpy as np

from matplotlib import pyplot as plt
from matplotlib import rc, font_manager
from obspy import read

def plot_data(fog, ramp, result, srate, title):

    sec = np.asarray(range(len(fog)))*(1./srate)
    ###########################################################################
    # bind in Latex to matplotlib
    params = {'text.usetex': True, 
            'text.latex.preamble': [r'\usepackage{cmbright}', r'\usepackage{amsmath}']}
    plt.rcParams['figure.figsize'] = 12, 7
    plt.rcParams.update(params)
    sizeOfFont = 16
    fontProperties = {'family':'sans-serif','sans-serif':['Helvetica'],
    'weight' : 'bold', 'size' : sizeOfFont}
    rc('font',**fontProperties)

    # define colors

    c_fog = (0,0,0)
    c_ramp = (0,0,0)

    c_result = (0.2,0.7,0)

    l_fog = 'raw data'
    l_result = 'deramped data'

    fig = plt.figure()
    ax0 = plt.subplot2grid((2, 1), (0, 0))
    ax1 = plt.subplot2grid((2, 1), (1, 0))
    
    line_ramp, = ax0.plot(sec, ramp, color=c_ramp)
    line_fog, = ax1.plot(sec, fog, color=c_fog, linewidth=2) 
    line_result, = ax1.plot(sec, result, color=c_result, linewidth=1)

    ax1.set_ylabel('fog [au]')
    ax0.set_ylabel('ramp [au]')
    ax1.set_xlabel('time [s]')

    ax0.set_title(title)

    #ax0.set_xlim(20, 30)
    #ax1.set_xlim(20, 30)

    # legend
    lines = (line_fog, line_result)
    labels = (l_fog, l_result)
    plt.legend(lines, labels,
            loc='upper left',
            ncol = 1)

    plt.subplots_adjust(
    top=0.89,
    bottom=0.083,
    left=0.095,
    right=0.90,
    hspace=0.2,
    wspace=0.2
    )

    plt.show()


# channel HJ1:
#st1 = read('../test_data/test_HJ1.mseed')
#stR1 = read('../test_data/test_YR1.mseed')
#stdeR1 = read('../test_data/deramped/XS.BS1..HJ1.D.2018.057.pro')

# channel HJ2:
#st1 = read('../test_data/test_HJ2.mseed')
#stR1 = read('../test_data/test_YR2.mseed')
#stdeR1 = read('../test_data/deramped/XS.BS1..HJ2.D.2018.057.pro')

# channel HJ3:
st1 = read('../test_data/test_HJ3.mseed')
stR1 = read('../test_data/test_YR3.mseed')
stdeR1 = read('../test_data/deramped/XS.BS1..HJ3.D.2018.057.pro')


fog = st1[0].data
ramp = stR1[0].data
result = stdeR1[0].data

srate = st1[0].stats.sampling_rate

title = 'channel: '+st1[0].stats.channel 

plot_data(fog, ramp, result, srate, title)
