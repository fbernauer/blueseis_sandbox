#! /usr/bin/env python

import sys, os
import numpy as np

from matplotlib import pyplot as plt
from matplotlib import rc, font_manager

def plot_data(fog, ramp, result, srate):

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

    fig = plt.figure()
    ax0 = plt.subplot2grid((3, 1), (0, 0))
    ax1 = plt.subplot2grid((3, 1), (1, 0))
    ax2 = plt.subplot2grid((3, 1), (2, 0))
    
    line_fog, = ax0.plot(sec, fog, color=c_fog)
    
    line_ramp, = ax1.plot(sec, ramp, color=c_ramp)
    
    line_result, = ax2.plot(sec, result, color=c_result)

    ax0.set_ylabel('fog [au]')
    ax1.set_ylabel('ramp [au]')
    ax2.set_ylabel('fog [au]')
    ax2.set_xlabel('time [s]')

    ax0.set_title('fog (not corrected)')
    ax1.set_title('ramp')
    ax2.set_title('fog (corrected)')

    plt.subplots_adjust(
    top=0.845,
    bottom=0.083,
    left=0.095,
    right=0.90,
    hspace=0.4,
    wspace=0.2
    )

    plt.show()

def plot_model(fog, ramp, model, mean):
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
    c_fvsr = (0,0,0)
    c_model = (1,0,0)
    c_mean = (0,0,1)

    #define labels
    l_fvsr = 'fof vs. ramp (raw)'
    l_model = 'fog vs. ramp (model)'

    l_mean = 'fog mean value'



    fig = plt.figure()
    ax0 = plt.subplot2grid((1, 1), (0, 0))
    
    line_fvsr, = ax0.plot(ramp, fog, color=c_fvsr, marker='o', linestyle='None', label=l_fvsr, markersize=5)
    line_model, = ax0.plot(range(max(ramp)), model, color=c_model, marker='o', linestyle='None', label=l_model, markersize=2)
    line_mean, = ax0.plot(ramp, np.ones(len(ramp))*mean, color=c_mean, linestyle='-', label=l_mean)
    
    ax0.set_xlabel('ramp [au]')
    ax0.set_ylabel('fog [au]')

    # legend
    lines = (line_fvsr, line_model, line_mean)
    labels = (l_fvsr, l_model, l_mean)
    plt.legend(lines, labels,
            loc='upper left',
            ncol = 1)


    plt.show()


