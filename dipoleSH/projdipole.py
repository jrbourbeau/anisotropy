#!/usr/bin/env python

#============================================================================
# File Name     : almTibet.py
# Description   : Optimization code for spherical harmonic fitting of
#                 relint skymaps from IceCube using Tibet data
# Creation Date : 03-04-2016
# Last Modified : Thu 03 Mar 2016 04:16:49 PM CST
# Created By    : James Bourbeau
#============================================================================

import numpy as np
import matplotlib.pyplot as plt
import argparse
import os, re, sys
import healpy as hp

from mapfunctions.mapFunctions import getMap, multifit, smoothMap
from mapfunctions.projFunctions import *

if __name__ == "__main__":

    public_path = '/home/jbourbeau/public_html/figures/2Dfit'
    parser = argparse.ArgumentParser(description='Calculate anisotropy')
    parser.add_argument('--nbins', dest='nbins', type=int, default=24, \
        help='Number of RA bins to be used in projection')
    parser.add_argument('-l','--lmax', dest='lmax', type=int, default=3, \
        help='Maximum l-value to use during 2D fit')
    parser.add_argument('--out', dest='out', default=False, \
        action='store_true', help='Save plots')
    parser.add_argument('-f','--file', dest='file',\
        default='/data/user/jbourbeau/anisotropy/maps/merged/IC_24H_sid_4-4.25GeV.fits',
        help='.fits file to be analyzed')

    args = parser.parse_args()
    kwargs = vars(args).copy()
    # Apply default values if arguments not in input dictionary
    defaults = {'smooth':5, 'stype':'tophat', 'swindow':3,\
                'verbose':False, 'decmin':-90., 'decmax':-25.}
    opts = {k:kwargs[k] for k in kwargs if k not in defaults}
    opts.update({k:defaults[k] for k in defaults if k not in opts})
    relint = getMap(*[opts['file']], mapName='relint', **opts)
    relerr = getMap(*[opts['file']], mapName='relerr', **opts)

    # Get chi2 minimized alm parameters
    lmax = opts['lmax']
    opts['verbose'] = True

    # Get fitmap proj relint data points that Tibet would see
    opts.update({'lmax':1,'ramin':0.,'ramax':360.})
    # opts.update({'lmax':1,'ramin':0.,'ramax':360.,'decmin':-30.,'decmax':90.})
    ra, ri, ra_err, ri_err = getRIRAProj(relint, relerr, **opts)
    print('ri = {}'.format(ri*1e4))
    print('ra = {}'.format(ra))
    popt, perr, chi2 = getHarmonicFitParams(ra, ri, 1)
    amp, amp_err, phase, phase_err = getProjDipole(relint, **opts)
    # projRI = cosFit(ra,*popt[:3])
    # print('amp = {}'.format(amp))
    # print('phase = {}'.format(phase))
    # print('chi2 = {}'.format(chi2))

    # Plot proj relint vs RA
    fig = plt.figure()
    ax = plt.gca()
    tPars = {'fontsize':16}
    ax.set_xlim(0.0,360.)
    # ax.set_ylim(-0.0015,0.0015)
    #ax.set_ylim(-0.0025,0.0025)
    ax.set_xlabel(r'Right Ascension $[^{}]$'.format('{\circ}'), **tPars)
    ax.set_ylabel(r'Relative Intensity',**tPars)
    ax.invert_xaxis()
    plt.plot(ra,cosFit(ra,*popt[:3]),color='blue',label='Dipole Fit')
    plt.errorbar(ra,ri,xerr=ra_err,yerr=ri_err,marker='.',linestyle='None',color='blue', label='IceCube Proj. RI')
    # ax.annotate(r'Dipole = {:1.1e}$\pm${:1.1e}'.format(amp,amp_err), xy=(125, -0.0010))
    # ax.annotate(r'Phase = {:2.1f}$\pm${:1.1f}'.format(phase,phase_err), xy=(125, -0.00115))
    plt.legend(loc=4)
    # plt.show()
    plt.savefig('/home/jbourbeau/public_html/figures/collaboration-meeting-2016/IC_projrelint_24bins_24H_sid_4-4.25GeV.png', dpi=300, bbox_inches='tight')
