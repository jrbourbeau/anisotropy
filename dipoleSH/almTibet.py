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
    parser.add_argument('--chi2', dest='chi2', default='RI', \
        help='Chi-squared function to minimize [standard, d1d2, RI, tibetonly]')
    parser.add_argument('--dipole', dest='dipole', default=False, \
        action='store_true', help='Only use dipole components')
    parser.add_argument('--nbins', dest='nbins', type=int, default=18, \
        help='Number of RA bins to be used in projection')
    parser.add_argument('-l','--lmax', dest='lmax', type=int, default=3, \
        help='Maximum l-value to use during 2D fit')
    parser.add_argument('--out', dest='out', default=False, \
        action='store_true', help='Save plots')
    parser.add_argument('-f','--file', dest='file',\
        default='/data/user/jbourbeau/anisotropy/maps/merged/IC_24H_sid_4-4.25GeV.fits',
        help='.fits file to be analyzed')
    parser.add_argument('--step', dest='step', type=float, \
        help='TMinuit initial step size or uncertainty')
    parser.add_argument('--init', dest='init', type=float, \
        help='TMinuit initial parameter value')

    args = parser.parse_args()
    kwargs = vars(args).copy()
    # Apply default values if arguments not in input dictionary
    defaults = {'smooth':5, 'stype':'tophat', 'swindow':3,\
                'verbose':False, 'decmin':-90., 'decmax':-25.}
    opts = {k:kwargs[k] for k in kwargs if k not in defaults}
    opts.update({k:defaults[k] for k in defaults if k not in opts})
    print('opts = {}'.format(opts))
    # sys.exit()
    # data, bg, local = hp.read_m   ap(opts['file'], range(3), verbose=False)
    # data = getMap(*[opts['file']], mapName='data', **opts)
    # bg = getMap(*[opts['file']], mapName='bg', **opts)

    # Get chi2 minimized alm parameters
    lmax = opts['lmax']
    opts['verbose'] = True
    # Useful breakdown of l, m values to be used
    nsph = sum([2*l_i+1 for l_i in range(lmax+1)])
    lvals = [[l_i]*(2*l_i+1) for l_i in range(lmax+1)]
    mvals = [[m for m in range(-l_i, l_i+1)] for l_i in range(lmax+1)]
    lvals = [item for sublist in lvals for item in sublist]
    mvals = [item for sublist in mvals for item in sublist]
    # fitparams, fiterrparams, p = multifit(lmax, data, bg, 1/20., **opts)
    fitparams = ['Y(%i,%i)' % (lvals[i], mvals[i]) for i in range(nsph)]
    fiterrparams = ['dY(%i,%i)' % (lvals[i], mvals[i]) for i in range(nsph)]
    p = np.load('SHCoeff/SHCoeff1e-15_1e-06/{}chi2_lmax_{}.npy'.format(opts['chi2'],lmax))
    # p = np.load('SHcoeff/{}chi2_lmax_{}.npy'.format(opts['chi2'],lmax))
    p = p.item()
    # np.save('SHCoeff{}_{}/{}chi2_lmax_{}.npy'.format(opts['init'], opts['step'], opts['chi2'],lmax), p)
    # sys.exit()

    # Create SH fit to relint map from alm values
    # npix = len(data)
    # nside = hp.npix2nside(npix)
    nside = 64
    npix = hp.nside2npix(nside)
    vx, vy, vz = hp.pix2vec(nside, [i for i in range(npix)])
    normedSH = np.load('/home/jbourbeau/anisotropy/dev/normedSH.npy')
    normedSH = normedSH.item()
    fitmap = np.zeros(npix)
    fiterrmap = np.zeros(npix)
    for i in range(nsph):
        if opts['dipole']:
            if fitparams[i]=='Y(1,1)' or fitparams[i]=='Y(1,-1)':
                print('fitparams[{}] = {}'.format(i,fitparams[i]))
                fitmap += p[fitparams[i]] * normedSH[fitparams[i]]
                fiterrmap += p[fiterrparams[i]] * normedSH[fitparams[i]]
            else:
                continue
        else:
            fitmap += p[fitparams[i]] * normedSH[fitparams[i]]
            fiterrmap += p[fiterrparams[i]] * normedSH[fitparams[i]]


    # Write fitmap to file
    fitsfile = 'fitmap_lmax{}_{}chi2.fits'.format(lmax,opts['chi2'])
    hp.write_map('fitmapsOptimal/'+fitsfile,fitmap,coord='C')
    sys.exit()

    # Get fitmap proj relint data points that Tibet would see
    opts.update({'lmax':1,'ramin':0.,'ramax':360.,'decmin':-30.,'decmax':90.})
    #ra, ri, ra_err, ri_err = getRIRAProj(fitmap, fiterrmap, **opts)
    ra, ri, ra_err = getRIRAProj(fitmap, **opts)
    # Get fitmap dipole amp and phase that Tibet would see
    #amp, amp_err, phase, phase_err = getProjDipole(fitmap,fiterrmap, **opts)
    amp, amp_err, phase, phase_err = getProjDipole(fitmap, **opts)
    popt, perr, chi2 = getHarmonicFitParams(ra, ri, 1)
    projRI = cosFit(ra,*popt[:3])
    print('amp = {}'.format(amp))
    print('phase = {}'.format(phase))
    print('chi2 = {}'.format(chi2))
    # Get Tibet proj relint data points
    tibetdata = np.loadtxt('Tibet_data.dat')
    tibetRA = tibetdata[:,0]
    tibetRAerr = [10.]*18
    tibetRI = tibetdata[:,1]-1.
    tibetRIerr = np.array([2.8e-5]*18)
    tibetpopt, tibetperr, tibetchi2 = getHarmonicFitParams(tibetRA, tibetRI, 1, tibetRIerr)
    tibetprojRI = cosFit(tibetRA,*tibetpopt[:3])

    # Calculate chi-squared for RI data points
    Chi2 = sum((ri - tibetRI)**2/(tibetRIerr)**2)
    redChi2 = Chi2/18.
    print('\nChi2 = {}'.format(Chi2))
    print('redChi2 = {}\n'.format(redChi2))

    # Plot proj relint vs RA
    fig = plt.figure()
    ax = plt.gca()
    #ax = fig.axes[0]
    tPars = {'fontsize':16}
    ax.set_xlim(0.0,360.)
    ax.set_ylim(-0.0015,0.0015)
    #ax.set_ylim(-0.0025,0.0025)
    ax.set_xlabel(r'Right Ascension', **tPars)
    ax.set_ylabel(r'Relative Intensity',**tPars)
    ax.invert_xaxis()
    plt.errorbar(ra,ri,xerr=ra_err,marker='.',linestyle='None',color='blue',
        label='Fitmap Proj. RI')
    plt.plot(ra,cosFit(ra,*popt[:3]),color='blue',label='Dipole Fit')
    ax.annotate(r'{} $\chi^2$'.format(opts['chi2']), xy=(125, -0.00085))
    ax.annotate(r'Dipole = {:1.1e}$\pm${:1.1e}'.format(amp,amp_err), xy=(125, -0.0010))
    ax.annotate(r'Phase = {:2.1f}$\pm${:1.1f}'.format(phase,phase_err), xy=(125, -0.00115))
    plt.errorbar(tibetRA,tibetRI,xerr=tibetRAerr,yerr=tibetRIerr,marker='.',
        linestyle='None',color='green',label='Tibet Proj. RI')
    plt.plot(tibetRA,cosFit(tibetRA,*tibetpopt[:3]),color='green',
        label='Tibet Dipole Fit')
    plt.legend()
    plt.show()
    #ax.annotate(r'Standard $\chi^2$', xy=(150, -0.008))
    #ax.annotate(r'Dipole = {:1.1e}$\pm${:1.1e}'.format(amp,amp_err), xy=(150, -0.010))
    #ax.annotate(r'Phase = {:2.1f}$\pm${:1.1f}'.format(phase,phase_err), xy=(150, -0.012))
    # plt.savefig('/home/jbourbeau/public_html/figures/2Dfit/tibet_relint_lmax{}_{}chi2.png'.format(lmax,opts['chi2']), dpi=300, bbox_inches='tight')
