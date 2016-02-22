#!/usr/bin/env python

#==============================================================================
# File Name     : projFunctions.py
# Description   : Store functions commonly used in RA projection and  
#                 harmonic fitting for calcuating dipole amplitude/phase
# Creation Date : 02-19-2016
# Last Modified : Mon 22 Feb 2016 12:17:01 PM CST
# Created By    : James Bourbeau
#==============================================================================

import numpy as np
import healpy as hp
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from mapFunctions import getMap

class re_order_errorbarHandler(mpl.legend_handler.HandlerErrorbar):
    def create_artists(self, *args, **kwargs):
        a = mpl.legend_handler.HandlerErrorbar.create_artists(self, 
                *args, **kwargs)
        a = a[-1:] + a[:-1]
        return a


def getRIRAProj(file, **opts):
    
    # Get relint and relerr maps
    relint = getMap(*file, mapName='relint', **opts)
    relerr = getMap(*file, mapName='relerr', **opts)
    # Setup right-ascension bins
    deg2rad = np.pi / 180
    ramin = opts['ramin'] * deg2rad
    ramax = opts['ramax'] * deg2rad
    rabins, rabinwidth = np.linspace(ramin, ramax, opts['nbins']+1,retstep=True)
    # Calculate phi for each pixel
    npix  = len(relint)
    nside = hp.npix2nside(npix)
    theta, phi = hp.pix2ang(nside, range(npix))
    # Bin in right ascension
    phiBins = np.digitize(phi, rabins) - 1
    # UNSEEN cut
    pass_UNSEEN_cut = (relint != hp.UNSEEN)

    ri, ri_err = np.zeros((2,opts['nbins']))
    for i in range(opts['nbins']):
        pass_phiBin_cut = (phiBins == i)
        pass_all_cuts = pass_UNSEEN_cut * pass_phiBin_cut
        ri[i] = np.mean(relint[pass_all_cuts])
        ri_err[i] = np.sqrt(np.sum(relerr[pass_all_cuts]**2))/pass_all_cuts.sum()
        #ri_err = np.sqrt(data * (bg + 1./20 * data)/bg**3)
    dx = rabinwidth/2.
    ra = np.linspace(ramin+dx, ramax-dx, opts['nbins']) / deg2rad
    ra_err = dx * np.ones(opts['nbins']) / deg2rad

    return (ra, ri, ra_err, ri_err)


def lineFit(x, y, sigmay):
    p = np.polyfit(x, y, 0, w=ri_err)
    fit = np.poly1d(p)
    fitVals = fit(np.asarray(y))
    chi2 = (1. / (len(y)-1)) * sum((y - fitVals)**2 / sigmay**2)
    return fitVals, p, chi2

def cosFit(x, *p):
    deg2rad = 2*np.pi / 360
    return sum([p[2*i+1] * np.cos(deg2rad*(i+1)*(x-p[2*i+2])) \
            for i in range(len(p)/2)]) + p[0]

def getHarmonicFitParams(x, y, l, sigmay):
    
    # Guess at best fit parameters
    amplitude = (3./np.sqrt(2)) * np.std(y)
    phase     = 0
    parm_init = [amplitude, phase]*l
    # Reduce amplitude as we go to higher l values
    for i in range(0,len(parm_init)/2):
        parm_init[2*i] *= 2.**(-i)
    parm_init = [0] + parm_init
    
    deg2rad = 2*np.pi / 360
    fitfunc = lambda x, *p: sum([p[2*i+1] * np.cos(deg2rad*(i+1)*(x-p[2*i+2])) \
            for i in range(len(p)/2)]) + p[0]
    # Do best fit
    popt, pcov = curve_fit(fitfunc, x, y, p0=parm_init, sigma=sigmay)
    fitVals = fitfunc(x, *popt)
    ndof  = len(popt)
    chi2 = (1. / (len(y)-ndof)) * sum((y - fitVals)**2 / sigmay**2)
    perr = np.sqrt(np.diag(pcov))
    
    return popt, perr, chi2


def getProjDipole(file, **opts):

    # Project the relint map in RA
    ra, ri, ra_err, ri_err = getRIRAProj(file,**opts)
    # Fit projected RI to cos harmonic functions
    popt, perr, chi2 = getHarmonicFitParams(ra, ri, opts['lmax'], ri_err)
    
    # Extract dipole amp/phase and errors
    a = np.reshape(popt[1:], (-1,2))
    amp, phase = a[0]
    e = np.reshape(perr[1:], (-1,2))
    amp_err, phase_err = e[0]

    return amp, amp_err, phase, phase_err

def calcAvgBkg(ax, datamap, basename, **opts):

    nside = hp.get_nside(datamap)
    npix = hp.nside2npix(nside)

    counts = zeros(4*nside, dtype=float)
    norm = zeros(4*nside, dtype=float)

    ringId = hp.pix2ring(nside, np.array(range(npix)))
    for i, ring in enumerate(ringId):
        counts[ring] += datamap[i]
        norm[ring] += 1.

    bgmap = zeros(npix, dtype=float)
    mapmax = max(datamap)
    for i, ring in enumerate(ringId):
        bgmap[i] = counts[ring] / norm[ring]

    ra, ri, sigmax, sigmay = returnRI(bgmap, datamap, **opts)
    ax.errorbar(ra, ri, xerr=0*sigmax, yerr=sigmay, marker='o', fmt='.',\
            capsize=4, label=basename+" (avg)", linewidth=2, markersize=8, mew=0)

