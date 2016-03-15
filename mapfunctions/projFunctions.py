#!/usr/bin/env python

#==============================================================================
# File Name     : projFunctions.py
# Description   : Store functions commonly used in RA projection and
#                 harmonic fitting for calcuating dipole amplitude/phase
# Creation Date : 02-19-2016
# Last Modified : Thu 03 Mar 2016 04:16:49 PM CST
# Created By    : James Bourbeau
#==============================================================================

import numpy as np
import healpy as hp
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib import rc

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

class re_order_errorbarHandler(mpl.legend_handler.HandlerErrorbar):
    def create_artists(self, *args, **kwargs):
        a = mpl.legend_handler.HandlerErrorbar.create_artists(self,
                *args, **kwargs)
        a = a[-1:] + a[:-1]
        return a


def getRIRAProj(relint, relerr=None, **opts):

    # Cut to desired dec range (equiv to healpy zenith range)
    deg2rad = np.pi / 180
    npix  = len(relint)
    nside = hp.npix2nside(npix)
    theta, phi = hp.pix2ang(nside, range(npix))
    thetamax = (90 - opts['decmin']) * deg2rad
    thetamin = (90 - opts['decmax']) * deg2rad
    pass_dec_cut = (theta <= thetamax) * (theta >= thetamin)
    # Setup right-ascension bins
    ramin = opts['ramin'] * deg2rad
    ramax = opts['ramax'] * deg2rad
    rabins, rabinwidth = np.linspace(ramin, ramax, opts['nbins']+1,retstep=True)
    # Calculate phi for each pixel
    theta, phi = hp.pix2ang(nside, range(npix))
    # Bin in right ascension
    phiBins = np.digitize(phi, rabins) - 1
    # UNSEEN cut
    pass_UNSEEN_cut = (relint != hp.UNSEEN)
    pass_INF_cut = np.array([True]*len(pass_UNSEEN_cut))
    if relerr is not None:
        pass_INF_cut = (relerr != np.inf)

    ri, ri_err = np.zeros((2,opts['nbins']))
    for i in range(opts['nbins']):
        pass_phiBin_cut = (phiBins == i)
        pass_all_cuts = pass_UNSEEN_cut*pass_phiBin_cut*pass_INF_cut*pass_dec_cut
        ri[i] = np.mean(relint[pass_all_cuts])
        if relerr is not None:
            ri_err[i] = np.sqrt(np.sum(relerr[pass_all_cuts]**2))/pass_all_cuts.sum()
        #ri_err = np.sqrt(data * (bg + 1./20 * data)/bg**3)
    dx = rabinwidth/2.
    ra = np.linspace(ramin+dx, ramax-dx, opts['nbins']) / deg2rad
    ra_err = dx * np.ones(opts['nbins']) / deg2rad

    if relerr is None:
        return (ra, ri, ra_err)
    else:
        return (ra, ri, ra_err, ri_err)


def cosFit(x, *p):
    deg2rad = 2*np.pi / 360
    return sum([p[2*i+1] * np.cos(deg2rad*(i+1)*(x-p[2*i+2])) \
            for i in range(len(p)/2)]) + p[0]

def getHarmonicFitParams(x, y, l, sigmay=None):

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
    if sigmay is not None:
        chi2 = (1. / (len(y)-ndof)) * sum((y - fitVals)**2 / sigmay**2)
    else:
        chi2 = (1. / (len(y)-ndof)) * sum((y - fitVals)**2)
    perr = np.sqrt(np.diag(pcov))

    return popt, perr, chi2


def getProjDipole(relint, relerr=None, **opts):

    # Project the relint map in RA
    if relerr is None:
        ra, ri, ra_err = getRIRAProj(relint, relerr, **opts)
    else:
        ra, ri, ra_err, ri_err = getRIRAProj(relint, relerr, **opts)
    # Fit projected RI to cos harmonic functions
    if relerr is None:
        popt, perr, chi2 = getHarmonicFitParams(ra, ri, opts['lmax'])
    else:
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
