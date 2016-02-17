#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import os, argparse

from icecube import photospline
from icecube.photospline import spglam as glam
from icecube.photospline import splinefitstable

from llhtools import getEbins
from load_sim import load_sim

import myGlobals as my
from useful import getMedian, getMids


def zfix(z, bintype='logdist'):

    my.setupShowerLLH(verbose=False)
    inFile = '%s/IT73_sim/Zfix_%s.fits' % (my.llh_data, bintype)
    tab = photospline.I3SplineTable(inFile)
    fit = np.array([tab.eval(cosz) for cosz in np.cos(z)])
    return fit


def makeTable(bintype='logdist', plot=False):

    # Starting parameters
    my.setupShowerLLH(verbose=False)
    s = load_sim(bintype=bintype)
    outFile = '%s/IT73_sim/Zfix_%s.fits' % (my.llh_data, bintype)
    nbins = 100
    thetaMax = 40.
    minE = 4.0

    thetaMax *= np.pi / 180.
    zbins = np.linspace(1, np.cos(thetaMax), nbins+1)[::-1]

    ebins = getEbins()
    t = np.log10(s['MC_energy'])
    r = np.log10(s['ML_energy'])
    z = np.pi - s['zenith']

    # Calculate cut values
    c0 = s['cuts']['llh']
    ecut = r >= minE
    c0 *= ecut

    # Store median and standard deviation info
    x1 = np.cos(z)[c0]
    if x1.min() > zbins.min():
        zbins = zbins[zbins >= x1.min()]
    y = (r - t)[c0]
    medians, sigL, sigR, vars = getMedian(x1, y, zbins)

    w = 1/vars
    nknots = 30
    step_scale = 2/5.
    step = (zbins.max() - zbins.min()) * step_scale
    mids = (zbins[1:] + zbins[:-1]) / 2.
    axes = [mids]
    knots = [np.linspace(zbins.min()-step, zbins.max()+step, nknots)]

    tab = glam.fit(medians, w, axes, knots, order=(4), penalties={2:1e-4})
    if os.path.exists(outFile):
        os.remove(outFile)
    splinefitstable.write(tab, outFile)

    if plot:
        fig, ax = plt.subplots()
        ax.set_title('Energy Resolution vs Reconstructed Zenith', fontsize=18)
        ax.set_xlabel('Cos(zenith)', fontsize=16)
        ax.set_ylabel('Ereco - Etrue (median)', fontsize=16)
        lw = 2
        ms = 7*lw
        pltParams = dict(fmt='.', lw=lw, ms=ms)
        # Energy resolution vs zenith
        x = getMids(zbins)
        ax.errorbar(x, medians, yerr=(sigL,sigR), **pltParams)
        # Spline fit
        fitx = np.linspace(x.min(), x.max(), len(x)*3)
        fit = glam.grideval(tab, [fitx])
        ax.plot(fitx, fit)
        plt.show()



if __name__ == "__main__":

    p = argparse.ArgumentParser(description='Creates spline fit for zenith')
    p.add_argument('-b', '--bintype', dest='bintype',
            default='logdist',
            help='Likelihood binning type')
    p.add_argument('-p', '--plot', dest='plot',
            default=False, action='store_true',
            help='Plot histogram with spline fit')
    args = p.parse_args()

    makeTable(args.bintype, args.plot)

