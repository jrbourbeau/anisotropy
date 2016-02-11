#!/usr/bin/env python

################################################################################
# GLAM B-spline fit to the EBL optical depth using Gilmore et al., 2009.
################################################################################

import numpy as np
np.set_printoptions(threshold='nan')
import matplotlib as mpl
import matplotlib.pyplot as plt
import os, argparse

import simFunctions_IC as simFunctions
from plotFunctions import cmap_discretize
from analysis import getSimBase, getEbins

from icecube.photospline import spglam as glam
from icecube.photospline import splinefitstable


def spline(config, inFile, outFile, plot=False):

    # Load input median file
    d = np.load(inFile)
    d = d.item()
    xbins, ybins = d['xbins'], d['ybins']
    energies = d['medians']
    vars = d['var']
    vars[vars<1e-15] = np.inf   # Variances sufficiently close to 0 aren't real
    w = 1/vars
    emin, emax = energies[energies!=0].min(), energies.max()
    ebins = [emin] + getEbins() + [emax]

    axes, knots = [],[]
    binList = [xbins, ybins]
    nknots = 30
    step_scale = 2/5.
    for bins in binList:
        mids = (bins[1:]+bins[:-1])/2.
        axes += [mids]
        step = (bins.max() - bins.min()) * step_scale
        knots += [np.linspace(bins.min()-step, bins.max()+step, nknots)]

    tab = glam.fit(energies, w, axes, knots, order=(4), \
            penalties={2:1e-3})

    if plot:

        # Look at spline fit with finer binning
        fitaxes = [[],[]]
        fitaxes[0] = np.linspace(xbins.min(),xbins.max(),len(xbins)*3)
        fitaxes[1] = np.linspace(ybins.min(),ybins.max(),len(ybins)*3)

        fit = glam.grideval(tab, fitaxes)
        #err_fit = 10**glam.grideval(err_tab,axes)

        fig = plt.figure(figsize=(17,6))
        mpl.rc("font", family="serif")
        X, Y = np.meshgrid(axes[0], axes[1])
        fitX, fitY = np.meshgrid(fitaxes[0], fitaxes[1])

        # Setup custom colormap
        cmap = plt.cm.jet
        cmap = cmap_discretize(cmap, ebins)
        cmap.set_under('white')

        ax = fig.add_subplot(121)
        p = ax.pcolor(fitX, fitY, fit.T, cmap=cmap, vmin=emin, vmax=emax)
        cb = fig.colorbar(p, ax=ax)
        ax.set_title('Median Energy')
        ax.set_xlabel('cos(zenith)')
        ax.set_ylabel('log10(Nchannel)')
        ax.set_xlim(xbins.min(), xbins.max())
        ax.set_ylim(ybins.min(), ybins.max())

        ax = fig.add_subplot(122)
        p = ax.pcolor(X, Y, energies.T, cmap=cmap, vmin=emin, vmax=emax)
        ax.set_title('Median Energy')
        ax.set_xlabel('cos(zenith)')
        ax.set_ylabel('log10(Nchannel)')
        ax.set_xlim(xbins.min(), xbins.max())
        ax.set_ylim(ybins.min(), ybins.max())
        plt.show()

    if outFile:
        if os.path.exists(outFile):
            os.remove(outFile)
        splinefitstable.write(tab, outFile)



if __name__ == "__main__":

    p = argparse.ArgumentParser(
            description='Created splined tables for IC energy reconstruction')
    p.add_argument('-c', '--config', dest='config',
            help='Detector configuration')
    p.add_argument('-i', '--inFile', dest='inFile',
            help='Option for input median file to read')
    p.add_argument('-o', '--out', dest='out',
            default=False, action='store_true',
            help='Option to produce output file')
    p.add_argument('-p', '--plot', dest='plot',
            default=False, action='store_true',
            help='Option to plot table and splined version')
    args = p.parse_args()

    simBase = getSimBase(args.config)
    if not args.inFile:
        args.inFile = '%s_median.npy' % simBase
    if args.out:
        args.out = '%s_spline.fits' % simBase

    spline(args.config, args.inFile, args.out, plot=args.plot)




