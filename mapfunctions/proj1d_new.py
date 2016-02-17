#!/usr/bin/env python

from anisotropy.mapfunctions.mapFunctions import getMap
from useful import getMids

import numpy as np
import healpy as hp
import argparse

import matplotlib.pyplot as plt

if __name__ == "__main__":

    p = argparse.ArgumentParser()
    p.add_argument('-d', '--decmin', dest='decmin', type=float,
            help='minimum declination')
    p.add_argument('-D', '--decmax', dest='decmax', type=float,
            help='maximum declination')
    p.add_argument('-n', '--nbins', dest='nbins', type=int,
            default=24,
            help='number of bins in right ascension')
    p.add_argument('-f', '--file', dest='file',
            help='Data file to be analyzed')

    args = p.parse_args()
    opts = vars(args).copy()
    opts['mapName'] = 'relint'

    if not args.file:
        raise SystemExit('\nFilename not entered.\n')

    if not args.decmax:
        opts['mask'] = True

    deg2rad = np.pi/180.

    fig, ax = plt.subplots()
    # Get relative intensity map
    relint = getMap(args.file, **opts)
    # Setup right ascension bins
    rabins = np.linspace(0, 2*np.pi, args.nbins+1)
    print('rabins = {}'.format(rabins))
    # Calculate phi for each pixel
    npix = len(relint)
    nside = hp.npix2nside(npix)
    theta, phi = hp.pix2ang(nside, range(npix))
    # Bin in right ascension
    phiBins = np.digitize(phi, rabins) - 1
    # UNSEEN cut
    cut = (relint != hp.UNSEEN)

    x = getMids(rabins) / deg2rad
    y, yerr = np.zeros((2, len(x)))
    for i in range(len(x)):
        phiCut = (phiBins == i)
        c0 = cut * phiCut
        y[i] = np.mean(relint[c0])
        yerr[i] = np.sqrt(np.var(relint[c0]))

    plt.errorbar(x, y, yerr=yerr, fmt='.')
    tPars = {'fontsize':16}
    ax.set_xlabel(r'Right Ascension', **tPars)
    ax.set_ylabel(r'Relative Intensity',**tPars)
    ax.set_xlim(0.,360.)
    ax.invert_xaxis()

    plt.show()
