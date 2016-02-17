#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

from llhtools import getEbins

def forceAspect(ax, aspect=1):
    im = ax.get_images()
    extent = im[0].get_extent()
    ax.set_aspect(abs((extent[1]-extent[0]) / (extent[3]-extent[2]))/aspect)

""" 2-D histogram of ML energy vs MC energy """
def hist2d(s, cut, log=True):

    # Setup 2D histogram of Ereco (y) vs Etrue (x) for ML reconstruction
    fig, ax = plt.subplots()
    ax.set_title('NStations vs Ereco (actual composition)')
    ax.set_xlabel('NStations')
    ax.set_ylabel('log10(LLH/GeV)')

    ## x translates to the y-axis and vice versa ##
    x = np.log10(s['ML_energy'][cut])
    y = s['NStations'][cut]
    Ebins = getEbins()
    Nbins = range(1, 82)

    h, xedges, yedges = np.histogram2d(x, y, bins=(Ebins, Nbins),
            normed=False, weights=None)
    extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
    if log:
        h = np.log10(h)

    im = ax.imshow(h, extent=extent, origin='lower', interpolation='None')
    fig.colorbar(im)
    forceAspect(ax)
    plt.show()
