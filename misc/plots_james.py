#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

from plotFunctions import cmap_discretize
from useful import getMids

def eres2(config, out=False, batch=False):

    # Basic setup
    fig, ax = plt.subplots()
    f = np.load('/home/fmcnally/anisotropy/icesim/%s_hists.npy' % config)
    #lw = 2
    #ms = 7*lw
    #pltParams = {'fmt':'.', 'lw':lw, 'ms':ms}

    # Energy binning information
    ebins = np.arange(2.75, 9.01, 0.05)
    x = getMids(ebins)
    for i, y in enumerate(f):
        ntot = float(y.sum())
        ax.step(x, y/ntot, label=i+1)

    ax.set_xlim(2.75, 8.5)
    tPars = {'fontsize':16}
    ax.set_xlabel(r'True Energy ($\log_{10}(E/\mathrm{GeV})$)', **tPars)
    ax.set_ylabel('Fraction of Events', **tPars)
    #ax.set_title('Energy Distributions for Cuts', fontsize=16)
    plt.legend(loc='upper right')

    if out != False:
        plt.savefig(out, dpi=300, bbox_inches='tight')
    if not batch:
        plt.show()


if __name__ == "__main__":

    eres2('IC59', out=False)
