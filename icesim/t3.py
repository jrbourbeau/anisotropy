#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import glob, os

import sys
sys.path.append('/home/fmcnally/useful')
from useful import getMids

def original():

    files = glob.glob('IC*_hists.npy')
    files.sort()

    for file in files:
        fig, ax = plt.subplots()
        f = np.load(file)
        ebins = np.arange(2.75, 9.01, 0.05)
        x = getMids(ebins)
        for i, y in enumerate(f):
            ntot = float(y.sum())
            ax.step(x, y/ntot, label=i)

        ax.set_title(os.path.basename(file)[:-4])
        ax.set_ylabel('Fraction of Events')
        ax.set_xlabel(r'$\log_{10}(E/\mathrm{GeV})$')
        #ax.legend(loc='upper right')
        #ax.set_yscale('log')

        plt.show()

def new():

    fig, ax = plt.subplots()
    f = np.load('IC59_hists_test.npy')
    f = f.item()
    zbins = np.linspace(0, 1, 151)
    x = getMids(zbins)
    for i, y in enumerate(f['zenith']):
        ntot = float(y.sum())
        ax.step(x, y/ntot, label=i)

    ax.set_ylabel('Fraction of Events')
    ax.set_ylabel('Cos(zenith)')
    ax.legend(loc='upper left')

    plt.show()

if __name__ == "__main__":

    new()
