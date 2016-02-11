#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

from plotFunctions import cmap_discretize
from analysis import getSimBase, getEbins, readDist
from useful import getMids


""" Show the histogram used for reconstructing energy """
def reco_energy_plot(config, out=False, batch=False):

    # Load median file information
    simBase = getSimBase(config)
    inFile = '%s_median.npy' % simBase
    d = np.load(inFile)
    d = d.item()
    xbins, ybins, ebins = d['xbins'], d['ybins'], d['ebins']

    # Calculate median energy value from bin
    energies = d['medians']
    energies = energies.T

    # Trim down the top of the grid
    ycts = energies.sum(axis=1)
    yidx = np.where(ycts==0)[0][0] + 2
    ybins = ybins[:yidx]
    energies = energies[:yidx]

    emin, emax = energies[energies!=0].min(), energies[energies!=0].max()
    emin, emax = 3.75, 8
    ebins = [emin] + getEbins() + [emax]

    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111)
    X, Y = np.meshgrid(xbins, ybins)
    cmap = plt.cm.jet
    cmap = cmap_discretize(cmap, ebins)
    cmap.set_under('white')
    tPars = {'fontsize':16}

    p = ax.pcolor(X, Y, energies, cmap=cmap, vmin=emin, vmax=emax)
    cb = fig.colorbar(p, ax=ax, ticks=ebins)
    cb.ax.set_yticklabels(['%.2f' % ebin for ebin in ebins])
    cb.set_label(r'$\mathrm{log}_{10}(E/\mathrm{GeV})$',
            rotation=270, labelpad=20, **tPars)
    #ax.set_title('Median Energy vs Zenith and Nchannel')
    ax.set_xlabel(r'$\mathrm{cos}(\theta_\mathrm{reco})$', **tPars)
    ax.set_ylabel(r'$\mathrm{log}_{10}(N_\mathrm{channel})$', **tPars)
    ax.set_xlim(xbins.min(), xbins.max())
    ax.set_ylim(ybins.min(), ybins.max())
    if out != False:
        #outFile = '%s/%s_Median_Energy' % (outPrefix, config)
        plt.savefig(out, dpi=300, bbox_inches='tight')
    if not batch:
        plt.show()


def eres(config, out=False, batch=False, old=False):

    # Basic setup
    fig, ax = plt.subplots()
    lw = 2
    ms = 7*lw
    pltParams = {'fmt':'.', 'lw':lw, 'ms':ms}

    d = {}
    configs = ['IC59','IC79','IC86','IC86-II','IC86-III','IC86-IV']
    if config in configs:
        configs = [config]

    # Energy binning information
    eList = getEbins() + [100]
    eStrings = ['%s' % float(i) for i in eList]
    ePairs = [[eStrings[i], eStrings[i+1]] for i in range(len(eStrings)-1)]
    minE, maxE = np.asarray(eList[:-1]), np.asarray(eList[1:])
    if maxE[-1] == 100:
        maxE[-1] = maxE[-2] + 1
    emids = (minE + maxE) / 2.
    xL = emids - minE
    xR = maxE - emids
    xR[-1] = 0

    # Read in distribution information
    for cfg in configs:
        d[cfg] = {}
        for key in ['median','sigL','sigR']:
            d[cfg][key] = []
        for emin, emax in ePairs:
            median, sigL, sigR = readDist(cfg, emin, emax, old=old)
            d[cfg]['median'] += [median]
            d[cfg]['sigL']   += [sigL]
            d[cfg]['sigR']   += [sigR]

        if len(configs) > 1:
            pltParams['label'] = cfg
            ax.errorbar(emids, d[cfg]['median'],
                    yerr=[d[cfg]['sigL'], d[cfg]['sigR']], **pltParams)
            emids += .02

        else:
            #ax.errorbar(emids, d[cfg]['median'], xerr=[xL, xR], 
            #        yerr=[d[cfg]['sigL'], d[cfg]['sigR']], **pltParams)
            #ax.arrow(emids[-1], d[cfg]['median'][-1], 0.5, 0.0,
            #        fc='b', ec='b', head_width=.1, head_length=.2)
            ax.errorbar(emids, d[cfg]['median'], 
                    yerr=[d[cfg]['sigL'], d[cfg]['sigR']], **pltParams)

    # Create combined distribution (start unweighted)
    ave_median = np.average([d[key]['median'] for key in d.keys()], axis=0)
    ave_sigL = np.average([d[key]['sigL'] for key in d.keys()], axis=0)
    ave_sigR = np.average([d[key]['sigR'] for key in d.keys()], axis=0)

    if len(configs) > 1:
        ax.errorbar(emids, ave_median, yerr=[ave_sigL, ave_sigR],
                fmt='kx', label='Total')

    ax.set_xlim(3.5, 8)
    tPars = {'fontsize':16}
    #ax.set_xlabel(r'Reconstructed Energy ($log_{10}(E/\mathrm{GeV})$)',
    #        **tPars)
    ax.set_xlabel(r'Energy of Bin Center ($log_{10}(E/\mathrm{GeV})$)',
            **tPars)
    ax.set_ylabel(r'True Energy ($log_{10}(E/\mathrm{GeV})$)', **tPars)
    #ax.set_title('Energy Distributions for Cuts', fontsize=16)
    if len(configs) > 1:
        plt.legend(loc='lower right')

    if out != False:
        plt.savefig(out, dpi=300, bbox_inches='tight')
    if not batch:
        plt.show()


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

    #eres('IC', old=False)
    reco_energy_plot('IC86')
