#!/usr/bin/env python

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import matplotlib as mpl


def load_sim(config):

    prefix = '/data/user/fmcnally/anisotropy/sim'
    simFile = '%s/%s_sim.npy' % (prefix, config)

    s = np.load(simFile)
    s = s.item()

    return s


def getWeights(s):

    # Weighting modules from offline software
    #from icecube.weighting.fluxes import Hoerandel, Hoerandel_IT
    import icecube.weighting.fluxes as fluxes
    from icecube.weighting.weighting import from_simprod

    simList = sorted(list(set(s['sim'])))
    simList = [int(sim) for sim in simList]
    if 9166 in simList:
        simList[simList.index(9166)] = 9122

    print 'Making generator...'
    generator = from_simprod(simList[0])[1] * 0
    for sim in simList:
        nfiles, gen = from_simprod(sim)
        generator += nfiles * gen

    print 'Calculating weights...'
    flux = fluxes.Hoerandel_IT()
    #flux = fluxes.Hoerandel()
    weights = flux(s['MC_energy'], s['MC_type']) / \
              generator(s['MC_energy'], s['MC_type'])

    return weights


def weighted_percentile(data, weights, percentile):

    data = np.asarray(data)
    weights = np.asarray(weights, np.float)
    percentile *= 0.01

    n = data.shape[0]
    i = np.argsort(data)
    sd = np.take(data, i, axis=0)
    sw = np.take(weights, i, axis=0)
    aw = np.add.accumulate(sw)
    if not aw[-1] > 0:
        raise ValueError, "Nonpositive weight sum"

    w = (aw-0.5*sw) / aw[-1]
    idx = np.searchsorted(w, percentile)

    if idx == 0:
        o = sd[0]
    elif idx == n:
        o = sd[n-1]
    else:
        f1 = (w[idx] - percentile)/(w[idx] - w[idx-1])
        f2 = (percentile - w[idx-1])/(w[idx] - w[idx-1])
        assert f1>=0 and f2>=0 and f1<=1 and f2<=1
        assert abs(f1+f2-1.0) < 1e-6
        o = sd[idx-1]*f1 + sd[idx]*f2

    return o


def readDist(config, nmin, nmax):

    # Read in distribution information
    infoFile = '/data/user/fmcnally/anisotropy/sim/eDist_IT.txt'
    with open(infoFile, 'r') as f:
        lines = f.readlines()
    table = [line.strip().split(' ') for line in lines]

    # Limit to desired data point
    data = [line[1:] for line in table if line[0]==config]
    if config == 'IT':
        data = [line[1:] for line in table if line[0][:2]==config]
    nmin, nmax = str(nmin), str(nmax)
    data = np.array([l for l in data if l[:2]==[nmin,nmax]], dtype='float')

    median = data[:,2]
    sigL   = data[:,3]
    sigR   = data[:,4]

    # Create combined distribution
    ave_median = np.average(median)
    ave_sigL = np.average(sigL)
    ave_sigR = np.average(sigR)

    return ave_median, ave_sigL, ave_sigR


def edist(s, nbins=20):

    c0 = s['nstation'] >= 8
    values  = s['MC_energy'][c0]
    weights = getWeights(s)[c0]

    # Plot the energy distribution
    fig, ax = plt.subplots()
    ax.hist(np.log10(values), bins=nbins, weights=weights, histtype='step')
    plt.show()

    # Print information about the distribution
    # Weighted average
    mean = np.average(values, weights=weights)
    # Biased weighted variance and standard deviation
    var  = np.average((values-mean)**2, weights=weights)
    std  = np.sqrt(var)
    # Unbiased weighted variance and standard deviation
    var = (sum(weights * (values-mean)**2) /
          (sum(weights) - sum(weights**2)/sum(weights)))
    std = np.sqrt(var)
    print 'Weighted Mean & Std. Dev.:'
    print '%.3e +- %.3e' % (mean, std)
    print 'Weighted Median & 68% Containment:'
    median = weighted_percentile(values, weights, 50)
    sigL = median - weighted_percentile(values, weights, 16)
    sigR = weighted_percentile(values, weights, 84) - median
    print '%.3e - %.3e + %.3e' % (median, sigL, sigR)


def plotDist(s, nbins=100):

    c0 = s['nstation'] >= 8
    energies = np.log10(s['MC_energy'][c0])
    weights  = getWeights(s)[c0]
    bins = np.linspace(energies.min(), energies.max(), nbins+1)

    fig, ax = plt.subplots()
    ax.hist(energies, bins=bins, weights=weights, histtype='step', label='Total')
    for comp in ['P','He','O','Fe']:
        c1 = (s['comp'] == comp)[c0]
        ax.hist(energies[c1], bins=bins, weights=weights[c1], histtype='step',
                label=comp)

    ax.set_yscale('log')
    ax.set_ylim(1e-5, 1)
    plt.legend()
    plt.show()

