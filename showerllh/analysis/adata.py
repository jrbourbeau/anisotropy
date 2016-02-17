#!/usr/bin/env python

##=============================================================================
## Designed to be similar to data.py
##  - modules should bin in fine histograms instead of calculating actual
## median, variance, etc
##  - modules should run over multiple npy files
##=============================================================================

import numpy as np
import matplotlib.pyplot as plt
import glob, os

import myGlobals as my
from useful import getMids
from goodrunlist.goodRunReader import getRunTime

from load_data import load_data
from load_sim import load_sim
from llhtools import *
import bakh, sim, powerFit
from prob import getProbs
from eff import getEff
#from duration import duration
from unfold import unfold
from zfix import zfix


""" Generic function for producing finer bins for histograms """
def fineBins(bins, mult=100):
    min, max = bins.min(), bins.max()
    nbins = len(bins) - 1
    return np.linspace(min, max, mult*nbins+1)


def getDataList(config='IT81', bintype='logdist'):

    my.setupShowerLLH(verbose=False)
    files = glob.glob('%s/%s_data/*_%s.npy' % (my.llh_data, config, bintype))
    if config == 'all':
        files = glob.glob('%s/*_data/*_%s.npy' % (my.llh_data, bintype))
    files.sort()

    configs = [f.split('/')[-2].split('_')[0] for f in files]
    dates = [os.path.basename(f).split('_')[1] for f in files]

    return zip(configs, dates)


""" Energy distribution """
def distro(config, bintype='logdist', cut='llh', xaxis='energy', weight=False):

    # General setup
    labelDict = {'energy':r'$\log_{10}(E/\mathrm{GeV})$',
                 'zenith':r'$\cos(\theta)$',
                 'core':'Distance from center (m)'}
    binDict = {'energy':getEbins(),
               'zenith':np.linspace(0.8, 1, 41),
               'core':np.linspace(0, 700, 71)}

    dataList = getDataList(config, bintype)
    bins = binDict[xaxis]
    xlabel = labelDict[xaxis]
    #fbins = fineBins(bins)

    # Build histograms of desired information
    for cfg, date in dataList[:1]:

        d = load_data(cfg, date, bintype)
        c0 = d['cuts'][cut]
        w = d['weights'][c0] if weight else None
        if xaxis == 'energy':
            y = np.log10(d['ML_energy'])
        if xaxis == 'zenith':
            y = np.cos(d['zenith'])
        if xaxis == 'core':
            y = np.sqrt(d['ML_x']**2 + d['ML_y']**2)

        counts = np.histogram(y[c0], bins=bins, weights=w)[0]
        try: h += counts
        except NameError:
            h = counts

    # Plot
    fig, ax = plt.subplots()
    x = getMids(bins)
    width = bins[1] - bins[0]
    ax.plot(x, h, ls='steps')
    ax.set_xlabel(xlabel)
    ax.set_ylabel('Counts')
    ax.set_yscale('log')

    plt.show()


## Plot the pure number of counts before unfolding ##
def counts(config, cut='llh', bintype='logdist', weight=False, zcorrect=False):

    dataList = getDataList(config, bintype)
    bins = getEbins(reco=True)

    # Build histograms of desired information
    N, Err = {},{}
    for cfg, date in dataList[:2]:

        d = load_data(cfg, date, bintype)
        eList = getComps(d)
        c0 = d['cuts'][cut]
        r = np.log10(d['ML_energy'])
        if zcorrect:
            r -= zfix(d['zenith'], bintype=bintype)

        # Total counts
        w = d['weights'][c0] if weight else None
        w2 = d['weights'][c0]**2 if weight else None
        counts = np.histogram(r[c0], bins=bins, weights=w)[0]
        errors = np.sqrt(np.histogram(r[c0], bins=bins, weights=w2)[0])
        try:
            N['All']   += counts
            Err['All'] += errors
        except KeyError:
            N['All']   = counts
            Err['All'] = errors

        # Counts by composition
        for e in eList:
            ecut = d['llh_comp'] == e
            c1 = c0 * ecut
            w = d['weights'][c1] if weight else None
            w2 = d['weights'][c1]**2 if weight else None
            counts = np.histogram(r[c1], bins=bins, weights=w)[0]
            errors = np.sqrt(np.histogram(r[c1], bins=bins, weights=w2)[0])
            try:
                N[e]   += counts
                Err[e] += errors
            except KeyError:
                N[e]   = counts
                Err[e] = errors

    fig, ax = plt.subplots()
    ax.set_xlabel(r'$\log_{10}(E/\mathrm{GeV})$')
    ax.set_ylabel('Counts')

    # Plot reconstructions
    for e in eList + ['All']:
        pnt = getColor(e) + '.'
        ax.errorbar(getMids(bins), N[e], yerr=Err[e], fmt=pnt, label=e)

    ax.set_yscale('log')
    ax.legend(loc='lower left')

    plt.show()


""" Get Bakhtiyar's spectrum from txt file """
def bakhPlot():

    values = [[] for i in range(8)]

    with open('bakh.txt', 'r') as f:
        for line in f:
            st = 0
            for i in range(8):
                end = line.find(' ',st)
                values[i].append(line[st:end])
                st = end + 1

    for i in range(len(values)):
        for j in range(len(values[i])):
            values[i][j] = float(values[i][j])

    for i in range(len(values)):
        values[i] = np.asarray(values[i])
    stbin, endbin, num, flux, stat, sysup, sysdn, mult = values

    mids = (stbin + endbin)/2.0
    y = {}
    y['flux'] = flux
    errup, errdn = stat+sysup, stat+sysdn
    y['relup'], y['reldn'] = errup/flux, errdn/flux

    y['flux'] *= 10**(-1*mult)
    y['flux'] /= 10**mids

    return mids, y['flux'], y['relup'], y['reldn']


## Convert the number of particles into a flux ##
def NumToFlux(counts, effarea, t):

    Ebins = getEbins(reco=True)
    angle = 0.4*np.pi
    ebin_width = (10**(Ebins[1:]) - 10**(Ebins[:-1]))
    with np.errstate(divide='ignore'):
        flux = counts / (ebin_width * effarea * angle * t)
    return flux


## Plot the flux after unfolding ##
def flux(config, niter, spec=0, smooth=True, zcorrect=False, weight=False):

    # Starting information
    bintype = 'logdist'
    dataList = getDataList(config, bintype)
    dataList = dataList[:2]
    Ebins = getEbins(reco=True)
    Emids = getMids(Ebins)
    cutName = 'llh'
    scale = (10**Emids)**spec

    # Load simulation information
    s = load_sim(bintype=bintype)
    # Relative errors
    effarea, sigma, relerr = getEff(s, s['cuts']['llh'])
    #effarea, sigma, relerr = getEff(s, scut, smooth=smooth)
    # Due to unfolding
    #fl = open('unfold_err.pkl', 'rb')
    #chi2, unrel = pickle.load(fl)
    #fl.close()


    # Get detector runtime
    configs = sorted(list(set([data[0] for data in dataList])))
    t = 0.
    for cfg in configs:
        dates = sorted([data[1] for data in dataList if data[0]==cfg])
        mindate = int(dates[0] + '01')
        next = str(int(dates[-1]) + 1)
        if dates[-1][-2:] == '12':
            next = str(int(date) + 100)
        maxdate = int(next + '01')
        t += getRunTime(cfg, minDate=mindate, maxDate=maxdate)

    # Build histograms of desired information
    N, Err = {},{}
    for cfg, date in dataList:

        d = load_data(cfg, date, bintype)
        eList = getComps(d)
        dcut = d['cuts'][cutName]
        r = np.log10(d['ML_energy'])
        if zcorrect:
            r -= zfix(d['zenith'], bintype='logdist')

        # Create starting arrays
        w = d['weights'][dcut] if weight else None
        #N_passed = float(np.histogram(r[dcut], bins=Ebins, weights=w)[0].sum())
        for e in eList:
            ecut = d['llh_comp'] == e
            w = d['weights'][dcut*ecut] if weight else None
            try: N[e] += np.histogram(r[dcut*ecut], bins=Ebins, weights=w)[0]
            except KeyError:
                N[e] = np.histogram(r[dcut*ecut], bins=Ebins, weights=w)[0]

    # Get probabilities for unfolding
    p = getProbs(s, s['cuts'][cutName], zcorrect=zcorrect)

    Nun, Rel, Flux, Err = {},{},{},{}
    # Get starting probabilities
    N_passed = float(np.sum([N[e] for e in eList]))
    for e in eList:
        p['R'+e] = N[e] / N_passed
        p['T'+e] = powerFit.powerFit(-2.7)

    # Setup plot
    fig, ax = plt.subplots()
    ax.set_title('Energy Spectum using '+cutName+' cut')
    ax.set_xlabel('Log10(Energy/GeV)')
    ax.set_ylabel('Flux (counts/(m^2 s ster GeV)) * E^' + str(spec))

    # Bayesian unfolding
    for i in range(niter):
        p = unfold(p)
        Nun['All'] = np.sum([p['T'+e] for e in eList], axis=0) * N_passed
        #Rel['All'] = np.sqrt(1/Nun['All'] + relerr**2 + unrel[i]**2)
        with np.errstate(divide='ignore'):
            Rel['All'] = np.sqrt(1/Nun['All'] + relerr**2)
        Flux['All'] = NumToFlux(Nun['All'], effarea, t) * scale
        Err['All']  = Flux['All'] * Rel['All']
        #ax.errorbar(Emids, Flux['All'], yerr=Err['All'], fmt='k.', label='Unfolded')
        if i < niter-1:
            for e in eList:
                p['T'+e] = smoother(p['T'+e])

    # Find bin values and errors
    for e in eList:
        Nun[e]  = p['T'+e] * N_passed
        Rel[e]  = np.sqrt(1/Nun[e] + relerr**2)
        Flux[e] = NumToFlux(Nun[e], effarea, t) * scale
        Err[e]  = Flux[e] * Rel[e]

    # Plot
    for e in eList + ['All']:
        pnt = getColor(e) + '.'
        ax.errorbar(Emids, Flux[e], yerr=Err[e], fmt=pnt, label=e)

    # plot original (not unfolded) spectrum
    O_N = np.sum([N[e] for e in eList], axis=0)
    O_relerr = np.sqrt(1/O_N + relerr**2)
    O_flux = NumToFlux(O_N, effarea, t) * scale
    O_err = O_flux * O_relerr
    ax.errorbar(Emids, O_flux, yerr=O_err, fmt='kx', label='Orig')

    # plot Bakhtiyar's spectrum
    B_mids, B_flux, B_relup, B_reldn = bakhPlot()
    B_flux *= ((10**B_mids)**spec)
    B_errup = (B_flux * B_relup)
    B_errdn = (B_flux * B_reldn)
    ax.errorbar(B_mids, B_flux, yerr=(B_errup, B_errdn), fmt='gx', label='Bakhtiyar')

    # plot IT26 spectrum
    #IT26_data = bakh.points['unfolded_twocomponent_th0-110412-shifted']
    #IT26_relerr = bakh.geterrorbars('unfolded_twocomponent_th0-110412-shifted')
    #IT26_mids = np.log10(IT26_data['E'])
    #IT26_flux = np.asarray(IT26_data['dN/dE'])
    #IT26_flux *= ((10**IT26_mids)**spec)
    #IT26_err = IT26_flux * IT26_relerr
    #ax.errorbar(IT26_mids, IT26_flux, yerr=IT26_err, fmt='rx', label='IT-26 Two-Component')

    ax.set_yscale('log')
    ax.legend(loc='lower left')
    #ax.set_xlim((6, 9.5))
    #ax.set_ylim((10**(3), 10**(5)))
    #ax.set_ylim((10**(-22), 10**(-10)))

    #plt.savefig('collab/pics/test.png')
    plt.show()


