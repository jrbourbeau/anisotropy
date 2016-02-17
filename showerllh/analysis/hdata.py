#!/usr/bin/env python

##=============================================================================
## Designed to be similar to data.py
##  - modules use pre-made histograms for speed
##  - for info on building histograms see showerllh/run_data
##=============================================================================

import numpy as np
import matplotlib.pyplot as plt
import glob
from copy import copy
from scipy.special import gammaln

import myGlobals as my
from useful import getMids
from goodrunlist.goodRunReader import getRunTime

from llhtools import *
from data import bakhPlot
import powerFit
from prob import getProbs_fast
from eff import getEff, getEff_fast
from unfold import unfold

def bayes2(x1, x2):

    # Calculate totals
    n1 = np.sum(x1)
    n2 = np.sum(x2)
    output = (sum(gammaln(x1+1) + gammaln(x2+1) - gammaln(x1+x2+2))
            + gammaln(n1+n2+2) - gammaln(n1+1) - gammaln(n2+1))

    return output


def loadHists(config=None):

    # Setup global paths
    my.setupShowerLLH(verbose=False)

    # Load existing information
    histfiles = glob.glob('%s/*_data/IT*_hists.npy' % my.llh_data)
    if config != None:
        histfiles = glob.glob('%s/%s_data/hists/*.npy' % (my.llh_data, config))

    h = {}
    for histfile in histfiles:
        htemp = np.load(histfile)
        htemp = htemp.item()
        for key in htemp.keys():
            h[key] = htemp[key]

    return h


def histReader(h, x=None, e=False,
        configs=None, months=None,
        decmin=None, decmax=None,
        ramin=None, ramax=None,
        w=False, z=False, err=False):

    # Require x-axis name
    if x == None:
        raise SystemError('Requires x-value [energy|zenith|core]')

    # Default values for optional cut parameters
    if configs == None:
        configs = ['IT73','IT81','IT81-II','IT81-III','IT81-IV']
    if months == None:
        months = ['%02i' % i for i in range(1,13)]
    decmin = 0 if decmin==None else decmin
    decmax = 180 if decmax==None else decmax
    ramin = 0 if ramin==None else ramin
    ramax = 360 if ramax==None else ramax

    # Fixed binning
    binDict = {'energy':getEbins(reco=True), 'zenith':np.linspace(0.8, 1, 81),
            'core':np.linspace(0, 700, 141)}

    # Build list of desired parameters
    myParams = []
    myParams += [x]
    if w: myParams += ['w']
    if z: myParams += ['z']
    if err: myParams += ['err']

    st, end = 2, -3
    if e:
        myParams += [e]
        end = -2

    # Reduce keylist to desired keys
    keys = [k for k in h.keys() if k.split('_')[st:end]==myParams]
    keys = [k for k in keys if any([c in k.split('_') for c in configs])]
    keys = [k for k in keys if any([m in k.split('_')[1][4:6] for m in months])]
    keys = [k for k in keys if int(k.split('_')[-2].split('-')[0]) >= decmin]
    keys = [k for k in keys if int(k.split('_')[-2].split('-')[-1]) <= decmax]
    keys = [k for k in keys if int(k.split('_')[-1].split('-')[0]) >= ramin]
    keys = [k for k in keys if int(k.split('_')[-1].split('-')[-1]) <= ramax]

    hvalues = np.sum([h[k] for k in keys], axis=0)
    bins = binDict[x]

    if err:
        hvalues = np.sqrt(hvalues)

    return hvalues, bins


""" Energy distribution """
def distro(configs=None, xaxis='energy', weight=False, zcorrect=False):

    # General setup
    labelDict = {'energy':r'$\log_{10}(E/\mathrm{GeV})$',
                 'zenith':r'$\cos(\theta)$',
                 'core':'Distance from center (m)'}

    h = loadHists()
    eList = ['p','h','o','f']
    if configs == None:
        configs = sorted(list(set([k.split('_')[0] for k in h.keys()])))

    # Total counts
    counts, bins = histReader(h, xaxis, configs=configs, 
            w=weight, z=zcorrect)
    xlabel = labelDict[xaxis]

    if xaxis == 'core':
        binArea = np.pi * (bins[1:]**2 - bins[:-1]**2)
        counts /= binArea

    # Plot
    fig, ax = plt.subplots()
    x = getMids(bins)
    width = bins[1] - bins[0]
    ax.plot(x, counts, ls='steps')
    ax.set_xlabel(xlabel)
    ax.set_ylabel('Counts')
    ax.set_yscale('log')

    plt.show()


## Plot the pure number of counts before unfolding ##
def counts(configs=None, weight=False, zcorrect=False):

    h = loadHists()
    eList = ['p','h','o','f']
    if configs == None:
        configs = sorted(list(set([k.split('_')[0] for k in h.keys()])))
    hParams = {'configs':configs, 'w':weight, 'z':zcorrect}

    # Build histograms of desired information
    N, Err = {},{}

    # Total counts
    N['All'], bins = histReader(h, 'energy', **hParams)
    Err['All'], bins = histReader(h, 'energy', err=True, **hParams)

    # Counts by composition
    for e in eList:
        N[e], bins = histReader(h, 'energy', e=e, **hParams)
        Err[e], bins = histReader(h, 'energy', e=e, err=True, **hParams)

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


def myEnergy():

    h = loadHists()
    configs = sorted(list(set([k.split('_')[0] for k in h.keys()])))
    hParams = {'configs':configs, 'w':False, 'z':True}
    N = {}
    for e in ['p','h','o','f']:
        N[e], bins = histReader(h, 'energy', e=e, **hParams)
        emids = getMids(bins)
        ecut = (emids >= 6.2)
        N[e] = N[e][ecut]
        emids = emids[ecut]
        print e, N[e].sum()
        ntemp, i = 0, 0
        check = False
        print N[e]
        while not check:
            ntemp += N[e][i]
            if ntemp >= sum(N[e])/2:
                check = True
                eval = emids[i]
            i += 1
        print eval


def myCounts(spec=0, out=None):

    h = loadHists()
    configs = sorted(list(set([k.split('_')[0] for k in h.keys()])))
    hParams = {'configs':configs, 'w':False, 'z':True}
    tParams = {'fontsize':14}

    # Minimum energy
    emin = 6.2

    fig, ax = plt.subplots()
    ax.set_xlabel(r'$\log_{10}(E/\mathrm{GeV})$', **tParams)
    ax.set_ylabel('Normalized Counts', **tParams)

    # Spectrum for deficit (60 - 120)
    ND, bins = histReader(h, 'energy', ramin=60, ramax=120, **hParams)
    emids = getMids(bins)
    ecut = (emids >= emin)
    emids = emids[ecut]
    scale = (10**emids)**spec
    ND = ND[ecut]
    Err = np.sqrt(ND)
    Err[ND == 0] = 2.3
    ntot = float(ND.sum())
    NDN = ND * scale / ntot
    Err = Err * scale / ntot
    # Custom errorbars
    Errmin = copy(Err)
    for i in range(len(Errmin)):
        if NDN[i] - Errmin[i] < 1e-8:
            Errmin[i] = NDN[i] - 1e-9
    ax.errorbar(emids, NDN, yerr=(Errmin, Err), fmt='.', label='Deficit')

    # Spectrum for other (0 - 60 + 120 - 360)
    N0, bins = histReader(h, 'energy', ramin=0, ramax=60, **hParams)
    N1, bins = histReader(h, 'energy', ramin=120, ramax=360, **hParams)
    N = (N0 + N1)[ecut]
    Err = np.sqrt(N)
    Err[N == 0] = 2.3
    ntot = float(N.sum())
    NN = N * scale / ntot
    Err = Err * scale / ntot
    # Custom errorbars
    Errmin = copy(Err)
    for i in range(len(Errmin)):
        if NN[i] - Errmin[i] < 1e-8:
            Errmin[i] = NN[i] - 1e-9
    ax.errorbar(emids, NN, yerr=(Errmin, Err), fmt='.', label='Rest of sky')

    if spec == 0:
        ax.set_ylim((1e-8, 1))
    ax.set_yscale('log')
    ax.legend()
    print 'Bayes factor:', bayes2(NDN, NN)
    if out != None:
        plt.savefig(out, dpi=300, bbox_inches='tight')
    plt.show()


#def effCorrect()

## Convert the number of particles into a flux ##
def NumToFlux(counts, effarea, t):

    Ebins = getEbins(reco=True)
    angle = 0.4*np.pi
    ebin_width = (10**(Ebins[1:]) - 10**(Ebins[:-1]))
    with np.errstate(divide='ignore'):
        flux = counts / (ebin_width * effarea * angle * t)
    return flux


def split_flux(splitName, **kwargs):

    defaultArgs = {'niter':5, 'spec':0, 'smooth':True, 'zcorrect':True,
                   'weight':True, 'orig':False, 'bakh':False, 'comps':False,
                   'emin':None, 'linear':False,
                   'months':None, 'configs':None, 'ramin':None, 'ramax':None}
    if splitName in defaultArgs:
        defaultArgs.pop(splitName)
    if splitName == 'dec':
        defaultArgs.pop('decmin')
        defaultArgs.pop('decmax')
    if splitName == 'ra':
        defaultArgs.pop('ramin')
        defaultArgs.pop('ramax')
    for k in defaultArgs:
        if k not in kwargs:
            kwargs[k] = defaultArgs[k]

    fig, ax = plt.subplots()

    m = [['12','01','02'],['03','04','05'],['06','07','08'],['09','10','11']]
    #m = [['%02i' % month] for month in range(1,13)]
    #m = [['%02i' % month] for month in range(1,7)]
    #m = [['%02i' % month] for month in range(7,13)]
    if splitName == 'months':
        for months in m:
            label = '-'.join([months[0], months[-1]])
            flux(months=months, tax=ax, tlabel=label, **kwargs)

    if splitName == 'configs':
        for cfg in ['IT73','IT81','IT81-II','IT81-III','IT81-IV']:
            flux(configs=[cfg], tax=ax, tlabel=cfg, **kwargs)

    decbins = [0,12,24,40]
    if splitName == 'dec':
        for i in range(len(decbins) - 1):
            decmin, decmax = decbins[i], decbins[i+1]
            label = '%s-%s' % (decmin, decmax)
            flux(decmin=decmin, decmax=decmax, tax=ax, tlabel=label, **kwargs)

    rabins = range(0,361,60)
    if splitName == 'ra':
        for i in range(len(rabins) - 1):
            ramin, ramax = rabins[i], rabins[i+1]
            label = '%s-%s' % (ramin, ramax)
            flux(ramin=ramin, ramax=ramax, tax=ax, tlabel=label, **kwargs)

    ax.set_title('Energy Spectum')
    ax.set_xlabel('Log10(Energy/GeV)')
    ax.set_ylabel('Flux (counts/(m^2 s ster GeV)) * E^' + str(kwargs['spec']))
    ax.set_yscale('log')
    #ax.legend(loc='lower left')
    if kwargs['emin']:
        ax.set_xlim((kwargs['emin'], 9.5))

    plt.show()


## Plot the flux after unfolding ##
def flux(niter=5, configs=None, spec=0, smooth=True, 
        zcorrect=False, weight=False,
        orig=True, bakh=True,
        emin=None, linear=False,
        tax=None, tlabel=None,
        comps=True, months=None,
        decmin=None, decmax=None,
        ramin=None, ramax=None):

    # Starting information
    Ebins = getEbins(reco=True)
    Emids = getMids(Ebins)
    scale = (10**Emids)**spec
    h = loadHists()
    hParams = {'configs':configs, 'months':months, 
               'decmin':decmin, 'decmax':decmax,
               'ramin':ramin, 'ramax':ramax,
               'w':weight, 'z':zcorrect}
    eList = ['p','h','o','f']
    if configs == None:
        configs = sorted(list(set([k.split('_')[0] for k in h.keys()])))

    # Load simulation information
    effarea, sigma, relerr = getEff_fast(linear)    # Effective area
    p = getProbs_fast(zcorrect=zcorrect)            # Probs for unfolding
    # Relative error due to unfolding...
    #fl = open('unfold_err.pkl', 'rb')
    #chi2, unrel = pickle.load(fl)
    #fl.close()

    # Get detector runtime
    t = 0.
    for cfg in configs:
        t += getRunTime(cfg, months=months)

    # Total counts
    N, Err = {},{}
    N['All'], bins = histReader(h, x='energy', **hParams)
    Err['All'], bins = histReader(h, x='energy', err=True, **hParams)
    # Counts by composition
    for e in eList:
        N[e], bins = histReader(h, x='energy', e=e, **hParams)
        Err[e], bins = histReader(h, x='energy', e=e, err=True, **hParams)

    # Get starting probabilities
    Nun, Rel, Flux, Err = {},{},{},{}
    N_passed = float(np.sum([N[e] for e in eList]))
    for e in eList:
        p['R'+e] = N[e] / N_passed
        p['T'+e] = powerFit.powerFit(-2.7)

    # Setup plot
    if tax == None:
        fig, ax = plt.subplots()
        ax.set_title('Energy Spectum')
        ax.set_xlabel('Log10(Energy/GeV)')
        ax.set_ylabel('Flux (counts/(m^2 s ster GeV)) * E^' + str(spec))
    else:
        ax = tax

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
    pltList = ['All']
    if comps:
        pltList += eList
    for e in pltList:
        pnt = getColor(e)+'.' if tax==None else '-'
        label = e if tlabel==None else tlabel
        ax.errorbar(Emids, Flux[e], yerr=Err[e], fmt=pnt, label=label)

    # plot original (not unfolded) spectrum
    if orig:
        O_N = np.sum([N[e] for e in eList], axis=0)
        O_relerr = np.sqrt(1/O_N + relerr**2)
        O_flux = NumToFlux(O_N, effarea, t) * scale
        O_err = O_flux * O_relerr
        ax.errorbar(Emids, O_flux, yerr=O_err, fmt='kx', label='Orig')

    # plot Bakhtiyar's spectrum
    if bakh:
        B_mids, B_flux, B_relup, B_reldn = bakhPlot()
        B_flux *= ((10**B_mids)**spec)
        B_errup = (B_flux * B_relup)
        B_errdn = (B_flux * B_reldn)
        ax.errorbar(B_mids, B_flux, yerr=(B_errup, B_errdn), 
                fmt='gx', label='Bakhtiyar')

    # plot IT26 spectrum
    #IT26_data = bakh.points['unfolded_twocomponent_th0-110412-shifted']
    #IT26_relerr = bakh.geterrorbars('unfolded_twocomponent_th0-110412-shifted')
    #IT26_mids = np.log10(IT26_data['E'])
    #IT26_flux = np.asarray(IT26_data['dN/dE'])
    #IT26_flux *= ((10**IT26_mids)**spec)
    #IT26_err = IT26_flux * IT26_relerr
    #ax.errorbar(IT26_mids, IT26_flux, yerr=IT26_err, fmt='rx', label='IT-26 Two-Component')

    if tax == None:
        ax.set_yscale('log')
        ax.legend(loc='lower left')
        if emin:
            ax.set_xlim((emin, 9.5))
        #ax.set_ylim((10**(3), 10**(5)))
        #ax.set_ylim((10**(-22), 10**(-10)))

        #plt.savefig('collab/pics/test.png')
        plt.show()


