#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

from llhtools import *
import bakh, sim, powerFit
from prob import getProbs
from eff import getEff
from duration import duration
from unfold import unfold

""" Energy distribution """
def distro(d, cut, xaxis='energy', log=True, returnValues=False):

    if xaxis == 'energy':
        x = np.log10(d['ML_energy'])
        bins = getEbins()
        xlabel = r'$\log_{10}(E/\mathrm{GeV})$'
    if xaxis == 'zenith':
        x = np.cos(d['zenith'])
        zmin = np.cos(d['zenith'][cut].max())
        bins = np.linspace(zmin, 1, 41)
        xlabel = r'$\cos(\theta)$'
    if xaxis == 'core':
        x = np.sqrt(d['ML_x']**2 + d['ML_y']**2)
        bins = np.linspace(0, x[cut].max(), 41)
        xlabel = 'Distance from center (m)'

    fig, ax = plt.subplots()
    counts, bins, bars = ax.hist(x[cut], bins=bins, histtype='step')
    ax.set_xlabel(xlabel)
    ax.set_ylabel('Counts')
    if log:
        ax.set_yscale('log')

    if returnValues:
        return ax, counts, bins

    plt.show()


def test(d):

    plt.ioff()
    ax, counts, bins = distro(d, d['cuts']['llh'], returnValues=True)
    xlabel = ax.get_xlabel
    ylabel = ax.get_ylabel
    plt.close(plt.gcf())
    ax0, counts0, bins0 = distro(d, d['cuts']['llh'], returnValues=True)
    plt.close(plt.gcf())
    if all(bins == bins0):
        counts += counts0
    plt.ion()


## Plot the pure number of counts before unfolding ##
#def counts(d, r, cutName):
def counts(d, cut):

    r = np.log10(d['ML_energy'])
    #if zcorrect:
    #    r = zfix(d)
    w = d['weights']
    Emids = getEmids()
    #cut = d['cuts'][cutName]
    eList = ['p','h','o','f']
    colorDict = {'p':'b', 'h':'y', 'o':'c', 'f':'r', 'All':'k'}

    fig, ax = plt.subplots()
    #ax.set_title('Counts using '+cutName+' cut')
    ax.set_title('Counts using cut')
    ax.set_xlabel('Log10(Energy/GeV)')
    ax.set_ylabel('Counts')

    # Get the counts and error
    N, Err = {},{}
    N['All']   = Nfinder(r, cut, w=w)
    Err['All'] = Efinder(r, cut, w=w)
    for e in eList:
        ecut = d['llh_comp'] == e
        N[e]   = Nfinder(r, cut*ecut, w=w)
        Err[e] = Efinder(r, cut*ecut, w=w)

    # Plot reconstructions
    for e in eList + ['All']:
        pnt = colorDict[e] + '.'
        ax.errorbar(Emids, N[e], yerr=Err[e], fmt=pnt, label=e)

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
def NumToFlux(counts, effarea, t, st):

    Ebins = getEbins()
    angle = 0.4*pi
    ebin_width = (10**(Ebins[1:]) - 10**(Ebins[:-1]))[st:]
    flux = counts / (ebin_width * effarea * angle * t)
    return flux


## Plot the flux after unfolding ##
def flux(d, s, niter, spec=0, smooth=True):

    # Starting information
    r = np.log10(d['ML_energy'])
    t = duration(d)
    Emids = getEmids()
    w = d['weights']
    cutName = 'llh'
    dcut, scut = d['cuts'][cutName], s['cuts'][cutName]
    eList = ['p', 'f']
    colorDict = {'p':'b', 'h':'y', 'o':'c', 'f':'r', 'All':'k'}

    # Setup plot
    f = plt.figure()
    ax1 = f.add_subplot(1,1,1)
    ax1.set_title('Energy Spectum using '+cutName+' cut')
    ax1.set_xlabel('Log10(Energy/GeV)')
    ax1.set_ylabel('Flux (counts/(m^2 s ster GeV)) * E^' + str(spec))

    # Create starting arrays
    N = {}
    N_passed = Nfinder(r, dcut, w=w).sum()
    for e in eList:
        ecut = d['llh_comp'] == e
        N[e] = Nfinder(r, dcut*ecut, w=w)

    # Get probabilities for unfolding
    p = getProbs(s, scut)
    st = len(Emids) - len(p['Rf|Tf'][0])
    Emids = Emids[st:]
    scale = (10**Emids)**spec

    # Relative errors
    effarea, sigma, relerr = getEff(s, scut, smooth=smooth)
    effarea, sigma, relerr = effarea[st:], sigma[st:], relerr[st:]
    # Due to unfolding
    #fl = open('unfold_err.pkl', 'rb')
    #chi2, unrel = pickle.load(fl)
    #fl.close()

    # Get starting probabilities
    for e in eList:
        p['R'+e] = N[e] / N_passed
        p['T'+e] = powerFit.powerFit(-2.7, st=st)

    # Bayesian unfolding
    Nun, Rel, Flux, Err = {},{},{},{}
    seterr(invalid='ignore')
    for i in range(niter):
        p = unfold(p)
        Nun['All'] = sum([p['T'+e] for e in eList], axis=0) * N_passed
        #Rel['All'] = np.sqrt(1/Nun['All'] + relerr**2 + unrel[i]**2)
        Rel['All'] = np.sqrt(1/Nun['All'] + relerr**2)
        Flux['All'] = NumToFlux(Nun['All'], effarea, t, st) * scale
        Err['All']  = Flux['All'] * Rel['All']
        #ax1.errorbar(Emids, Flux['All'], yerr=Err['All'], fmt='k.', label='Unfolded')
        if i < niter-1:
            for e in eList:
                p['T'+e] = smoother(p['T'+e])
    seterr(invalid='warn')

    # Find bin values and errors
    for e in eList:
        Nun[e]  = p['T'+e] * N_passed
        Rel[e]  = np.sqrt(1/Nun[e] + relerr**2)
        Flux[e] = NumToFlux(Nun[e], effarea, t, st) * scale
        Err[e]  = Flux[e] * Rel[e]

    # Plot
    for e in eList + ['All']:
        pnt = colorDict[e] + '.'
        ax1.errorbar(Emids, Flux[e], yerr=Err[e], fmt=pnt, label=e)

    # plot original (not unfolded) spectrum
    O_N = (sum([N[e] for e in eList], axis=0))[st:]
    O_relerr = np.sqrt(1/O_N + relerr**2)
    O_flux = NumToFlux(O_N, effarea, t, st) * scale
    O_err = O_flux * O_relerr
    ax1.errorbar(Emids, O_flux, yerr=O_err, fmt='kx', label='Orig')

    # plot Bakhtiyar's spectrum
    #B_mids, B_flux, B_relup, B_reldn = bakhPlot()
    #B_flux *= ((10**B_mids)**spec)
    #B_errup = (B_flux * B_relup)
    #B_errdn = (B_flux * B_reldn)
    #ax1.errorbar(B_mids, B_flux, yerr=(B_errup, B_errdn), fmt='gx', label='Bakhtiyar')

    # plot IT26 spectrum
    #IT26_data = bakh.points['unfolded_twocomponent_th0-110412-shifted']
    #IT26_relerr = bakh.geterrorbars('unfolded_twocomponent_th0-110412-shifted')
    #IT26_mids = np.log10(IT26_data['E'])
    #IT26_flux = np.asarray(IT26_data['dN/dE'])
    #IT26_flux *= ((10**IT26_mids)**spec)
    #IT26_err = IT26_flux * IT26_relerr
    #ax1.errorbar(IT26_mids, IT26_flux, yerr=IT26_err, fmt='rx', label='IT-26 Two-Component')

    ax1.set_yscale('log')
    ax1.legend(loc='lower left')
    #ax1.set_xlim((6, 9.5))
    #ax1.set_ylim((10**(3), 10**(5)))
    #ax1.set_ylim((10**(-22), 10**(-10)))

    #plt.savefig('collab/pics/test.png')
    plt.show()


