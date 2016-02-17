#!/usr/bin/env python

################################################################################
## Look at the cut efficiency (requires load_sim.py)                          ##
################################################################################

import numpy as np
from scipy import optimize
from llhtools import *
import matplotlib.pyplot as plt
import time

from useful import getMedian, getMids
degree = np.pi / 180.

##===========================================================================##
## Plotting functions

""" Return the median, upper, and lower 34% containment for an array """
def getMed(x):
    if len(x) == 0:
        x = np.array([0])
    median = np.percentile(x, 50)
    sigma_min = median - np.percentile(x, 16)
    sigma_max = np.percentile(x, 84) - median
    return median, sigma_min, sigma_max


""" Correction for energy resolution vs zenith angle """
def zfix(s, cut='llh', nbins=100, thetaMax=40., minE=5.0, 
        plot=False, out=False):

    from icecube.photospline import spglam as glam
    thetaMax *= degree
    zbins = np.linspace(1, np.cos(thetaMax), nbins+1)[::-1]

    ebins = getEbins()
    t = np.log10(s['MC_energy'])
    r = np.log10(s['ML_energy'])
    z = np.pi - s['zenith']

    # Calculate cut values
    c0 = np.logical_not(np.isnan(r))
    ecut = r >= minE
    c0 *= ecut
    if cut != None:
        c0 *= s['cuts'][cut]

    # Store median and standard deviation info
    x1 = np.cos(z)[c0]
    if x1.min() > zbins.min():
        zbins = zbins[zbins >= x1.min()]
    y = (r - t)[c0]
    medians, sigL, sigR, vars = getMedian(x1, y, zbins)

    w = 1/vars
    nknots = 30
    step_scale = 2/5.
    step = (zbins.max() - zbins.min()) * step_scale
    mids = (zbins[1:] + zbins[:-1]) / 2.
    axes = [mids]
    knots = [np.linspace(zbins.min()-step, zbins.max()+step, nknots)]

    tab = glam.fit(medians, w, axes, knots, order=(4), penalties={2:1e-4})

    if plot:
        fig, ax = plt.subplots()
        #ax.set_title('Energy Resolution vs Reconstructed Zenith', fontsize=18)
        ax.set_xlabel(r'$\cos(\mathrm{\theta_{reco}})$', fontsize=16)
        ax.set_ylabel(r'$\log_{10}(E_{\mathrm{LLH}}/E_{\mathrm{true}})$', fontsize=16)
        lw = 2
        ms = 7*lw
        pltParams = dict(fmt='.', lw=lw, ms=ms)
        # Energy resolution vs zenith
        zbins = np.linspace(1, np.cos(thetaMax), 21)[::-1]
        x = getMids(zbins)
        medians, sigL, sigR, vars = getMedian(x1, y, zbins)
        ax.errorbar(x, medians, yerr=(sigL,sigR), **pltParams)
        # Spline fit
        fitx = np.linspace(x.min(), x.max(), len(x)*3)
        fit = glam.grideval(tab, [fitx])
        ax.plot(fitx, fit)
        ax.set_xlim(0.8,1)
        if out:
            plt.savefig(out)
        plt.show()

    fit = glam.grideval(tab, [np.cos(z)])
    return r - fit



""" Get the energy resolution """
def eres(s, cut=None, xaxis='energy', nbins=20, minE=4.0, rt='t', 
         thetaMax=55., zcorrect=False, title=True, out=False, ax=None):

    titleDict = {'zenith':'Zenith Angle', 'energy':'Energy',
                 'core':'Core Position'}
    xDict = {'zenith':r'$\cos(\mathrm{\theta_{true}})$',
             'energy':r'$\log_{10}(E/\mathrm{GeV})$',
             'core':'Core Position (m)'}
    rtDict = {'r':'Reconstructed', 't':'True'}

    # Setup plot
    if ax == None:
        fig, ax = plt.subplots()
    if title:
        ax.set_title('Energy Resolution vs %s %s' % \
                (rtDict[rt], titleDict[xaxis]), fontsize=18)
    ax.set_xlabel(xDict[xaxis], fontsize=16)
    ax.set_ylabel(r'$\log_{10}(E_{\mathrm{LLH}}/E_{\mathrm{true}})$',
            fontsize=16)
    lw = 2
    ms = 7*lw
    pltParams = dict(fmt='.', lw=lw, ms=ms)

    if xaxis == 'energy':
        Ebins = getEbins()
        bins = np.linspace(Ebins.min(), Ebins.max(), nbins+1)
    if xaxis == 'zenith':
        bins = np.linspace(1, np.cos(thetaMax * degree), nbins+1)[::-1]
    if xaxis == 'core':
        bins = np.linspace(0, 1000, nbins+1)

    t = np.log10(s['MC_energy'])
    r = np.log10(s['ML_energy'])
    if zcorrect:
        r = zfix(s)

    if rt == 't':
        e = t
        z = s['MC_zenith'] 
        cx, cy = s['MC_x'], s['MC_y']
    if rt == 'r':
        e = r
        z = np.pi - s['zenith'] 
        cx, cy = s['ML_x'], s['ML_y']

    # Calculate cut values
    c0 = np.logical_not(np.isnan(r))
    #ecut = (t >= minE) if rt=='t' else (r >= minE)
    ecut = (r >= minE)
    c0 *= ecut
    if cut != None:
        c0 *= s['cuts'][cut]

    # Store median and standard deviation info
    x = getMids(bins)
    x1 = {'energy':e, 'zenith':np.cos(z), 'core':np.sqrt(cx**2+cy**2)}
    x1 = x1[xaxis]
    y = r - t
    medians, sigL, sigR, vars = getMedian(x1[c0], y[c0], bins)
    #if xaxis == 'zenith':
    #    x = (np.arccos(x) / degree)[::-1]
    #    medians, sigL, sigR = medians[::-1], sigL[::-1], sigR[::-1]

    ax.errorbar(x, medians, yerr=(sigL,sigR), **pltParams)

    if out:
        plt.savefig(out)
    #if ax == None:
    plt.show()


""" Get the core resolution """
def core_res(s, cut=None, nbins=40, minE=4, title=True, zcorrect=False, 
        out=False):

    # Setup plot
    fig, ax = plt.subplots()
    if title:
        ax.set_title('Core Resolution vs True Energy', fontsize=18)
    ax.set_xlabel(r'$\log_{10}(E/\mathrm{GeV})$', fontsize=16)
    ax.set_ylabel(r'$\vec{x}_{\mathrm{LLH}} - \vec{x}_{\mathrm{true}} [m]$', 
            fontsize=16)
    lw = 2
    ms = 7*lw
    pltParams = dict(fmt='.', lw=lw, ms=ms)

    # Group into larger bins in energy
    Ebins = getEbins()
    ebins = np.linspace(Ebins.min(), Ebins.max(), nbins+1)
    if ebins[-1] != Ebins[-1]:
        ebins.append(Ebins[-1])

    # Plot Ereco histograms for each energy bin in Etrue
    r = np.log10(s['ML_energy'])
    if zcorrect:
        r = zfix(s)
    c0 = np.logical_not(np.isnan(s['ML_energy']))
    ecut = (r >= minE)
    c0 *= ecut
    if cut != None:
        c0 *= s['cuts'][cut]

    t = np.log10(s['MC_energy'])[c0]
    tx, ty = s['MC_x'][c0], s['MC_y'][c0]
    rx, ry = s['ML_x'][c0], s['ML_y'][c0]

    # Store median and standard deviation info
    x = getMids(ebins)
    y = np.sqrt((rx-tx)**2 + (ry-ty)**2)
    medians, sigL, sigR, vars = getMedian(t, y, ebins)
    ax.errorbar(x, medians, yerr=(sigL,sigR), **pltParams)

    if out:
        plt.savefig(out)
    plt.show()


""" Plot the counts """
def counts(s, cut, name=""):

    # Plot the pure number of counts before unfolding
    r = np.log10(s['ML_energy'])
    t = np.log10(s['MC_energy'])
    Emids = getEmids()
    compDict = {'p':'P', 'h':'He', 'o':'O', 'f':'Fe'}
    colorDict = {'p':'b', 'h':'y', 'o':'c', 'f':'r', 'All':'k'}

    fig, ax = plt.subplots()
    ax.set_title('Counts using '+name+' cut')
    ax.set_xlabel('Log10(Energy/GeV)')
    ax.set_ylabel('Counts')

    # Get the counts
    N = {}
    N['All'] = Nfinder(r, cut)
    for e in compDict.keys():
        ecut = s['llh_comp'] == e
        N[e] = Nfinder(r, cut*ecut)

    # Plot reconstructions
    for e in compDict.keys() + ['All']:
        pnt = colorDict[e] + '.'
        ax.errorbar(Emids, N[e], yerr=np.sqrt(N[e]), fmt=pnt, label=e)

    # MC True values
    for e in compDict.keys():
        truecut = s['comp'] == compDict[e]
        pnt = colorDict[e] + 'x'
        ax.plot(Emids, Nfinder(t, cut*truecut), pnt)
    ax.plot(Emids, Nfinder(t, cut), 'kx')

    ax.set_yscale('log')
    ax.set_ylim([10**(0),10**(6)])
    ax.legend(loc='upper right')

    plt.show()


""" 2-D histogram of ML energy vs MC energy """
def hist2d(s, cut, log=True, out=False, batch=False, zcorrect=False):

    # Setup 2D histogram of Ereco (y) vs Etrue (x) for ML reconstruction
    fig, ax = plt.subplots()
    #ax.set_title('Etrue vs Ereco (actual composition)')
    axParams = {'fontsize':16}
    ax.set_xlabel(r'$\log_{10}(E_{\mathrm{true}}/\mathrm{GeV})$', **axParams)
    ax.set_ylabel(r'$\log_{10}(E_{\mathrm{LLH}}/\mathrm{GeV})$', **axParams)

    ## x translates to the y-axis and vice versa ##
    x = np.log10(s['ML_energy'][cut])
    if zcorrect:
        x = zfix(s)[cut]
    y = np.log10(s['MC_energy'][cut])
    rEbins   = getEbins(reco=True)
    tEbins   = getEbins()

    h, xedges, yedges = np.histogram2d(x, y, bins=(rEbins, tEbins), 
            normed=False, weights=None)
    ntot = np.sum(h, axis=0).astype('float')
    ntot[ntot==0] = 1.
    h /= ntot
    print h.max(), h.min()
    extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
    if log:
        h = np.log10(h)

    im = ax.imshow(h, extent=extent, origin='lower', interpolation='none')
    fig.colorbar(im)
    if out:
        plt.savefig(out)
    if not batch:
        plt.show()


""" Plot a given function for multiple energies """
def e_plot(s, func, cut, nplt, erange=False, out=False):

    r = np.log10(s['MC_energy'])
    Ebins = getEbins()
    n = len(Ebins)/nplt

    for i in range(nplt):
        if erange:
            e_start, e_end = erange
        else:
            e_start, e_end = Ebins[i*n], Ebins[(i+1)*n]
        temp_cut = (r>=e_start) * (r<e_end)
        title = '%s to %s' % (e_start, e_end)
        if out:
            temp_out = '%s_%s_%s.png' % (out, e_start, e_end)
        func(s, cut * temp_cut, title, out=temp_out)



