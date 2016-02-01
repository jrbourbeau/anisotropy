#!/usr/bin/env python

################################################################################
## Look at the cut efficiency (requires load_sim.py)                          ##
################################################################################

import numpy as np
from scipy import optimize
from useful import *
import matplotlib.pyplot as plt
import time


##===========================================================================##
## Plotting functions

""" Return the median, upper, and lower 34% containment for an array """
def getMed(x):
    if len(x)==0:
        return 0, 0, 0
    x.sort()
    med_index  = len(x)/2
    sigma_step = len(x)*34/100
    med = x[med_index]
    sigma_min = med - x[med_index - sigma_step]
    sigma_max = x[med_index + sigma_step] - med
    return med, sigma_min, sigma_max


""" Get the energy resolution """
def eres(s, r, cut, st=40, ndiv=5, ideal=False, out=False):

    # Plot Ereco histograms for each energy bin in Etrue
    t = np.log10(s['MC_energy'])
    #r = np.log10(s['ML_energy'])
    Ebins, Emids = getEbins(), getEmids()
    if ideal:
        r = (np.log10(s['pML_energy']*s['P'] + s['hML_energy']*s['He'] + 
                s['oML_energy']*s['O'] + s['fML_energy']*s['Fe']))

    # Group into larger bins in energy
    nplt = len(Emids[st:])/ndiv
    x, med, min, max = np.zeros((4,nplt))

    for j in range(nplt):
        start = j*ndiv + st
        end = (j+1)*ndiv + st
        e_cut = (t >= Ebins[start]) * (t < Ebins[end])
        c0 = e_cut * cut
        y = (r - t)[c0]

        # Store median and standard deviation info
        x[j] = (Ebins[end] + Ebins[start]) / 2.
        med[j], min[j], max[j] = getMed(y)

    # Setup plot
    f = plt.figure()
    ax1 = f.add_subplot(1,1,1)
    ax1.set_title('Energy Resolution vs True Energy', fontsize=18)
    ax1.set_xlabel('Log10(Energy/GeV)', fontsize=16)
    ax1.set_ylabel('Ereco - Etrue (median)', fontsize=16)

    # Plot parameters
    lw = 2
    ms = 7*lw

    # Plot energy resolution
    plt.errorbar(x, med, yerr=(min,max), fmt='.', lw=lw, ms=ms)
    if out:
        plt.savefig(out)
    plt.show()


""" Get the energy resolution as a function of zenith angle """
def eres_z(s, cut, minE=6.0, zbins=5, thetaMax=np.arccos(0.8)*180/np.pi, out=False):

    # Plot Ereco histograms for each energy bin in Etrue
    t = np.log10(s['MC_energy'])
    r = np.log10(s['ML_energy'])
    rz = s['ShowerPlane_zenith']
    tz = s['MC_zenith']
    #cut = s['cuts']['llh']

    # Apply minimum energy cut
    #e_cut = (r >= minE)
    e_cut = (t >= minE)
    c0 = cut * e_cut

    # Group into equal area zenith bins
    theta0 = 0
    thetaMax *= np.pi / 180.
    Zbins = np.array([theta0])
    sigma = 2*pi*(np.cos(theta0) - np.cos(thetaMax)) / zbins
    while len(Zbins) < zbins+1:
        theta  = np.arccos(np.cos(theta0) - sigma/(2*pi))
        Zbins  = append(Zbins, theta)
        theta0 = theta

    x, med, min, max = np.zeros((4,zbins))
    for j in range(zbins):
        z_cut = (rz >= Zbins[j]) * (rz < Zbins[j+1])
        #z_cut = (tz >= Zbins[j]) * (tz < Zbins[j+1])
        c1 = z_cut * c0
        y = (r - t)[c1]

        # Store median and standard deviation info
        x[j] = (Zbins[j] + Zbins[j+1]) / 2.
        x[j] *= 180/pi
        med[j], min[j], max[j] = getMed(y)

    # Setup plot
    f = plt.figure()
    ax1 = f.add_subplot(1,1,1)
    ax1.set_title('Energy Resolution vs Zenith Angle', fontsize=18)
    ax1.set_xlabel('Zenith (deg)', fontsize=16)
    ax1.set_ylabel('Ereco - Etrue (median)', fontsize=16)
    #ax1.set_xlim(0,40)
    #ax1.set_ylim(-.2,.2)

    # Plot parameters
    lw = 2
    ms = 7*lw

    # Plot energy resolution
    plt.errorbar(x, med, yerr=(min,max), fmt='.', lw=lw, ms=ms)
    if out:
        plt.savefig(out)
    plt.show()


""" Get the core resolution """
def core_res(s, r, cut, st=40, ndiv=5, out=False):

    # Plot Ereco histograms for each energy bin in Etrue
    t = np.log10(s['MC_energy'][cut])
    tx, ty = s['MC_x'][cut], s['MC_y'][cut]
    #rx, ry = s['ML_x'][cut], s['ML_y'][cut]
    rx, ry = s[r+'_x'][cut], s[r+'_y'][cut]
    Ebins, Emids = getEbins(), getEmids()

    # Group into larger bins in energy
    nplt = len(Emids[st:])/ndiv
    x, med, min, max = np.zeros((4,nplt))

    for j in range(nplt):
        start = j*ndiv + st
        end = (j+1)*ndiv + st
        e_cut = (t >= Ebins[start]) * (t < Ebins[end])
        y = sqrt((rx-tx)**2 + (ry-ty)**2)[e_cut]

        # Store median and standard deviation info
        x[j] = (Ebins[end] + Ebins[start]) / 2.
        med[j], min[j], max[j] = getMed(y)

    # Setup plot
    f = plt.figure()
    ax1 = f.add_subplot(1,1,1)
    ax1.set_title('Core Resolution vs True Energy', fontsize=18)
    ax1.set_xlabel('Log10(Energy/GeV)', fontsize=16)
    ax1.set_ylabel('Reco_core - True_core (median)', fontsize=16)

    # Plot parameters
    lw = 2
    ms = 7*lw

    # Plot energy resolution
    plt.errorbar(x, med, yerr=(min,max), fmt='.', lw=lw, ms=ms)
    if out:
        plt.savefig(out)
    plt.show()


""" Look at the zenith distribution for a variety of energies  """
def z_check(s, r, cut):

    f = plt.figure()
    ax1 = f.add_subplot(1,1,1)
    ax1.set_title('Zenith Distribution for Various Energy Bins')
    ax1.set_xlabel('Reconstructed Zenith (Degrees)')
    ax1.set_ylabel('Counts')

    t = np.log10(s['MC_energy'])
    z = s['ShowerPlane_zenith'] * 180./pi

    # Setup plot

    zbins = arange(0,41,0.5)
    ebins = arange(5, 7.6, 0.5)

    for i in range(len(ebins)-1):
        reco_ecut = (r>=ebins[i]) * (r<ebins[i+1])
        true_ecut = (t>=ebins[i]) * (t<ebins[i+1])
        c0 = reco_ecut * cut
        c1 = true_ecut * cut
        label = '%s-%s' % (ebins[i], ebins[i+1])
        reco_label = 'reco ' + label
        true_label = 'true ' + label

        plt.hist(z[c0], bins=zbins, histtype='step', label=reco_label)
        plt.hist(z[c1], bins=zbins, histtype='step', label=true_label)

    plt.legend()
    plt.show()



""" Plot the counts """
def counts(s, cut, name=""):

    # Plot the pure number of counts before unfolding
    r = np.log10(s['ML_energy'])
    t = np.log10(s['MC_energy'])
    Emids = getEmids()
    #compDict = {'p':'P','h':'He','o':'O','f':'Fe'}
    compDict = {'p':'P', 'f':'Fe'}
    colorDict = {'p':'b', 'h':'y', 'o':'c', 'f':'r', 'All':'k'}

    plt.title('Counts using '+name+' cut')
    plt.xlabel('Log10(Energy/GeV)')
    plt.ylabel('Counts')

    # Get the counts
    N = {}
    N['All'] = Nfinder(r, cut)
    for e in compDict.keys():
        N[e] = Nfinder(r, cut*s[e])

    # Plot reconstructions
    for e in compDict.keys() + ['All']:
        pnt = colorDict[e] + '.'
        plt.errorbar(Emids, N[e], yerr=sqrt(N[e]), fmt=pnt, label=e)

    # MC True values
    for e in compDict.keys():
        true = compDict[e]
        pnt = colorDict[e] + 'x'
        plt.plot(Emids, Nfinder(t, cut*s[true]), pnt)
    plt.plot(Emids, Nfinder(t, cut), 'kx')

    plt.yscale('log')
    plt.legend(loc='upper right')
    plt.ylim([10**(0),10**(6)])

    plt.show()


""" 2-D histogram of ML energy vs MC energy """
def hist2d(s, cut, log=True):

    # Setup 2D histogram of Ereco (y) vs Etrue (x) for ML reconstruction
    f = plt.figure()
    ax1 = f.add_subplot(1,1,1)
    ax1.set_title('Etrue vs Ereco (actual composition)')
    ax1.set_xlabel('log10(MC/GeV)')
    ax1.set_ylabel('log10(LLH/GeV)')

    ## x translates to the y-axis and vice versa ##
    x = np.log10(s['ML_energy'][cut])
    y = np.log10(s['MC_energy'][cut])
    Ebins   = getEbins()
    myBins  = len(Ebins)-1
    myRange = (4, 9.5)

    h, xedges, yedges = histogram2d(x, y, bins=myBins, range=(myRange,myRange), normed=False, weights=None)
    extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
    if log:
        h = np.log10(h)

    im = ax1.imshow(h, extent=extent, origin='lower')
    f.colorbar(im)
    plt.show()


""" Look for separation in true composition using the llh ratio """
def llhratio(s, cut, title, out=False):

    f = plt.figure()
    ax1 = f.add_subplot(1,1,1)
    ax1.set_title('Likelihood Ratio from ' + title, fontsize=18)
    ax1.set_xlabel('Log-Likelihood Difference (pLLH-fLLH)', fontsize=16)
    ax1.set_ylabel('Counts', fontsize=16)

    r = np.log10(s['ML_energy'])
    pf_ratio = s['pLLH'] - s['fLLH']
    colorDict = {'P':'b','Fe':'r'}

    for true_comp in ['P','Fe']:

        c0 = s[true_comp] * cut
        if c0.sum()==0:
            continue
        bins = linspace(-5,5,100)
        plt.hist(pf_ratio[c0], bins=bins, histtype='step', label=true_comp, linewidth=2, color=colorDict[true_comp])

    plt.legend()
    if out:
        plt.savefig(out)
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


""" Look at a given parameter split for both iron and proton """
def pfplot(s, varName, bins=10, xlog=False, ylog=False):

    f = plt.figure()
    ax1 = f.add_subplot(1,1,1)
    ax1.set_title(varName)
    ax1.set_xlabel(varName)
    ax1.set_ylabel('counts')

    cut = s['cuts']['llh']
    nancut = (s[varName]==s[varName])
    cut *= nancut

    for e in ['P', 'H', 'O', 'F']:
        c0 = cut * s[e]
        x = s[varName][c0]
        h, myBins, temp = ax1.hist(x, bins=bins, histtype='step', label=e)

    ax1.legend(loc='lower left')
    if xlog:
        ax1.set_xscale('log')
    if ylog:
        ax1.set_yscale('log')

    plt.show()





