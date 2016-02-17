#!/usr/bin/env python

##=============================================================================
## Contains unfinished functions removed from sim.py
##=============================================================================

import numpy as np
import matplotlib.pyplot as plt

""" Show the splits in counts based on energy cuts """
def esplit(s, st=0, ndiv=5, out=False):

    # Plot Ereco histograms for each energy bin in Etrue
    t = np.log10(s['MC_energy'])
    r = np.log10(s['ML_energy'])
    cut = s['rrc']
    Ebins, Emids = getEbins(), getEmids()

    # Setup plot
    fig, ax = plt.subplots()
    ax.set_title('Counts vs True Energy for Reco Energy Cuts', fontsize=18)
    ax.set_xlabel('Log10(Energy/GeV)', fontsize=16)
    ax.set_ylabel('Counts', fontsize=16)

    # Group into larger bins in energy
    nplt = len(Emids[st:80])/ndiv
    h = np.zeros((nplt, len(Emids)))
    labels = []
    for j in range(nplt):
        start = j*ndiv + st
        end = (j+1)*ndiv + st
        e_cut = (r >= Ebins[start]) * (r < Ebins[end])
        c0 = e_cut * cut
        labels += ['%s to %s' % (Ebins[start], Ebins[end])]
        h[j] += np.histogram(t[c0], bins=Ebins)[0]
        #ntot = float(h[0].sum())
        #plt.plot(Emids, h[0]/ntot, '.', label=label)

    nEtrue = h.sum(axis=0)
    for j in range(nplt):
        ax.plot(Emids, h[j]/nEtrue, label=labels[j])

    ax.set_yscale('log')
    ax.legend(loc='lower right')
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
    z = s['ShowerPlane_zenith'] * 180./np.pi

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
        plt.hist(pf_ratio[c0], bins=bins, histtype='step', 
                label=true_comp, linewidth=2, color=colorDict[true_comp])

    plt.legend()
    if out:
        plt.savefig(out)
    plt.show()


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


""" Get the energy resolution as a function of llh value """
def eres_llh(s, nbins=20, minE=4.0, thetaMax=55, out=False):

    # Plot Ereco histograms for each energy bin in Etrue
    t = np.log10(s['MC_energy'])
    r = np.log10(s['ML_energy'])
    rz = s['ShowerPlane_zenith']
    llh = s['fLLH']

    # Apply cuts
    e_cut = (r >= minE)
    z_cut = (rz <= thetaMax * np.pi/180.)
    c0 = e_cut * z_cut

    bins = np.linspace(llh.min(), llh.max(), nbins+1)
    x, med, min, max = np.zeros((4,nbins))

    for j in range(nbins):
        llh_cut = (llh >= bins[j]) * (llh < bins[j+1])
        c1 = llh_cut * c0
        y = (r - t)[c1]

        # Store median and standard deviation info
        x[j] = (bins[j] + bins[j+1]) / 2.
        med[j], min[j], max[j] = getMed(y)

    # Setup plot
    f = plt.figure(figsize=(17,6))
    lw = 2
    ms = 7*lw

    ax = f.add_subplot(121)
    ax.set_title('Energy Resolution vs ShowerLLH Likelihood', fontsize=18)
    ax.set_xlabel('LLH', fontsize=16)
    ax.set_ylabel('Ereco - Etrue (median)', fontsize=16)
    ax.errorbar(x, med, yerr=(min,max), fmt='.', lw=lw, ms=ms)

    ax = f.add_subplot(122)
    ax.set_title('Energy Resolution vs ShowerLLH Likelihood', fontsize=18)
    ax.set_xlabel('LLH', fontsize=16)
    ax.set_ylabel('Counts', fontsize=16)
    ax.hist(llh, bins=bins, histtype='step', log=True)

    # Plot energy resolution
    if out:
        plt.savefig(out)
    plt.show()

