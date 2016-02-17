#!/usr/bin/env python

###############################################################################
# Functions associated with finding and analyzing the effective area
###############################################################################

import argparse
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt

from load_sim import load_sim
from simFunctions_IT import getErange
from llhtools import *

import myGlobals as my
from useful import getMids


""" Get the effective area for a given composition and cut """
def getEff(s, cut, comp='joint', reco=True):

    eff, sig, relerr = {},{},{}
    a = np.log10(s['MC_energy'])
    Ebins = getEbins()
    Emids = getMids(Ebins)
    erangeDict = getErange()

    c0 = cut
    if comp != 'joint':
        compcut = s['comp'] == comp
        c0 = cut * compcut

    # Set radii for finding effective area
    rDict = {}
    keys = ['low', 'mid', 'high']
    for key in keys:
        rDict[key] = np.array([600, 800, 1100, 1700, 2600, 2900])
    rDict['low'][1] = 600
    Ebreaks = np.array([4, 5, 6, 7, 8, 9])
    rgrp = np.digitize(Emids, Ebreaks) - 1

    for key in keys:

        # Get efficiency and sigma
        simcut = np.array([sim in erangeDict[key] for sim in s['sim']])
        k = np.histogram(a[c0*simcut], bins=Ebins)[0]
        #k = Nfinder(a, c0*simcut)
        n = s['MC'][comp][key].astype('float')
        eff[key], sig[key], relerr[key] = np.zeros((3, len(k)))
        with np.errstate(divide='ignore', invalid='ignore'):
            eff[key] = k / n
            var = (k+1)*(k+2)/((n+2)*(n+3)) - (k+1)**2/((n+2)**2)
        sig[key] = np.sqrt(var)

        # Multiply by throw area
        r = np.array([rDict[key][i] for i in rgrp])
        eff[key] *= np.pi*(r**2)
        sig[key] *= np.pi*(r**2)

        # Deal with parts of the arrays with no information
        for i in range(len(eff[key])):
            if n[i] == 0:
                eff[key][i] = 0
                sig[key][i] = np.inf

    # Combine low, mid, and high energy datasets
    eff_tot = (np.sum([eff[key]/sig[key] for key in keys], axis=0) /
            np.sum([1/sig[key] for key in keys], axis=0))
    sig_tot = np.sqrt(1 / np.sum([1/sig[key]**2 for key in keys], axis=0))
    with np.errstate(divide='ignore'):
        relerr  = sig_tot / eff_tot

    # UGH : find better way to do this
    if reco:
        eff_tot = eff_tot[20:]
        sig_tot = sig_tot[20:]
        relerr  = relerr[20:]

    return eff_tot, sig_tot, relerr


""" Get effective area information fast for preset cuts """
def getEff_fast(linear=False):
    my.setupShowerLLH(verbose=False)
    eFile = '%s/EffectiveArea.npy' % my.llh_resource
    temp = np.load(eFile)
    temp = temp.item()
    effarea, sigma, relerr = temp['effarea'], temp['sigma'], temp['relerr']
    if linear:
        effarea = np.ones(len(effarea), dtype=float)
    return effarea, sigma, relerr


""" Flat line fit for effective area at high energies  """
def line_fit(s, cut, comp='joint', st=45):

    lineFit = lambda x, b: b
    Emids = getMids(getEbins())
    eff, sigma, relerr = getEff(s, cut, comp=comp)
    Emids, eff, sigma, relerr = Emids[st:], eff[st:], sigma[st:], relerr[st:]
    x = Emids.astype('float64')
    y = eff
    popt, pcov = optimize.curve_fit(lineFit, x, y, sigma=sigma)
    yfit = lineFit(x, *popt)
    yfit = np.array([yfit for i in range(len(x))])
    chi2 = np.sum(((yfit - y)/sigma) **2) / len(x)
    print 'Chi2:', chi2
    return yfit, relerr



""" Define sigmoid function """
def sigmoid(x, x0, k, a):
    y = a / (1 + exp(-k*(x-x0)))
    return y


""" Fit a sigmoid function to the effective area """
def sig_fit(s, cut, comp='joint'):

    Emids = getEmids()
    eff, sigma, relerr = getEff(s, cut, comp=comp)
    x = Emids.astype('float64')
    y = eff

    p0 = [5.75, 6.75, 493000]
    yguess = sigmoid(x, *p0)

    # Fit with error bars
    popt, pcov = optimize.curve_fit(sigmoid, x, y, p0, sigma)
    yfit = sigmoid(x, *popt)
    chi = np.sum(((yfit - y)/sigma) ** 2) / len(x)
    print "x0, k, a =", popt
    print "chi2 = %e" % chi
    # Fit ignoring error bars
    popt2, pcov2 = optimize.curve_fit(sigmoid, x, y, p0)
    yfit2 = sigmoid(x, *popt2)
    chi2 = np.sum(((yfit2 - y)/sigma) ** 2) / len(x)
    print "x02, k2, a2 =", popt2
    print "chi2 = %e" % chi2

    plt.title('Fitting Effective Area ('+comp+')')
    plt.xlabel('Log10(Energy/GeV)')
    plt.ylabel('Effective Area (m**2)')

    plt.errorbar(x, y, yerr=sigma, fmt='.', label='data')
    plt.plot(x, yfit, 'k', label='fit')
    plt.plot(x, yfit2, 'g', label='fit2')
    plt.plot(x, yguess, 'r', label='guess')
    plt.legend(loc='upper left')
    plt.grid(True)
    #plt.yscale('log')

    plt.show()


""" Look at chi2 for flat line fit for a variety of cuts """
def chi2(s, cut):

    lineFit = lambda x, b: b
    st = 44
    Emids = getEmids()

    chi2List = []
    QTable = arange(1, 10, 0.25)
    for i in range(len(QTable)):
        qcut = cut * (s['Q1']>QTable[i])
        eff, sigma, relerr = getEff(s, qcut)

        x = Emids.astype('float64')[st:]
        y = eff[st:]

        popt, pcov = optimize.curve_fit(lineFit, x, y, sigma=sigma)
        yfit = lineFit(x, *popt)
        yfit = [yfit for i in range(len(x))]
        chi2 = np.sum(((yfit - y)/sigma) ** 2) / len(x)
        chi2List.append(chi2)

    plt.title('Chi2 vs Charge')
    plt.xlabel('Maxcharge > (VEM)')
    plt.ylabel('Chi2')
    plt.plot(QTable, chi2List, '.')
    plt.show()


""" Show the effective area """
def plotter(s, cut, comp='joint', emin=6.2, ndiv=False, out=False, fit=False):

    fig, ax = plt.subplots()
    #ax.set_title('Effective Area vs Energy')
    ax.set_xlabel(r'$\log_{10}(E/\mathrm{GeV})$')
    ax.set_ylabel(r'Effective Area ($m^2$)')
    Emids = getMids(getEbins(reco=True))

    eff, sigma, relerr = getEff(s, cut, comp=comp)
    #eff2, relerr2 = line_fit(s, cut, comp=comp, st=st)

    lineFit = lambda x, b: b
    x = Emids.astype('float64')
    y = eff
    c0 = x >= emin
    popt, pcov = optimize.curve_fit(lineFit, x[c0], y[c0], sigma=sigma[c0])
    yfit = lineFit(x[c0], *popt)
    yfit = np.array([yfit for i in range(len(x[c0]))])
    chi2 = np.sum(((yfit - y[c0])/sigma[c0]) **2) / len(x[c0])
    print 'Chi2:', chi2
    eff2 = yfit

    # Give the option for combined bins
    if ndiv:
        eff_joint, en_joint = [], []
        for j in range(len(eff)/ndiv):
            start = ndiv*j
            end = ndiv*(j+1)
            eff_joint.append(mean(eff[start:end]))
            en_joint.append(mean(Emids[start:end]))
        ax.plot(en_joint, eff_joint, 'o', label=comp)
    else:
        ax.errorbar(Emids, eff, yerr=sigma, fmt='.', label=comp)
        if fit:
            ax.plot(Emids[c0], eff2)

    #ax.legend(loc='upper left')
    #ax.set_yscale('log')
    if out:
        plt.savefig(out)
    plt.show()



if __name__ == "__main__":

    # Setup global paths
    my.setupShowerLLH(verbose=False)

    p = argparse.ArgumentParser()
    p.add_argument('-c', '--config', dest='config',
            default='IT73',
            help='Detector configuration [IT73|IT81]')
    args = p.parse_args()

    s = load_sim(config=args.config, bintype='logdist')
    d = {}
    d['effarea'], d['sigma'], d['relerr'] = getEff(s, s['cuts']['llh'])

    outFile = '%s/%s_EffectiveArea.npy' % (my.llh_resource, args.config)
    np.save(outFile, d)














