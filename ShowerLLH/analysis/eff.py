#!/usr/bin/env python

###############################################################################
# Functions associated with finding and analyzing the effective area
###############################################################################

from numpy import *
from useful import *
from scipy import optimize
import matplotlib.pyplot as plt


""" Get the effective area for a given composition and cut """
def getEff(s, cut, comp='joint', zbin=False):

    # Suppress divide by 0 warnings
    inv_err = geterr()['invalid']
    div_err = geterr()['divide']
    seterr(invalid='ignore')
    seterr(divide='ignore')

    eff, sig, relerr = {},{},{}
    a = log10(s['MC_energy'])
    Emids = getEmids()

    c0 = cut
    if comp != 'joint':
        c0 = cut * s[comp]
    if zbin:
        c0 *= s['cuts']['z'+zbin]

    # Set radii for finding effective area
    rDict = {}
    keys = ['low', 'med', 'high']
    for key in keys:
        rDict[key] = array([600, 800, 1100, 1700, 2600, 2900])
    rDict['low'][1] = 600
    Ebreaks = array([4, 5, 6, 7, 8, 9])
    rgrp = digitize(Emids, Ebreaks) - 1

    for key in keys:

        # Get efficiency and sigma
        k = Nfinder(a, c0*s[key+'E'])
        n = s['MC'][comp][key]
        if zbin:
            n = s['MC'][comp][key+'_'+zbin]
        eff[key], sig[key], relerr[key] = zeros((3, len(k)))
        eff[key] = k / n
        var = (k+1)*(k+2)/((n+2)*(n+3)) - (k+1)**2/((n+2)**2)
        sig[key] = sqrt(var)

        # Multiply by throw area
        r = array([rDict[key][i] for i in rgrp])
        eff[key] *= pi*(r**2)
        sig[key] *= pi*(r**2)

        # Deal with parts of the arrays with no information
        for i in range(len(eff[key])):
            if n[i] == 0:
                eff[key][i] = 0
                sig[key][i] = inf

    # Combine low, med, and high energy datasets
    eff_tot = (sum([eff[key]/sig[key] for key in keys], axis=0) /
            sum([1/sig[key] for key in keys], axis=0))
    sig_tot = sqrt(1 / sum([1/sig[key]**2 for key in keys], axis=0))
    relerr  = sig_tot / eff_tot

    seterr(invalid=inv_err)
    seterr(divide=div_err)

    return eff_tot, sig_tot, relerr


""" Flat line fit for effective area at high energies  """
def line_fit(s, cut, comp='joint', st=45):

    lineFit = lambda x, b: b
    Emids = getEmids()
    eff, sigma, relerr = getEff(s, cut, comp=comp)
    Emids, eff, sigma, relerr = Emids[st:], eff[st:], sigma[st:], relerr[st:]
    x = Emids.astype('float64')
    y = eff
    popt, pcov = optimize.curve_fit(lineFit, x, y, sigma=sigma)
    yfit = lineFit(x, *popt)
    yfit = array([yfit for i in range(len(x))])
    chi2 = sum(((yfit - y)/sigma) **2) / len(x)
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
    chi = sum(((yfit - y)/sigma) ** 2) / len(x)
    print "x0, k, a =", popt
    print "chi2 = %e" % chi
    # Fit ignoring error bars
    popt2, pcov2 = optimize.curve_fit(sigmoid, x, y, p0)
    yfit2 = sigmoid(x, *popt2)
    chi2 = sum(((yfit2 - y)/sigma) ** 2) / len(x)
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
        chi2 = sum(((yfit - y)/sigma) ** 2) / len(x)
        chi2List.append(chi2)

    plt.title('Chi2 vs Charge')
    plt.xlabel('Maxcharge > (VEM)')
    plt.ylabel('Chi2')
    plt.plot(QTable, chi2List, '.')
    plt.show()


""" Show the effective area """
def plotter(s, cut, comp='joint', ndiv=False, st=45, zbin=False):

    plt.title('Effective Area vs Energy')
    plt.xlabel('Log10(Energy/GeV)')
    plt.ylabel('Effective Area (m**2)')
    Emids = getEmids()
    Emids = Emids[st:]

    eff, sigma, relerr = getEff(s, cut, comp=comp, zbin=zbin)
    eff, sigma, relerr = eff[st:], sigma[st:], relerr[st:]
    eff2, relerr2 = line_fit(s, cut, comp=comp, st=st)
    # Give the option for combined bins
    if ndiv:
        eff_joint, en_joint = [], []
        for j in range(len(eff)/ndiv):
            start = ndiv*j
            end = ndiv*(j+1)
            eff_joint.append(mean(eff[start:end]))
            en_joint.append(mean(Emids[start:end]))
        plt.plot(en_joint, eff_joint, 'o', label=comp)
    else:
        plt.plot(Emids, eff2)
        plt.errorbar(Emids, eff, yerr=sigma, fmt='.', label=comp)

    plt.legend(loc='upper left')
    #plt.savefig('collab/pics/effarea.png')
    plt.show()
















