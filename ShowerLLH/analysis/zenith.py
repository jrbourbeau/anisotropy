#!/usr/bin/env python

##############################################################################
## Analyze the effective area and flux as a function of zenith angle
## for a variety of compositions
##############################################################################

from numpy import *
import sim, powerFit
from eff import getEff
from data import NumToFlux
import matplotlib.pyplot as plt
from useful import *


# Define zenith cuts
def get_zcuts(q, nbins):
    zcuts = {}
    cosz = cos(q['zenith'])
    if nbins == 4:
        zcuts['40'] = (cosz >= 0.95)
        zcuts['41'] = (cosz >= 0.90) * (cosz < 0.95)
        zcuts['42'] = (cosz >= 0.85) * (cosz < 0.90)
        zcuts['43'] = (cosz >= 0.80) * (cosz < 0.85)
    elif nbins == 3:
        zcuts['30'] = (cosz >= 0.93)
        zcuts['31'] = (cosz >= 0.87) * (cosz < 0.93)
        zcuts['32'] = (cosz >= 0.80) * (cosz < 0.87)
    return zcuts


def effCheck(s, cut, st=45, comp='joint', ndiv=1, zbins=3):

    # Define zenith cuts
    zcuts = get_zcuts(s, zbins)

    # Setup plot
    f = plt.figure()
    ax1 = f.add_subplot(1,1,1)
    ax1.set_title('Effective Area for ' + comp)
    ax1.set_xlabel('Log10(Energy/GeV)')
    ax1.set_ylabel('Effective Area (m^2)')
    Emids = getEmids()
    Emids = Emids[st:]

    keys = zcuts.keys()
    keys.sort()
    for key in keys:
        tcut = cut * zcuts[key]
        eff, sigma, relerr = getEff(s, tcut, comp=comp, smooth=False, zbin=key)
        fit, sigma, relerr = getEff(s, tcut, comp=comp, smooth=True, zbin=key)
        # Give the option for combined bins
        nplt = len(eff)/ndiv
        eff_joint, rel_joint, err_joint, en_joint = zeros((4, nplt))
        for j in range(nplt):
            start = ndiv*j
            end = ndiv*(j+1)
            eff_joint[j] = mean(eff[start:end])
            rel_joint[j] = mean(relerr[start:end])
            en_joint[j]  = mean(Emids[start:end])

        err_joint = eff_joint * rel_joint
        ax1.errorbar(en_joint, eff_joint, yerr=err_joint, fmt='.', label=key)
        ax1.plot(Emids, fit)

    ax1.legend(loc='lower left')
    plt.show()


def countCheck(s, st=45, comp='joint', flux=False):

    Emids = getEmids()
    Emids = Emids[st:]
    cut = s['cuts']['llh']

    ekey = 'ML_energy'
    if comp != 'joint':
        ekey = comp.lower() + ekey
    r = log10(s[ekey])

    # Setup plot
    f = plt.figure()
    ax1 = f.add_subplot(1,1,1)
    ax1.set_title('Counts for ' + comp)
    ax1.set_xlabel('Log10(Energy/GeV)')
    ax1.set_ylabel('Counts')

    # Toy MC spectrum
    sDict = {}
    for key in ['P','H','O','F']:               # Mixed spectrum
        sDict[key] = [-0.25, [[7.5, -0.75]]]
    sDict['F'] = [2, [[7.5, 0]]]
    temp = {}
    for key in sDict.keys():
        s0 = sDict[key][0]
        sTable = sDict[key][1]
        temp[key] = powerFit.powerFit(s0, sTable=sTable, st=st)
    specCut = sim.fakeSpec(s, cut, temp)

    zcuts = get_zcuts(s, 3)
    keys = zcuts.keys()
    keys.sort()
    for key in keys:
        #counts = Nfinder(r, cut*zcuts[key])[st:]
        counts = Nfinder(r, cut*zcuts[key]*specCut)[st:]
        fit, sigma, relerr = getEff(s, cut*zcuts[key], comp=comp, smooth=True, zbin=key)
        rel = sqrt(1/counts + relerr**2)
        if flux:
            flux = counts / fit
            ax1.errorbar(Emids, flux, yerr=rel*flux, fmt='.', label=key)
        else:
            ax1.errorbar(Emids, counts, yerr=rel*counts, fmt='.', label=key)

    ax1.set_yscale('log')
    ax1.legend(loc='lower left')
    plt.show()
    

def fluxCheck(s, d, st=45, comp='joint', spec=2.7):

    monthList = ['201006']
    t = d['duration']
    wList = [d['w1'], d['w24'], d['w8']]
    Emids = getEmids()
    Emids = Emids[st:]
    scut = s['cuts']['llh']
    dcut = d['cuts']['llh']
    scale = (10**Emids)**spec

    ekey = 'ML_energy'
    if comp != 'joint':
        ekey = comp.lower() + ekey
    r = log10(d[ekey])

    # Setup plot
    f = plt.figure()
    ax1 = f.add_subplot(1,1,1)
    ax1.set_title('Flux for ' + comp)
    ax1.set_xlabel('Log10(Energy/GeV)')
    ax1.set_ylabel('Flux')

    dzcuts = get_zcuts(d, 3)
    szcuts = get_zcuts(s, 3)
    keys = dzcuts.keys()
    keys.sort()
    for key in keys:
        counts = Nfinder(r, dcut*dzcuts[key], wList)[st:]
        fit, sigma, relerr = getEff(s, scut*szcuts[key], comp=comp, smooth=True, zbin=key)
        rel = sqrt(1/counts + relerr**2)
        flux = NumToFlux(counts, fit, t, st) * scale
        ax1.errorbar(Emids, flux, yerr=rel*flux, fmt='.', label=key)
        #ax1.errorbar(Emids, counts, yerr=rel*counts, fmt='.', label=key)

    ax1.set_yscale('log')
    ax1.legend(loc='lower left')
    plt.show()


def compCheck(s, d, cut, st=45, zbins=4):

    r = log10(s['ML_energy'])
    Emids = getEmids()
    Emids = Emids[st:]

    zcuts = get_zcuts(s, zbins)
    pLike = (s['maxLLH_p'] > s['maxLLH_f'])
    fLike = logical_not(pLike)

    f = plt.figure()
    ax1 = f.add_subplot(1,1,1)
    ax1.set_title('# of protons / # of iron')
    ax1.set_xlabel('Log10(Energy/GeV)')
    ax1.set_ylabel('ratio (p/f)')

    keys = zcuts.keys()
    keys.sort()
    for key in keys:
        np = Nfinder(r, cut*pLike*zcuts[key])[st:]
        nf = Nfinder(r, cut*fLike*zcuts[key])[st:]
        ax1.plot(Emids, np/nf, label=key)

    ax1.legend(loc='lower left')
    plt.show()


def lnA(q, cut, st=45, zbins=4, mc=False):

    Emids = getEmids()
    Emids = Emids[st:]
    zcuts = get_zcuts(q, zbins)

    N = {}
    if not mc:
        rList = ['p','h','o','f']
        r = log10(q['ML_energy'])
    else:
        rList = ['P','H','O','F']
        r = log10(q['MC_energy'])
    wts = {'p':1, 'h':4, 'o':16, 'f':56}
    Ntot = Nfinder(r, cut)[st:]
    for e in rList:
        N[e] = Nfinder(r, cut*q[e])[st:]
        N[e] *= (wts[e.lower()] / Ntot)

    a = sum([N[e] for e in rList], axis=0)
    
    f = plt.figure()
    ax1 = f.add_subplot(1,1,1)
    ax1.set_title('Ln(A)')
    ax1.set_xlabel('Log10(Energy/GeV)')
    ax1.set_ylabel('ln(A)')

    ax1.plot(Emids, log(a), '.')

    ax1.set_ylim(2.6, 3.5)
    plt.show()









