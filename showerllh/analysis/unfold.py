#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

from llhtools import *
import sim, eff, powerFit, spline
from prob import getProbs, p2
from fakeSpec import fakeSpec
from useful import getMids

##===========================================================================##
## Useful functions

## Perform the Bayesian unfolding ##
def unfold(p, diag=False, reco=True):

    # Build composition list
    eList = list(set([k[1] for k in p if '|' in k and k[0]=='R']))
    tList = list(set([k[-1] for k in p if '|' in k and k[0]=='R']))
    eTable = [[r, t] for r in eList for t in tList]

    # Bayesian unfolding
    #p2(p, 'Rf|Tp', reco=reco)
    for r, t in eTable:
        rr, tt = 'R'+r, 'T'+t
        p[tt+'|'+rr] = ((p[rr+'|'+tt] * p[tt]).transpose() /
                np.nansum([p[rr+'|T'+k]*p['T'+k] for k in tList], axis=(0,2)))
    #p2(p, 'Tf|Rp', reco=reco)

    # Option to replace the probability tables with diagonal matrices
    # Used for testing purposes
    if diag:
        keyList = [key for key in p.keys() if '|' in key]
        for key in keyList:
            p[key] = np.ones(p[keyList[0]].shape)
            p[key] = np.diag(np.diag(p[key]))

    # Normalize
    for t in tList:
        pTable = [p['T'+t+'|R'+e]*p['R'+e] for e in eList]
        p['T'+t] = np.nansum(pTable, axis=(0,2))
        # Not sure yet if this is the right way to do it
        p['T'+t][np.isnan(p['T'+t])] = 0
    sumprobs = np.sum([p['T'+t] for t in tList])
    for t in tList:
        p['T'+t] /= sumprobs

    return p


##===========================================================================##
## Applying unfolding to simulation

## Plot the counts after unfolding ##
def counts(s, niter, sDict={'':[-0.25, [[7.5, -0.75]]]}, 
                        smooth=True, spl=False, out=False,
                        reco=True, zcorrect=False):

    r = np.log10(s['ML_energy'])
    rEbins = getEbins(reco=True)
    tEbins = getEbins(reco=reco)
    rEmids = getMids(rEbins)
    tEmids = getMids(tEbins)
    cutName = 'llh'
    cut = s['cuts'][cutName]
    eList = getComps(s)
    tList = getComps(s, reco=False)

    fig, ax = plt.subplots()
    ax.set_title('Energy Spectum using '+cutName+' cut')
    ax.set_xlabel('Log10(Energy/GeV)')
    ax.set_ylabel('Counts')

    # Load probability tables
    p = getProbs(s, cut, reco=reco, zcorrect=zcorrect)

    # Option for splining
    if spl:
        for key in p.keys():
            p[key] = 10**(spline.spline(s, p[key], nk=2, npoly=3))
        print sum(p['Rf|Tp']+p['Rp|Tp'], axis=0)
        print sum(p['Rf|Tf']+p['Rp|Tf'], axis=0)

    # Create our toy MC spectrum
    temp = {}
    for key in sDict.keys():
        s0 = sDict[key][0]
        sTable = sDict[key][1]
        temp[key] = powerFit.powerFit(s0, sTable=sTable, reco=reco)
    specCut = fakeSpec(s, cut, temp, reco=reco)
    #specCut = np.array([True for i in range(len(specCut))])

    # Create starting arrays
    N_passed = float(np.histogram(r[cut*specCut], bins=rEbins)[0].sum())
    N = {}
    for e in eList:
        recocut = s['llh_comp'] == e
        N[e] = np.histogram(r[cut*recocut*specCut], bins=rEbins)[0]
        p['R'+e] = N[e] / N_passed
        p['T'+e] = powerFit.powerFit(-2.7, reco=reco)

    # Get relative errors due to unfolding
    #chi2, unrel = load('unfold_err.npy')
    # Due to efficiency
    #effarea, sigma, relerr = eff.getEff(s, cut, smooth=smooth)
    effarea, sigma, relerr = eff.getEff(s, cut, reco=reco)

    # Bayesian unfolding
    for i in range(niter):
        p = unfold(p, reco=reco)
        #All_unfold = (p['Tf']+p['Tp']) * N_passed
        #ax.errorbar(Emids, All_unfold, yerr=unrel[i]*All_unfold, label=str(i))
        # Smooth prior before next iteration (except last time)
        if i < niter-1:
            for e in eList:
                p['T'+e] = smoother(p['T'+e])

    # Find bin values
    Nun, Rel, Err = {},{},{}
    for e in eList:
        Nun[e] = p['T'+e] * N_passed
    Nun['All'] = np.sum([Nun[e] for e in eList], axis=0)

    # Calculate errors
    for e in eList + ['All']:
        ## NOTE: I don't think you can use the relative errors like this ##
        #Rel[e] = np.sqrt(1/Nun[e] + relerr**2 + unrel[niter-1]**2) 
        Rel[e] = np.sqrt(1/Nun[e] + relerr**2) 
        Err[e] = Nun[e] * Rel[e]
        # And plot
        pnt = getColor(e) + '.'
        ax.errorbar(tEmids, Nun[e], yerr=Err[e], fmt=pnt, label=e)

    # Attempt to fit with broken power law
    #p0 = [1, 10**(-4.5), 7.5, -0.5]
    #yfit = powerFit.pow_fit(10**Emids, Nun['f'], Err['f'], p0)
    #ax.plot(Emids, yfit, label='test')

    #err = {}
    #err[0] = 1/np.sqrt(All_Nunfold)
    #err[1] = relerr
    #err[2] = unrel[niter-1]

    # Plot error bars nicely
    #errCalc = lambda i: np.sqrt(np.sum([err[j]**2 for j in range(i+1)], axis=0))
    #uperr, dnerr = {}, {}
    #for i in range(len(err)):
    #    uperr[i] = All_Nunfold * (1 + errCalc(i))
    #    dnerr[i] = All_Nunfold * (1 - errCalc(i))
    #    for j in range(len(dnerr[i])):
    #        if dnerr[i][j] < 0:
    #            dnerr[i][j] = 10**(-3)

    #ax.errorbar(Emids, All_Nunfold, yerr=All_err, fmt='gx', label='Unfold')
    #ax.plot(Emids, All_Nunfold, 'b.', label='Unfold')
    #for i in range(len(err)):
    #    ax.fill_between(Emids, dnerr[i], uperr[i], facecolor='blue', alpha=0.2)

    # Plot true spectrum
    MC = {}
    etrue = np.log10(s['MC_energy'])
    MC['All'] = np.histogram(etrue[cut*specCut], bins=tEbins)[0]
    for t in tList:
        truecut = s['comp'] == t
        MC[t] = np.histogram(etrue[cut*specCut*truecut], bins=tEbins)[0]
    for t in tList + ['All']:
        pnt = getColor(t) + 'x'
        ax.plot(tEmids, MC[t], pnt, label='MC_'+t)

    # plot original (not unfolded) spectrum
    O_N = (np.sum([N[e] for e in eList], axis=0))
    # Just for now...
    temperr = relerr if reco else relerr[20:]
    O_relerr = np.sqrt(1/O_N + temperr**2)
    O_err = O_N * O_relerr
    ax.errorbar(rEmids, O_N, yerr=O_err, fmt='gx', label='Original')

    ax.set_yscale('log')
    #ax.legend(loc='lower left')
    #ax.set_ylim((10**(-1), 10**(4)))
    if out:
        plt.savefig('collab/pics/'+out+'.png')
    plt.show()






