#!/usr/bin/env python

from numpy import *
import sim
from useful import *
import matplotlib.pyplot as plt
from prob import getProbs
from fakeSpec import fakeSpec
import powerFit, spline

##===========================================================================##
## Useful functions

## Perform the Bayesian unfolding ##
def unfold(p, diag=False):

    # Suppress divide by 0 warnings
    inv_err = geterr()['invalid']
    seterr(invalid = 'ignore')

    # Build composition list
    eList = []
    for key in p.keys():
        if key[1] not in eList:
            eList += [key[1]]
    eTable = [[r, t] for r in eList for t in eList]

    # Bayesian unfolding
    for r, t in eTable:
        rr, tt = 'R'+r, 'T'+t
        p[tt+'|'+rr] = ((p[rr+'|'+tt] * p[tt]).transpose() /
            sum(sum([p[rr+'|T'+e]*p['T'+e] for e in eList], axis=0), axis=1))

    # Option to replace the probability tables with diagonal matrices
    # Used for testing purposes
    if diag:
        keyList = [key for key in p.keys() if '|' in key]
        for key in keyList:
            p[key] = ones(p[keyList[0]].shape)
            p[key] = diag(diag(p[key]))

    # Normalize
    for t in eList:
        pTable = [p['T'+t+'|R'+e]*p['R'+e] for e in eList]
        p['T'+t] = nansum(nansum(pTable, axis=0), axis=1)
    sumprobs = sum([p['T'+e] for e in eList])
    for t in eList:
        p['T'+t] /= sumprobs

    seterr(invalid = inv_err)

    return p


##===========================================================================##
## Applying unfolding to simulation

## Plot the counts after unfolding ##
def counts(s, niter, sDict={'':[-0.25, [[7.5, -0.75]]]}, 
                        smooth=True, spl=False, out=False):

    r = log10(s['ML_energy'])
    Emids = getEmids()
    cutName = 'llh'
    cut = s['cuts'][cutName]
    eList = ['p', 'f']
    compDict = {'p':'P', 'f':'Fe', 'All':'All'}
    colorDict = {'p':'b', 'h':'y', 'o':'c', 'f':'r', 'All':'k'}

    f = plt.figure()
    ax1 = f.add_subplot(1,1,1)
    ax1.set_title('Energy Spectum using '+cutName+' cut')
    ax1.set_xlabel('Log10(Energy/GeV)')
    ax1.set_ylabel('Counts')

    # Load probability tables
    p = getProbs(s, cut)
    st = len(Emids) - len(p['Rf|Tf'][0])
    Emids = Emids[st:]

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
        temp[key] = powerFit.powerFit(s0, sTable=sTable, st=st)
    specCut = fakeSpec(s, cut, temp)
    specCut = array([True for i in range(len(specCut))])

    # Create starting arrays
    N_passed = Nfinder(r, cut*specCut).sum()
    N = {}
    for e in eList:
        N[e] = Nfinder(r, cut*s[e]*specCut)
        p['R'+e] = N[e] / N_passed
        p['T'+e] = powerFit.powerFit(-2.7, st=st)

    # Get relative errors due to unfolding
    #chi2, unrel = load('unfold_err.npy')
    # Due to efficiency
    effarea, sigma, relerr = sim.getEff(s, cut, smooth=smooth)
    relerr = relerr[st:]

    # Bayesian unfolding
    for i in range(niter):
        p = unfold(p)
        #All_unfold = (p['Tf']+p['Tp']) * N_passed
        #ax1.errorbar(Emids, All_unfold, yerr=unrel[i]*All_unfold, label=str(i))
        # Smooth prior before next iteration (except last time)
        if i < niter-1:
            for e in eList:
                p['T'+e] = smoother(p['T'+e])

    # Find bin values
    Nun, Rel, Err = {},{},{}
    for e in eList:
        Nun[e] = p['T'+e] * N_passed
    Nun['All'] = sum([Nun[e] for e in eList], axis=0)

    # Calculate errors
    for e in eList + ['All']:
        ## NOTE: I don't think you can use the relative errors like this ##
        #Rel[e] = sqrt(1/Nun[e] + relerr**2 + unrel[niter-1]**2) 
        Rel[e] = sqrt(1/Nun[e] + relerr**2) 
        Err[e] = Nun[e] * Rel[e]
        # And plot
        pnt = colorDict[e] + '.'
        ax1.errorbar(Emids, Nun[e], yerr=Err[e], fmt=pnt, label=e)

    # Attempt to fit with broken power law
    #p0 = [1, 10**(-4.5), 7.5, -0.5]
    #yfit = powerFit.pow_fit(10**Emids, Nun['f'], Err['f'], p0)
    #ax1.plot(Emids, yfit, label='test')

    #err = {}
    #err[0] = 1/sqrt(All_Nunfold)
    #err[1] = relerr
    #err[2] = unrel[niter-1]

    # Plot error bars nicely
    #errCalc = lambda i: sqrt(sum([err[j]**2 for j in range(i+1)], axis=0))
    #uperr, dnerr = {}, {}
    #for i in range(len(err)):
    #    uperr[i] = All_Nunfold * (1 + errCalc(i))
    #    dnerr[i] = All_Nunfold * (1 - errCalc(i))
    #    for j in range(len(dnerr[i])):
    #        if dnerr[i][j] < 0:
    #            dnerr[i][j] = 10**(-3)

    #ax1.errorbar(Emids, All_Nunfold, yerr=All_err, fmt='gx', label='Unfold')
    #ax1.plot(Emids, All_Nunfold, 'b.', label='Unfold')
    #for i in range(len(err)):
    #    ax1.fill_between(Emids, dnerr[i], uperr[i], facecolor='blue', alpha=0.2)

    # Plot true spectrum
    MC = {}
    t = log10(s['MC_energy'])
    MC['All'] = Nfinder(t, cut*specCut)[st:]
    for e in eList:
        MC[e] = Nfinder(t, cut*specCut*s[compDict[e]])[st:]
    for e in eList + ['All']:
        pnt = colorDict[e] + 'x'
        ax1.plot(Emids, MC[e], pnt, label='MC_'+compDict[e])

    # plot original (not unfolded) spectrum
    O_N = (sum([N[e] for e in eList], axis=0))[st:]
    O_relerr = sqrt(1/O_N + relerr**2)
    O_err = O_N * O_relerr
    ax1.errorbar(Emids, O_N, yerr=O_err, fmt='gx', label='Original')

    ax1.set_yscale('log')
    ax1.legend(loc='lower left')
    #ax1.set_ylim((10**(-1), 10**(4)))
    if out:
        plt.savefig('collab/pics/'+out+'.png')
    plt.show()






