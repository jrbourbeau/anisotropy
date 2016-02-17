#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

import myGlobals as my
import powerFit, prob
from llhtools import *
from unfold import unfold

##===========================================================================##
## Useful functions

""" Find the chi2 of our data through unfolding """
def getChi2(s, d, k=20, n=1000, test=False, smooth=True):

    # Setup initial information
    print 'Calculating unfolding variance...'
    #r = np.log10(s['ML_energy'])
    r = np.log10(d['ML_energy'])
    Emids = getEmids()
    cutName = 'llh'
    dcut, scut  = d['cuts'][cutName], s['cuts'][cutName]
    eList = ['p', 'f']

    # Load probability tables
    p = prob.getProbs(s, scut)
    st = len(Emids) - len(p['Rf|Tf'][0])
    Emids = Emids[st:]
    l2 = len(Emids)

    # Create a fake spectrum
    #s0 = -0.5
    #sDict = {'':powerFit.powerFit(s0, st=st)}
    #specCut = fakeSpec.fakeSpec(s, scut, sDict)
    #for e in eList:
    #    N[e] = Nfinder(r, scut*s[e]*specCut)

    # Calculate counts for data distribution
    N = {}
    for e in eList:
        N[e] = Nfinder(r, dcut*d[e])
    Ntot = sum([N[e] for e in eList])

    # Build a table of n "true" distributions after 'wiggling' observed
    # distribution in poisson errors
    Nnew = {}
    nTable = np.zeros((n,k,l2))
    for e in eList:
        Nnew[e] = np.random.poisson(N[e], (n, len(N[e])))
    # Total number of events for each "true" distribution
    totN = np.sum([Nnew[e] for e in eList], axis=(0,2)).astype('float')

    # Unfold every wiggled distribution k times
    for i in range(n):
        for e in eList:
            p['R'+e] = Nnew[e][i] / totN[i]
            p['T'+e] = powerFit.powerFit(-2.7, [[6, -3.0]], st=st)

        for j in range(k):
            p = unfold(p)
            nTable[i][j] = np.sum([p['T'+e] for e in eList], axis=0) * totN[i]
            if smooth:
                for e in eList:
                    p['T'+e] = smoother(p['T'+e])

    # Option for analyzing the behavior in a single test bin
    if test:

        f = plt.figure()
        ax1 = f.add_subplot(1,1,1)
        ax1.set_title('Distribution for log10(E/GeV) =' + str(Emids[test]))
        ax1.set_xlabel('Counts')
        ax1.set_ylabel('')

        counts = np.sum([N[e] for e in eList], axis=0)
        newcts = np.sum([Nnew[e] for e in eList], axis=0)
        print 'Counts:', (counts)[test]
        print 'Variance:', np.sqrt((counts)[test])
        print 'Variance (calc):', np.sqrt(np.var((newcts).transpose()[test]))
        print 'Mean:', np.mean((newcts).transpose()[test])
        ax1.hist((newcts).transpose()[test], bins=100, color='red')
        test -= st
        counts = (nTable.transpose()[test])
        for j in range(k):
            print 'Var:', np.sqrt(np.var(counts[j]))
            print 'Mean:', np.mean(counts[j])
            ax1.hist(counts[j], bins=100, color='blue', label=str(j))

        plt.show()
        return

    # Calculate the variance using nTable
    sigma = np.sqrt(np.var(nTable, axis=0))
    ## ALWAYS GIVING A 0 IN A SPECIFIC BIN. WHY IS THIS HAPPENING
    #print nTable
    #print nTable.shape

    # Unfold the original counts
    chi2 = np.zeros(k)
    relerr = []
    old = {}
    for e in eList:
        p['R'+e] = N[e] / Ntot
        p['T'+e] = powerFit.powerFit(-2.7, st=st)
        old[e] = p['T'+e]

    seterr(invalid='ignore')
    for j in range(1, k+1):
        p = unfold(p)
        oldsum = np.sum([old[e] for e in eList], axis=0)
        new = np.sum([p['T'+e] for e in eList], axis=0)
        chi2[j-1] = (1./(l2+1) * np.nansum(((new - oldsum)*Ntot / sigma[j-1])**2))
        relerr += [(sigma[j-1] / (new*Ntot))]

        for e in eList:
            old[e] = p['T'+e]
            if smooth:
                p['T'+e] = smoother(p['T'+e])

    return chi2, relerr


##===========================================================================##
## Write to output if called

if __name__ == "__main__":

    # Setup global variables
    my.setupShowerLLH(vebose=False)
    from load_sim import load_sim
    from load_data import load_data

    # Load simulation and data
    s = load_sim('IT73')
    d = load_data(config, '201006')

    chiDict = {}
    chi2, relerr = getChi2(s, d, k=30, n=10000, test=False)
    chiDict['chi2'] = chi2
    chiDict['relerr'] = relerr

    outFile = '%s/unfold_err.npy' % my.llh_resource
    np.save(outFile, chiDict)
















