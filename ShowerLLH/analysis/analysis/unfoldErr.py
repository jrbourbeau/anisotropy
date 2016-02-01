#!/usr/bin/env python

from numpy import *
import powerFit, prob, load_data, load_sim
from useful import *
from unfold import unfold
import matplotlib.pyplot as plt

##===========================================================================##
## Useful functions

""" Find the chi2 of our data through unfolding """
def getChi2(s, d, k=20, n=1000, test=False, smooth=True):

    # Setup initial information
    print 'Calculating unfolding variance...'
    #r = log10(s['ML_energy'])
    r = log10(d['ML_energy'])
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
    nTable = zeros((n,k,l2))
    for e in eList:
        Nnew[e] = random.poisson(N[e], (n, len(N[e])))
    # Total number of events for each "true" distribution
    totN = sum([Nnew[e] for e in eList], axis=(0,2)).astype('float')

    # Unfold every wiggled distribution k times
    for i in range(n):
        for e in eList:
            p['R'+e] = Nnew[e][i] / totN[i]
            p['T'+e] = powerFit.powerFit(-2.7, [[6, -3.0]], st=st)

        for j in range(k):
            p = unfold(p)
            nTable[i][j] = sum([p['T'+e] for e in eList], axis=0) * totN[i]
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

        counts = sum([N[e] for e in eList], axis=0)
        newcts = sum([Nnew[e] for e in eList], axis=0)
        print 'Counts:', (counts)[test]
        print 'Variance:', sqrt((counts)[test])
        print 'Variance (calc):', sqrt(var((newcts).transpose()[test]))
        print 'Mean:', mean((newcts).transpose()[test])
        ax1.hist((newcts).transpose()[test], bins=100, color='red')
        test -= st
        counts = (nTable.transpose()[test])
        for j in range(k):
            print 'Var:', sqrt(var(counts[j]))
            print 'Mean:', mean(counts[j])
            ax1.hist(counts[j], bins=100, color='blue', label=str(j))

        plt.show()
        return

    # Calculate the variance using nTable
    sigma = sqrt(var(nTable, axis=0))
    ## ALWAYS GIVING A 0 IN A SPECIFIC BIN. WHY IS THIS HAPPENING
    #print nTable
    #print nTable.shape

    # Unfold the original counts
    chi2 = zeros(k)
    relerr = []
    old = {}
    for e in eList:
        p['R'+e] = N[e] / Ntot
        p['T'+e] = powerFit.powerFit(-2.7, st=st)
        old[e] = p['T'+e]

    seterr(invalid='ignore')
    for j in range(1, k+1):
        p = unfold(p)
        oldsum = sum([old[e] for e in eList], axis=0)
        new = sum([p['T'+e] for e in eList], axis=0)
        chi2[j-1] = (1./(l2+1) * nansum(((new - oldsum)*Ntot / sigma[j-1])**2))
        relerr += [(sigma[j-1] / (new*Ntot))]

        for e in eList:
            old[e] = p['T'+e]
            if smooth:
                p['T'+e] = smoother(p['T'+e])

    return chi2, relerr


##===========================================================================##
## Write to output if called

if __name__ == "__main__":

    # Load simulation
    prefix = '/home/fmcnally/ShowerLLH/analysis/'
    config = 'IT73'
    s = load_sim.load_sim(config)

    # Load data
    d = load_data.load_data(config, '201006')

    chiDict = {}
    chi2, relerr = getChi2(s, d, k=30, n=10000, test=False)
    chiDict['chi2'] = chi2
    chiDict['relerr'] = relerr
    outFile = '/net/user/fmcnally/ShowerLLH/resources/unfold_err.npy'
    save(outFile, chiDict)
















