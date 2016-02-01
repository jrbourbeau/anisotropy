#!/usr/bin/env python

from numpy import *
from useful import *
import matplotlib.pyplot as plt

""" Calculate probability tables for a given cut  """
def getProbs(s, cut):

    print 'Calculating probability...'
    seterr(invalid='ignore')
    q = {}
    reco = log10(s['ML_energy'])
    true = log10(s['MC_energy'])
    Ebins = getEbins()

    # Calculate probabilities
    tList = ['P', 'Fe']
    rList = [t[0].lower() for t in tList]
    trTable = [[t, r] for t in tList for r in rList]
    for t, r in trTable:
        c0 = cut * s[t] * s[r]
        RvT, x, y = histogram2d(reco[c0], true[c0], bins=(Ebins,Ebins))
        q['RvT_'+t+r] = asarray(RvT, dtype='float')

    # Normalize
    Ntot = {}
    for t in tList:
        Ntot[t] = sum(sum([q['RvT_'+t+r] for r in rList], axis=0), axis=0)
    for t, r in trTable:
        q['RvT_'+t+r] /= Ntot[t]

    # Store in probability dictionary
    # p['Rr|Tt'][i][j] = probability comp t in bin j is reconstructed as r in i
    #  - sum each column to get 1
    prob = {}
    for t, r in trTable:
        prob['R'+r+'|T'+t[0].lower()] = q['RvT_'+t+r]

    # Setup so N_true for lowest energy is >= N_true for highest energy
    N_max = Nfinder(true, cut)[-1]
    min_bin = nonzero(Nfinder(true, cut) > N_max)[0][0]
    for key in prob.keys():
        prob[key] = prob[key][:,min_bin:]

    seterr(invalid='print')
    return prob


""" Look at our probability space """
def plotter(s, cut, name):

    r = log10(s['ML_energy'])
    Ebins, Emids = getEbins(), getEmids()

    f = plt.figure()
    ax1 = f.add_subplot(1,1,1)
    ax1.set_title(name)
    ax1.set_xlabel('true')
    ax1.set_ylabel('reco')

    p = getProbs(s, cut)
    st = len(Emids) - len(p[name][0])

    x, y = Emids, Emids[st:]
    X, Y = meshgrid(x, y)
    Z = log10(p[name])
    im = ax1.imshow(Z, interpolation='bilinear',
            origin='lower', extent=[Ebins[st], 9.5, Ebins[0], 9.5])
    f.colorbar(im)
    plt.show()


if __name__ == "__main__":

    print 'Loading simulation...'
    vars = {}
    execfile('load_sim.py', vars)
    s = vars['s']
    resourcedir = '/net/user/fmcnally/ShowerLLH/resources/'

    print 'Calculating probabilities...'
    p = getProbs(s, s['cuts']['llh'])
    outFile = resourcedir + 'ProbTables.npy'
    save(outFile, p)
