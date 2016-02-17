#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

import myGlobals as my
from load_sim import load_sim
from llhtools import *
from useful import getMids

""" Calculate probability tables for a given cut  """
def getProbs(s, cut, reco=True, zcorrect=False):

    from zfix import zfix
    print 'Calculating probability...'
    q = {}
    with np.errstate(invalid='ignore'):
        ereco = np.log10(s['ML_energy'])
        etrue = np.log10(s['MC_energy'])
    if zcorrect:
        z = np.pi - s['zenith']
        ereco -= zfix(z, bintype='logdist')

    # Option to limit both axes to reconstructed energy range
    tEbins = getEbins(reco=reco)
    rEbins = getEbins(reco=True)

    # Calculate probabilities
    rList = getComps(s)
    tList = getComps(s, reco=False)
    trTable = [[t, r] for t in tList for r in rList]
    for t, r in trTable:
        recocut = s['llh_comp'] == r
        truecut = s['comp'] == t
        c0 = cut * truecut * recocut
        RvT, x, y = np.histogram2d(ereco[c0], etrue[c0], bins=(rEbins,tEbins))
        q['RvT_'+t+r] = np.asarray(RvT, dtype='float')

    # Normalize
    Ntot = {}
    for t in tList:
        Ntot[t] = np.sum([q['RvT_'+t+r] for r in rList], axis=(0,1))
    for t, r in trTable:
        with np.errstate(invalid='ignore'):
            q['RvT_'+t+r] /= Ntot[t]

    # Store in probability dictionary
    # p['Rr|Tt'][i][j] = probability comp t in bin j is reconstructed as r in i
    #  - sum each column to get 1
    prob = {}
    for t, r in trTable:
        prob['R'+r+'|T'+t[0].lower()] = q['RvT_'+t+r]

    # Restrict reconstructed energy range to where we have counts
    #N_reco = Nfinder(ereco, cut)
    #min_bin = np.nonzero(N_reco)[0][0]
    #for key in prob.keys():
    #    prob[key] = prob[key][min_bin:]

    # Original cut - this was working but I don't know why...
    # Setup so N_true for lowest energy is >= N_true for highest energy
    #N_max = Nfinder(etrue, cut)[-1]
    #min_bin = np.nonzero(Nfinder(etrue, cut) > N_max)[0][0]
    #for key in prob.keys():
    #    prob[key] = prob[key][min_bin:]

    return prob


""" Get probability tables fast for preset cuts """
def getProbs_fast(zcorrect=False):
    my.setupShowerLLH(verbose=False)
    pFile = '%s/ProbTables.npy' % my.llh_resource
    if zcorrect:
        pFile = pFile.replace('.npy','_zfix.npy')
    p = np.load(pFile)
    p = p.item()
    return p


""" Look at our probability space """
def plotter(s, cut, name, reco=False, zcorrect=False):

    rEbins = getEbins(reco=True)
    tEbins = getEbins(reco=reco)

    fig, ax = plt.subplots()
    ax.set_title(name)
    ax.set_xlabel(r'$\log_{10}(E_{\mathrm{true}}/\mathrm{GeV})$')
    ax.set_ylabel(r'$\log_{10}(E_{\mathrm{reco}}/\mathrm{GeV})$')

    p = getProbs(s, cut, reco=reco, zcorrect=zcorrect)

    #st = len(Emids) - len(p[name])
    x, y = tEbins, rEbins
    X, Y = np.meshgrid(x, y)
    with np.errstate(divide='ignore'):
        Z = np.log10(p[name])
    c0 = (Z==Z) * (abs(Z) < np.inf)
    vmin, vmax = Z[c0].min(), Z[c0].max()
    pc = ax.pcolor(X, Y, Z, vmin=vmin, vmax=vmax)
    cb = fig.colorbar(pc, ax=ax)
    ax.set_xlim(tEbins.min(), tEbins.max())
    ax.set_ylim(rEbins.min(), rEbins.max())
    plt.show()


""" Look at our probability space """
def p2(p, name, reco=False):

    if name[0] == 'R':
        y = getEbins(reco=True)
        x = getEbins(reco=reco)
    if name[0] == 'T':
        x = getEbins(reco=True)
        y = getEbins(reco=reco)

    fig, ax = plt.subplots()
    ax.set_title(name)
    if name[0] == 'R':
        ax.set_xlabel(r'$\log_{10}(E_{\mathrm{true}}/\mathrm{GeV})$')
        ax.set_ylabel(r'$\log_{10}(E_{\mathrm{reco}}/\mathrm{GeV})$')
    if name[0] == 'T':
        ax.set_xlabel(r'$\log_{10}(E_{\mathrm{reco}}/\mathrm{GeV})$')
        ax.set_ylabel(r'$\log_{10}(E_{\mathrm{true}}/\mathrm{GeV})$')

    #st = len(Emids) - len(p[name])
    X, Y = np.meshgrid(x, y)
    with np.errstate(divide='ignore'):
        Z = np.log10(p[name])
    c0 = (Z==Z) * (abs(Z) < np.inf)
    vmin, vmax = Z[c0].min(), Z[c0].max()
    pc = ax.pcolor(X, Y, Z, vmin=vmin, vmax=vmax)
    cb = fig.colorbar(pc, ax=ax)
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(y.min(), y.max())
    plt.show()



if __name__ == "__main__":

    # Setup global variables
    my.setupShowerLLH(verbose=False)

    print 'Loading simulation...'
    s = load_sim(bintype='logdist')

    print 'Calculating probabilities...'
    p = getProbs(s, s['cuts']['llh'])
    outFile = '%s/ProbTables.npy' % my.llh_resource
    np.save(outFile, p)

    p = getProbs(s, s['cuts']['llh'], zcorrect=True)
    outFile = '%s/ProbTables_zfix.npy' % my.llh_resource
    np.save(outFile, p)


