#!/usr/bin/env python

import glob
import matplotlib.pyplot as plt
import numpy as np

def plotter(comp, x='E', y='D', binVals={}, count=False):

    if comp == 'G':
        prefix = '/net/user/zgriffith/ShowerLLH/resources'
    else:
        prefix = '/net/user/fmcnally/ShowerLLH/resources'

    compDict = {'P':['7006','7351','7579'], 'Fe':['7007','7394','7784'], 'G':['10687']}
    simList = compDict[comp]
    binDict = np.load('%s/ShowerLLH_bins.npy' % prefix)
    binDict = binDict.item()

    if not count:
        # Load llh tables
        d = np.load('%s/LLHTables.npy' % prefix)
        d = d.item()
        d = d[comp]

    if count:
        # Load count tables
        for i, sim in enumerate(simList):
            if comp == 'G':
                temp = np.load('%s/IT81/CountTable_%s.npy' % (prefix, sim))
            else:
                temp = np.load('%s/IT73/CountTable_%s.npy' % (prefix, sim))
            if i == 0:
                d = np.zeros(temp.shape)
            d += temp

    # Collapse array to desired size
    names = ['E','Z','S','D','C']
    labels = ['log(Energy) (GeV)', 'Zenith Angle', 'Snow Table',
              'Distance from Shower Core (m)', 'Charge (VEM)']
    fig, ax = plt.subplots()
    plt.xlabel(labels[names.index(x)])
    plt.ylabel(labels[names.index(y)])
    i, j = 0, len(names)
    index = 0
    while i < j:
        # Select desired bins
        if names[i] in binVals.keys():
            binCut = [slice(None, None, None) for k in range(j)]
            binCut[i] = binVals[names[i]]
            d = d[binCut]
            names.remove(names[i])
            j -= 1
        # Collapse other bins not needed
        elif names[i] not in [x, y]:
            d = d.sum(axis=i)
            names.remove(names[i])
            j -= 1
        else:
            i += 1

    # Order as desired
    wanted = [x, y]
    if names != wanted:
        names = names.reverse()
        d = d.T

    # Setup plot
    xbins, ybins = binDict[x+'bins'], binDict[y+'bins']
    # Eliminate infinities
    for bins in [xbins, ybins]:
        if bins[-1] == np.inf:
            bins[-1] = bins[-2] + (bins[-2]-bins[-3])
    X = (xbins[:-1]+xbins[1:])/2.
    Y = (ybins[:-1]+ybins[1:])/2.
    if y=='C':
        Y = np.log10(Y)
    print d.shape
    Z = d.T
    Z[Z==0] += 1
    if count:
        Z = np.log10(Z)
    print Z.shape
    pc = ax.pcolor(X, Y, Z)
    cb = fig.colorbar(pc, ax=ax)
    plt.show()

    
