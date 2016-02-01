#!/usr/bin/env python

#############################################################################
# Converts individual count tables into a single normalized LLH table
#############################################################################

import glob, sys
from numpy import *

if __name__ == "__main__":

    config = 'IT81'

    resourcedir = '/net/user/zgriffith/ShowerLLH/resources/'
    simDict, orig, norm = {},{},{}
    simDict = dict({'7351':'P', '7006':'P', '7579':'P'})
    simDict.update({'7483':'He','7241':'He','7263':'He','7791':'He'})
    simDict.update({'7486':'O', '7242':'O', '7262':'O', '7851':'O'})
    simDict.update({'7394':'Fe','7007':'Fe','7784':'Fe'})
    simDict.update({'10687':'G'})

    # Get the bins and their sizes
    n = {}
    binDict = load(resourcedir + 'ShowerLLH_bins.npy')
    binDict = binDict.item()
    for key in binDict.keys():
        n[key[0]]  = len(binDict[key]) - 1
    binLengths = (n['E'], n['Z'], n['S'], n['D'], n['C'])

    EZSDList = [[e, z, s, d] for e in range(n['E']) for z in range(n['Z']) \
                                for s in range(n['S']) for d in range(n['D'])]

    fileList = glob.glob(resourcedir + config+'/CountTable_*.npy')
    fileList.sort()

    for file in fileList:
        print 'Loading', file
        h = load(file)

        st = file.rfind('_') + 1
        end = file.rfind('.')
        sim = file[st:end]
        e = simDict[sim]
        if e not in orig.keys():
            orig[e] = zeros(binLengths, dtype=float)

        orig[e] += h

    for e in orig.keys():
        orig[e] += .1    # baseline so zeros don't give errors
        norm[e] = zeros(binLengths, dtype=float)
        for EZSD in EZSDList:
            # normalize and log table (log10(Hist/Hist.sum()))
            E,Z,S,D = EZSD
            norm[e][E,Z,S,D] =log10((orig[e][E,Z,S,D])/(orig[e][E,Z,S,D]).sum())

    # Write to file
    outFile = resourcedir + 'LLHTables.npy'
    save(outFile, norm)


