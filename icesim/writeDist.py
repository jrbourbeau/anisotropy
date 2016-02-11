#!/usr/bin/env python

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import os, glob, argparse
from analysis import getSimBase, getSplineEnergies
from saveMedian import histPercentile, histVar, getWeights
import myGlobals as my

if __name__ == "__main__":

    # Setup global path names
    my.setupAnisotropy(verbose=False)

    # Argument parser
    p = argparse.ArgumentParser(
            description="""Calculates energy distribution information for given
                           cuts and writes to file.""")
    p.add_argument('-c', '--config', dest='config',
            help='Detector configuration')
    p.add_argument('-e', '--elist', dest='eList', nargs='*', type=float,
            help='Energy values to bin in')
    p.add_argument('-o', '--outFile', dest='outFile',
            help='Option for output file name')
    p.add_argument('--overwrite', dest='overwrite',
            default=False, action='store_true',
            help='Option to overwrite existing energy distribution info')
    args = p.parse_args()

    # Default values
    if not args.outFile:
        args.outFile = '%s/eDist_IC.txt' % my.ani_sim
    if not args.eList:
        args.eList = np.arange(4, 8.1, .25)
        args.eList = np.insert(args.eList, 0, 0)
        args.eList = np.append(args.eList, 100)

    # Energy bins for histograms within reco energy bins
    estep = 0.05
    ebins = np.arange(2.75, 9+estep/2., estep)
    ebins = np.append(ebins, [0, 100])
    ebins.sort()
    ehists = np.zeros((len(args.eList)-1, len(ebins)-1))

    # Load simulation
    simBase = getSimBase(args.config)
    simList = ['%s.npy' % simBase]
    if not os.path.isfile(simList[0]):
        simList = glob.glob('%s_Part*.npy' % simBase)
    simList.sort()

    # Bad run file
    badRunFile = '%s/badFiles.txt' % (my.ani_sim)

    for sim in simList:

        print 'Loading', sim
        d = np.load(sim)
        d = d.item()

        # Only store necessary information
        if 'reco2' in d.keys():
            theta, phi = hp.pix2ang(1024, d['reco2'].astype('int'))
            d['dstZenith'] = theta
        x = np.cos(d['dstZenith'])
        y = np.log10(d['nchannel'])
        e = np.log10(d['energy'])
        weights = getWeights(d, args.config, badRunFile)

        # Apply cut
        qcut = (x >= 0.3)
        x, y, e = x[qcut], y[qcut], e[qcut]
        if weights != None:
            weights = weights[qcut]

        print 'Calculating spline energies...'
        esplines = getSplineEnergies(args.config, x, y)

        # Fill histogram
        ehists += np.histogramdd(np.transpose([esplines, e]), 
                bins=[args.eList,ebins], weights=weights)[0]

    # Read in existing distribution file
    lines = []
    if os.path.isfile(args.outFile):
        with open(args.outFile, 'r') as f:
            lines = f.readlines()
    lines = [line.strip() for line in lines]
    table = [line.split(' ') for line in lines]
    paramList = [i[:3] for i in table]

    # Info stored in [config, emin, emax, median energy, sigL, sigR, counts]
    newParams = [[args.config, emin, emax] \
            for emin in args.eList for emax in args.eList if emin < emax]

    for config, emin, emax in newParams:

        """ This doesn't seem to be working """
        # Check to see if the params exists already
        found = True
        try:  i = paramList.index([config, str(emin), str(emax)])
        except ValueError:
            found = False
        # If it exists, option to overwrite
        if found and not args.overwrite:
            print 'Info for', config, emin, emax, 'already exists...'
            continue
        if found and args.overwrite:
            table.remove(table[i])

        # Calculate information for energy range
        idx0 = np.where(args.eList==emin)[0][0]
        idx1 = np.where(args.eList==emax)[0][0]
        h = ehists[idx0:idx1].sum(axis=0)
        emids  = (ebins[:-1] + ebins[1:]) / 2.
        median = histPercentile(h, 50, emids)
        sigL   = histPercentile(h, 16, emids)
        sigR   = histPercentile(h, 84, emids)
        var = histVar(h, emids)
        counts = h.sum()

        # Append to file
        newLine = [config, emin, emax, median, sigL, sigR, var, counts]
        newLine = ['%s' % i for i in newLine]
        table.append(newLine)

        # Write as you go
        lines = [' '.join(line) for line in table]
        lines = [line + '\n' for line in lines]
        lines.sort()

        with open(args.outFile, 'w') as f:
            f.writelines(lines)




