#!/usr/bin/env python

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import sys, os, argparse, analysis
import myGlobals as my

def writeDist(config, infoFile, overwrite=False):

    # Load simulation
    s = analysis.load_sim(config)
    energy  = np.log10(s['MC_energy'])
    weights = analysis.getWeights(s)

    # Read in existing distribution file
    lines = []
    if os.path.isfile(infoFile):
        with open(infoFile, 'r') as f:
            lines = f.readlines()
    lines = [line.strip() for line in lines]
    table = [line.split(' ') for line in lines]
    paramList = [i[:3] for i in table]

    # Info stored in [config, nmin, nmax, median energy, sigL, sigR, counts]
    nList = range(3, 16)
    nList.append(100)
    newParams = [[config, i, j] for i in nList for j in nList if i < j]

    for config, nmin, nmax in newParams:

        # Check to see if the params exists already
        found = True
        try:  i = paramList.index([config, str(nmin), str(nmax)])
        except ValueError:
            found = False
        # If it exists, option to overwrite
        if found and not overwrite:
            print 'Info for', config, nmin, nmax, 'already exists...'
            continue
        if found and overwrite:
            table.remove(table[i])

        print 'Working on', config, nmin, nmax
        # Calculate information for energy range
        cut = (s['nstation'] >= nmin) * (s['nstation'] < nmax)
        values = energy[cut]
        w = weights[cut]
        # Weighted median and 68% containment
        median = analysis.weighted_percentile(values, w, 50)
        sigL   = median - analysis.weighted_percentile(values, w, 16)
        sigR   = analysis.weighted_percentile(values, w, 84) - median
        # Unbiased weighted variance
        mean = np.average(values, weights=w)
        var  = (sum(w * (values-mean)**2) /
               (sum(w) - sum(w**2)/sum(w)))
        counts = len(values)

        # Append to file
        newLine = [config, nmin, nmax, median, sigL, sigR, var, counts]
        newLine = ['%s' % i for i in newLine]
        table.append(newLine)

        # Write as you go
        lines = [' '.join(line) for line in table]
        lines = [line + '\n' for line in lines]
        lines.sort()

        with open(infoFile, 'w') as f:
            f.writelines(lines)


if __name__ == "__main__":

    # Setup global path names
    my.setupAnisotropy(verbose=False)

    p = argparse.ArgumentParser(
            description='Writes energy distribution information to txt file')
    p.add_argument('-c', '--config', dest='config', nargs='*',
            help='Detector configuration for simulation')
    p.add_argument('-o', '--out', dest='out',
            help='Name for output file')
    p.add_argument('--overwrite', dest='overwrite',
            default=False, action='store_true',
            help='Option to overwrite existing destination file')
    args = p.parse_args()

    if not args.out:
        args.out = '%s/sim/eDist_IT.txt' % (my.ani_data)

    for config in args.config:
        writeDist(config, args.out, overwrite=args.overwrite)

