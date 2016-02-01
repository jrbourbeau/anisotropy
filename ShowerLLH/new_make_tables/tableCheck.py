#!/usr/bin/env python

import glob, argparse
import numpy as np
import matplotlib.pyplot as plt

import myGlobals as my
from useful import getMids

def counts(d, param, logx, logy):

    binList = [d['bins'][i][0] for i in d['bins']]
    idx = binList.index(param)
    axes = [i for i in d['bins']]
    axes.remove(idx)
    axes = tuple(axes)

    bins = d['bins'][idx][1]
    x = getMids(bins)
    y = d['counts'].sum(axis=axes)

    fig, ax = plt.subplots()
    ax.plot(x, y, '.')
    if logx:
        ax.set_xscale('log')
    if logy:
        ax.set_yscale('log')
    plt.show()


if __name__ == "__main__":

    my.setupShowerLLH(verbose=False)

    p = argparse.ArgumentParser(
            description='Plot 1D versions of CountTables')
    p.add_argument('-s', '--sim', dest='sim', nargs='*',
            help='Simulation to load (default=all)')
    p.add_argument('-b', '--bintype', dest='bintype',
            default='standard',
            choices=['standard','nozenith','logdist'],
            help='Desired binning')
    p.add_argument('-p', '--param', dest='param',
            default='E',
            choices=['E','Z','S','D','C'],
            help='Parameter you want to view')
    p.add_argument('--logx', dest='logx',
            default=False, action='store_true',
            help='Log x-axis')
    p.add_argument('--logy', dest='logy',
            default=False, action='store_true',
            help='Log y-axis')
    args = p.parse_args()

    files = glob.glob('%s/CountTables/CountTable_*.npy' % my.llh_resource)
    files = [f for f in files if 'Part' not in f]
    files = [f for f in files if args.bintype in f]
    if args.sim != None:
        files = [f for f in files if any([s in f for s in args.sim])]

    for i, file in enumerate(files):
        q = np.load(file)
        q = q.item()
        if i == 0:
            d = q
        else:
            d['counts'] += q['counts']

    counts(d, args.param, args.logx, args.logy)
    
