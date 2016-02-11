#!/usr/bin/env python

import numpy as np
import glob

from useful import getMids
from saveMedian import histPercentile, histVar

def histMedian(h, bins):

    # Assume you only want to operate on the last dimension
    nx, ny = h.shape
    median, sigL, sigR, var = np.zeros((4,nx))
    binMids = getMids(bins, infvalue=100)
    for i in range(nx):
        median[i] = histPercentile(h[i], 50, binMids)
        sigL[i]   = histPercentile(h[i], 16, binMids)
        sigR[i]   = histPercentile(h[i], 84, binMids)
        var[i]    = histVar(h[i], binMids)

    return median, sigL, sigR, var


def percentOverlap():

    files = glob.glob('IC*_hists.npy')
    files.sort()

    for f in files:

        h = np.load(f)
        nbins = 3
        bins = np.arange(2.75, 9.01, 0.05)
        binMids = getMids(bins, infvalue=100)
        print f
        median, sigL, sigR, var = histMedian(h, bins)
        for i in range(1, len(h)):
            print 'Bin %i max = %f' % (i-1, sigR[i-1])
            ecut = (binMids > sigR[i-1])
            myfrac = float(h[i][ecut].sum()) / h[i].sum()
            print 'Fraction of events in bin %i > u.l.: %f' % (i, myfrac)
        #print 'Bin 0: %i, Bin 1: %i, Bin2: %i'%(sum(h[0]),sum(h[1]),sum(h[2]))
        #print '0/1 overlap: %i' % np.amin([h[0],h[1]], axis=0).sum()
        #print '1/2 overlap: %i' % np.amin([h[1],h[2]], axis=0).sum()


if __name__ == "__main__":

    #percentOverlap()
    f = 'IC86-III_hists.npy'
    h = np.load(f)
    textTable = []
    for row in h:
        textTable += [', '.join(['%s' % r for r in row]) + '\n']

    with open('spencerTable.txt', 'w') as f:
        f.writelines(textTable)

