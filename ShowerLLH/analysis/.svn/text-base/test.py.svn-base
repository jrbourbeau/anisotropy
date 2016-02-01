#!/usr/bin/env python
import numpy as np

def test(d, nbins):

    binList = [0]
    i, j, count = 0,0,0
    nstation = d['NStations']
    w = (3*d['weights']['w3'] + 1.5*d['weights']['w1.5'] + d['weights']['w1'])
    ntot = w.sum()
    #ntot = float(len(nstation))
    print 'Bin size:', ntot/nbins

    while i < nstation.max():
        cut = (nstation >= j) * (nstation < i)
        count = w[cut].sum()
        if count > (ntot/nbins):
            binList += [i]
            j = i
        i += 1

    binList += [np.inf]

    return binList

