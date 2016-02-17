#!/usr/bin/env python

import numpy as np

import myGlobals as my
from llhtools import inPoly

if __name__ == "__main__":

    # Load information
    my.setupShowerLLH(verbose=False)
    infile = '%s/IT73_data/DataPlot_201006_logdist_full.npy' % my.llh_data
    d = np.load(infile)
    d = d.item()

    c = {}
    c[0] = (d['ML_energy'] == d['ML_energy'])
    c[1] = (np.cos(d['zenith']) >= 0.8) * (d['zenith']!=0)
    c[2] = inPoly(d['ML_x'], d['ML_y'], 0, config='IT73')
    c[3] = (np.logical_not(d['LoudestOnEdge']))
    c[4] = (d['Q1'] >= 6)

    tempCut = np.array([True for i in c[0]])
    ntotal = float(tempCut.sum())
    for k in c.keys():
        tempCut *= c[k]
        print k
        print '%.4f' % (c[k].sum()/ntotal)
        print '%.4f' % (tempCut.sum()/ntotal)
