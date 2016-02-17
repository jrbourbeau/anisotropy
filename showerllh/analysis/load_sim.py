#!/usr/bin/env python

################################################################################
## Load essential information for simulation analysis                         ##
################################################################################


import numpy as np
import time, glob, re

from llhtools import inPoly
import myGlobals as my


def load_sim(config='IT73', bintype='logdist'):

    # Load simulation information
    my.setupShowerLLH(verbose=False)
    infile = '%s/%s_sim/SimPlot_%s.npy' % (my.llh_data, config, bintype)

    s = {}

    t0 = time.time()
    print 'Working on %s...' % infile
    s = np.load(infile)
    s = s.item()
    print 'Time taken to load:', time.time()-t0

    # Set array types
    bools = []
    bools += ['LoudestOnEdge']
    for cut in bools:
        s[cut] = s[cut].astype('bool')

    ## Quality Cuts ##
    c = {}

    # Adapted versions of Bakhtiyar's cuts
    c['llh1']  = np.logical_not(np.isnan(s['ML_energy']))
    c['llh2']  = (np.cos(np.pi - s['zenith']) >= 0.8)
    c['llh2a'] = (180./np.pi * (np.pi - s['zenith']) >= 40)
    c['llh3']  = inPoly(s['ML_x'], s['ML_y'], 0)
    #c['llh3a'] = inPoly(s['ML_x'], s['ML_y'], -50)
    #c['llh3t'] = inPoly(s['ML_x'], s['ML_y'], -90)
    c['llh4']  = np.logical_not(s['LoudestOnEdge'])
    c['llh5']  = (s['Q1'] >= 6)
    c['llh5a'] = np.sum([s['Q'+str(i)]>=6 for i in range(1,5)], axis=0)
    c['llh5a'] = c['llh5a'].astype('bool')

    # Final versions of cuts for testing
    c['llh']  = c['llh1'] * c['llh2'] * c['llh4']
    #c['llht'] = c['llh'] * c['llh3t'] * c['llh5']
    c['llha'] = c['llh'] * c['llh3'] * c['llh5a']
    c['llh']  = c['llh'] * c['llh3'] * c['llh5']
    s['cuts'] = c

    ## Load MC information ##
    print 'Loading MC information...'
    s['MC'] = np.load('%s/%s_sim/SimPlot_MC.npy' % (my.llh_data, config))
    s['MC'] = s['MC'].item()

    compList = s['MC'].keys()
    #if bintype == 'logdist':
    #    compList = ['P','Fe']
    s['MC']['joint'] = {}
    shape = s['MC']['P']['low'].shape
    for key in ['low', 'mid', 'high']:
        s['MC']['joint'][key] = np.zeros(shape, dtype=int)
        s['MC']['joint'][key] = np.sum([s['MC'][comp][key] 
                for comp in compList], axis=0)

    return s


if __name__ == "__main__":

    s = load_sim()




