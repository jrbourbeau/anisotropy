#!/usr/bin/env python

################################################################################
## Load essential information for IT73 simulation analysis                    ##
################################################################################


import numpy as np
import time, sys

from llhtools import inPoly
import myGlobals as my


def load_sim(config, bintype='standard'):

    # Load simulation information
    my.setupShowerLLH(verbose=False)
    prefix = '%s/%s_sim' % (my.llh_data, config)
    infile = '%s/SimPlot_%s.npy' % (prefix, bintype)
    t0 = time.time()
    s = np.load(infile)
    s = s.item()
    print 'Time taken to load:', time.time()-t0

    # Set array types
    bools = []
    #bools  = ['rrc', 'lowE', 'medE', 'highE', 'LoudestOnEdge']
    #bools += ['P', 'He', 'O', 'Fe', 'p', 'h', 'o', 'f']
    bools += ['P', 'Fe', 'p', 'f']
    for cut in bools:
        s[cut] = s[cut].astype('bool')

    ## Quality Cuts ##
    s['cuts'] = {}

    # Zenith cuts
    cosz = np.cos(s['MC_zenith'])
    s['cuts']['z40'] = (cosz >= 0.95)
    s['cuts']['z41'] = (cosz >= 0.90) * (cosz < 0.95)
    s['cuts']['z42'] = (cosz >= 0.85) * (cosz < 0.90)
    s['cuts']['z43'] = (cosz >= 0.80) * (cosz < 0.85)
    s['cuts']['z30'] = (cosz >= 0.93)
    s['cuts']['z31'] = (cosz >= 0.87) * (cosz < 0.93)
    s['cuts']['z32'] = (cosz >= 0.80) * (cosz < 0.87)

    # Bakhtiyar's cuts
    #s['cuts']['bakh1'] = (s['lap_fit_status'] == 0)
    s['cuts']['bakh2'] = (np.cos(s['lap_zenith']) >= 0.8)
    s['cuts']['bakh3'] = inPoly(s['lap_x'], s['lap_y'], -90)
    s['cuts']['bakh4'] = np.logical_not(s['LoudestOnEdge'])
    s['cuts']['bakh5'] = (s['Q1'] >= 6)
    s['cuts']['lapE']  = (s['lap_e_proton']==s['lap_e_proton'])
    s['cuts']['lapE'] *= (s['lap_e_proton']!=0)

    # Adapted versions of Bakhtiyar's cuts
    s['cuts']['llh1'] = s['rrc']
    s['cuts']['llh2'] = (np.cos(s['ShowerPlane_zenith']) >= 0.8)
    s['cuts']['llh3'] = inPoly(s['ML_x'], s['ML_y'], -50)
    s['cuts']['llh3t']= inPoly(s['ML_x'], s['ML_y'], -90)
    s['cuts']['llh4'] = np.logical_not(s['LoudestOnEdge'])
    s['cuts']['llh5'] = (s['Q1'] >= 6)

    # Final versions of cuts for testing
    s['cuts']['llh']  = s['cuts']['llh1'] * s['cuts']['llh2']
    s['cuts']['llh'] *= s['cuts']['llh4'] * s['cuts']['llh5']
    s['cuts']['llht'] = s['cuts']['llh'] * s['cuts']['llh3t']
    s['cuts']['llh'] *= s['cuts']['llh3']

    # Bakhtiyar's final cut
    s['cuts']['bakh']  = s['cuts']['bakh2']*s['cuts']['bakh3']
    s['cuts']['bakh'] *= s['cuts']['bakh4']*s['cuts']['bakh5']*s['cuts']['lapE']

    # Select the most common outliers
    ## Note : no longer really used
    s['cuts']['out'] = (abs(np.log10(s['ML_energy'] / s['MC_energy'])) > 0.5)
    s['cuts']['in']  = np.logical_not(s['cuts']['out'])
    s['cuts']['out'] *= (np.log10(s['MC_energy']) > 6)
    s['cuts']['in']  *= (np.log10(s['MC_energy']) > 6)


    ## Load MC information ##
    # Get the counts for MC showers
    s['MC'] = np.load('%s/%s_MC.npy' % (prefix, config))
    s['MC'] = s['MC'].item()

    compList = s['MC'].keys()
    s['MC']['joint'] = {}
    shape = s['MC']['P']['low'].shape
    for erange in ['low', 'med', 'high']:
        s['MC']['joint'][erange] = zeros(shape, dtype=int)
        for i in range(3,5):
            for j in range(i):
                k = str(i) + str(j)
                s['MC']['joint'][erange+'_'+k] = zeros(shape, dtype=int)
    for comp in compList:
        for key in s['MC'][comp].keys():
            s['MC']['joint'][key] += s['MC'][comp][key]

    return s


if __name__ == "__main__":

    config = sys.argv[1]
    s = load_sim(config)




