#!/usr/bin/env python

import numpy as np
import healpy as hp
import sys, glob
from os.path import basename


def merger(config, params, avoids=False):

    # Collect all fits files
    prefix = '/net/user/fmcnally/ShowerLLH/maps/'
    masterList = glob.glob(prefix + 'raw/'+config+'*.fits')
    masterList.sort()

    # Refine to files containing every parameter in params
    fList = []
    for f in masterList:
        f_params = basename(f).split('_')
        if all(x in f_params for x in params):
            fList.append(f)

    # Refine to files containing no parameter in avoids
    if avoids:
        fList = [f for f in fList if not any(x in basename(f) for x in avoids)]

    # Merge files
    if len(fList) == 0:
        print 'No files found'
        return
    for i in range(len(fList)):
        temp = hp.read_map(fList[i])
        if i == 0:
            combined_map = np.zeros(temp.shape)
        combined_map += temp

    # Write to file
    outBase = '_'.join(params)
    outFile = prefix + 'merged/'+config+'_'+outBase
    hp.write_map(outFile, combined_map)
    

if __name__ == "__main__":

    config = 'IT'

    prefix = '/net/user/fmcnally/ShowerLLH/maps/raw/'
    fList = glob.glob(prefix + '*.fits')
    fList.sort()

    # Calculate comprehensive list of all unique parameters
    paramList = []
    testList = [basename(f) for f in fList]
    for f in testList:
        f_params = f.split('_')[2:]     # first 2 params always config & date
        if f_params not in paramList:
            paramList.append(f_params)

    # Option for limiting to specific set by keyword ('sta','GeV','nocuts',etc)
    keys = ['NotSTA8']
    avoids = False
    if keys:
        paramList = [p for p in paramList if all(key in p for key in keys)]

    for params in paramList:
        print 'Reading files with', config, params
        merger(config, params, avoids=avoids)


