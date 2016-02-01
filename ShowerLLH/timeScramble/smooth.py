#!/usr/bin/env python

import numpy as np
import healpy as hp
import time, os, glob
from timeScramble import LMSignificance

# Takes in dictionary with keys 'bg' and 'data'
def smooth(orig_maps, smooth_radius):

    # Create mask
    nPix = len(orig_maps['bg'])
    nSide = hp.npix2nside(nPix)
    theta, phi = hp.pix2ang(nSide, range(nPix))
    theta = np.pi - theta
    mask = np.cos(theta)>=0.8
    masked_pix = np.logical_not(mask)

    maps = {}
    smooth_radius *= np.pi/180.
    for key in orig_maps.keys():
        maps[key] = np.zeros(nPix)

    # Smooth maps
    for i in range(nPix):
        vec = hp.pix2vec(nSide, i)
        neighbors = hp.query_disc(nSide, vec, smooth_radius)
        for key in orig_maps.keys():
            maps[key][i] += orig_maps[key][neighbors].sum()

    # Make significance and relative intensity maps
    alpha = 1/20.
    Non, Noff = maps['data'], maps['bg']
    maps['signal'] = LMSignificance(Non, Noff, alpha)
    maps['relint'] = (Non-Noff) / Noff
    maps['relint_err'] = (Non/Noff) * np.sqrt(1/Non + alpha/Noff)

    # Catch NaN's
    for key in ['relint', 'relint_err']:
        nancatch = (maps[key]!=maps[key])
        maps[key][nancatch] = 0

    # Mask maps
    for key in maps.keys():
        maps[key][masked_pix] = hp.UNSEEN

    return maps


def makeSmooth(keys=False, all_any='all'):

    prefix = '/net/user/fmcnally/ShowerLLH/maps'
    #inPrefix  = prefix + '/merged'
    inPrefix  = prefix + '/raw'
    outPrefix = prefix + '/smoothed'
    fList = glob.glob(inPrefix + '/*data.fits')
    fList.sort()

    # Calculate comprehensive list of all unique parameters
    paramList = []
    testList = [os.path.basename(f) for f in fList]
    for f in testList:
        f_params = f.split('_')[:-1]     # get rid of bg/data
        paramList.append(f_params)

    # Option for limiting to specific set by keyword ('NStations','nocuts',etc)
    if keys and all_any=='all':
        paramList = [p for p in paramList if all(key in p for key in keys)]
    if keys and all_any=='any':
        paramList = [p for p in paramList if any(key in p for key in keys)]

    for params in paramList:
        print 'Working on', params
        orig_maps = {}
        start = '_'.join(params)
        for key in ['bg', 'data']:
            inFile = '%s/%s_%s.fits' % (inPrefix, start, key)
            orig_maps[key] = hp.read_map(inFile)

        #deg_range = range(6, 61)       # Unnecessary in general
        deg_range = [6, 20, 40]
        for deg in deg_range:
            print deg, 'degree smoothing...'
            maps = smooth(orig_maps, deg)
            for key in ['relint', 'relint_err', 'signal']:
                out = '%s/%s_%sdeg_%s.fits' % (outPrefix, start, deg, key)
                hp.write_map(out, maps[key])


if __name__ == "__main__":

    prefix = '/net/user/fmcnally/ShowerLLH/maps/raw/'
    mList = glob.glob(prefix + 'IT??_*')
    mList = [os.path.basename(m) for m in mList]
    mList = [m.split('_')[1] for m in mList]
    mList = list(set(mList))
    mList.sort()

    for month in mList:

        keys = [month, 'nocuts']
        makeSmooth(keys=keys)







