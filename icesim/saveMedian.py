#!/usr/bin/env python

import numpy as np
import healpy as hp
import glob, argparse
import myGlobals as my
import simFunctions_IC as simFunctions
from analysis import getSimBase


def histPercentile(h, p, binMids):
    tot = h.sum()
    if tot == 0:
        return 0
    i, cts = 0,0
    while cts < (p*tot)/100.:
        cts += h[i]
        i += 1
    return binMids[i-1]

def histVar(h, binMids):
    if h.sum() == 0:
        return 0
    vals = np.asarray(binMids)
    ave = np.average(vals, weights=h)
    var = np.average((vals-ave)**2, weights=h)
    return var


def histMedian(h, bins):

    # Assume you only want to operate on the last dimension
    nx, ny, nz = h.shape
    median, sigL, sigR, var = np.zeros((4,nx,ny))
    binMids = (bins[:-1] + bins[1:]) / 2.
    for i in range(nx):
        for j in range(ny):
            median[i][j] = histPercentile(h[i][j], 50, binMids)
            sigL[i][j]   = histPercentile(h[i][j], 16, binMids)
            sigR[i][j]   = histPercentile(h[i][j], 84, binMids)
            var[i][j]    = histVar(h[i][j], binMids)

    return median, sigL, sigR, var


def getWeights(s, config, badRunFile):

    # IC86 simulation files are already weighted
    if config in ['IC86','IC86-II','IC86-III']:
        return None

    sim = int(simFunctions.cfg2sim(config))
    from icecube.weighting.fluxes import Hoerandel5
    from icecube.weighting.weighting import from_simprod

    # Load bad simulation files
    with open(badRunFile, 'r') as f:
        badFiles = f.readlines()
        badFiles = [l.strip() for l in badFiles if str(sim) in l.split('/')]
        nbad = len(badFiles)

    # Make generator
    nfiles, generator = from_simprod(sim)
    nfiles -= nbad
    generator *= nfiles
    flux = Hoerandel5()

    print 'Calculating weights...'
    weights = flux(s['energy'], s['type']) / \
            generator(s['energy'], s['type'])

    # Eliminate CNO and MgSiAl components that weighting can't deal with
    weights[weights!=weights] = 0
    weights[weights==np.inf]  = 0

    return weights


if __name__ == "__main__":

    # Setup global path names
    my.setupAnisotropy(verbose=False)

    p = argparse.ArgumentParser(
            description='Saves median information for energy reconstruction')
    p.add_argument('-c', '--config', dest='config',
            help='Detector configuration')
    p.add_argument('--nx', dest='nx',
            default=60,
            help='Number of bins in cos(zenith)')
    p.add_argument('--ny', dest='ny',
            default=40,
            help='Number of bins in log10(Nchannel)')
    p.add_argument('--estep', dest='estep',
            default=.05,
            help='Energy step size in log10(E/GeV)')
    p.add_argument('-o', '--out', dest='out',
            help='Option for outfile name')
    args = p.parse_args()

    if not args.out:
        simBase = getSimBase(args.config)
        outFile = '%s_median.npy' % simBase
    else:
        outFile = '%s/%s.npy' % (my.ani_sim, args.out)

    # Create histogram
    xbins = np.linspace(0.3, 1, args.nx+1)
    ybins = np.linspace(1, 3.2, args.ny+1)
    ebins = np.arange(2.75, 9+args.estep/2., args.estep)
    h = np.zeros((args.nx, args.ny, len(ebins)-1))
    #w = np.zeros((args.nx, args.ny, len(ebins)-1))

    # Get list of simulation files
    simBase = getSimBase(args.config)
    simList = glob.glob('%s.npy' % simBase)
    if len(simList) == 0:
        simList = glob.glob('%s_Part*.npy' % simBase)
    simList.sort()

    # Get bad run file
    badRunFile = '%s/badFiles.txt' % my.ani_sim

    for simFile in simList:

        # Load data
        print 'Loading %s...' % simFile
        d = np.load(simFile)
        d = d.item()

        if 'reco2' in d.keys():
            theta, phi = hp.pix2ang(1024, d['reco2'].astype('int'))
            d['dstZenith'] = theta

        x = np.cos(d['dstZenith'])
        y = np.log10(d['nchannel'])
        e = np.log10(d['energy'])
        #nfiles = int(sum(d['nFiles']))
        d['type'] = d['type'].astype('int')
        weight = getWeights(d, args.config, badRunFile)

        # Fill weighted histogram
        print 'Populating histogram...'
        #h += np.histogramdd(np.transpose([x,y,e]), bins=[xbins,ybins,ebins])[0]
        #w += np.histogramdd(np.transpose([x,y,e]), bins=[xbins,ybins,ebins],
        #        weights=weight)[0]
        h += np.histogramdd(np.transpose([x,y,e]), bins=[xbins,ybins,ebins],
                weights=weight)[0]

    # Setup data storage arrays
    t = {}
    histValues = ['medians','sigma_min','sigma_max','var']
    #histValues += ['wmedians','wsigma_min','wsigma_max','wvar']
    for value in histValues:
        t[value] = np.zeros((args.nx, args.ny))

    print 'Calculating mean and median values and errors...'
    # Unweighted cases
    median, sigL, sigR, var = histMedian(h, ebins)
    t['medians']   = median
    t['sigma_min'] = sigL
    t['sigma_max'] = sigR
    t['var']       = var
    # Weighted cases
    #median, sigL, sigR, var = histMedian(w, ebins)
    #t['wmedians']   = median
    #t['wsigma_min'] = sigL
    #t['wsigma_max'] = sigR
    #t['wvar']       = var

    # Additional information
    t['xbins'] = xbins
    t['ybins'] = ybins
    t['ebins'] = ebins

    print 'Writing to', outFile, '...'
    np.save(outFile, t)


