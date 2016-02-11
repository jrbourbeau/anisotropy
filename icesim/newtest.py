#!/usr/bin/env python

import numpy as np
import healpy as hp
import argparse
from analysis import getSimFiles, getSplineEnergies
from saveMedian import getWeights

def buildHists(config, outFile):

    recoEbins = np.array([4, 4.25, 4.5, 4.75, 5, 5.25, 5.5, 6, 6.5, 100])
    trueEbins = np.arange(2.75, 9.01, 0.05)
    nr, nt = len(recoEbins)-1, len(trueEbins)-1
    d = {}
    d['energy'] = np.zeros((nr, nt))
    zbins = np.linspace(0, 1, 151)
    nz = len(zbins)-1
    d['zenith'] = np.zeros((nr, nz))

    simFiles = getSimFiles(config)
    for i, simFile in enumerate(simFiles):

        # Load simulation file
        print 'Working on %s' % simFile
        s = np.load(simFile)
        s = s.item()

        # Calculate splined energies
        if 'reco2' in s.keys():
            theta, phi = hp.pix2ang(1024, s['reco2'].astype('int'))
            s['dstZenith']  = theta
            s['dstAzimuth'] = phi
        x = np.cos(s['dstZenith'])
        z = np.cos(s['zen'])
        with np.errstate(divide='ignore'):
            y = np.log10(s['nchannel'])
            e = np.log10(s['energy'])
        esp = getSplineEnergies(config, x, y)

        # Bin in histogram
        badRunFile = '/data/user/fmcnally/anisotropy/sim/badFiles.txt'
        s['type'] = s['type'].astype('int')
        weight = getWeights(s, config, badRunFile)
        d['energy'] += np.histogramdd(np.transpose([esp,e]), 
                bins=[recoEbins,trueEbins],
                weights=weight)[0]
        d['zenith'] += np.histogramdd(np.transpose([esp,z]),
                bins=[recoEbins,zbins],
                weights=weight)[0]

    if outFile != False:
        np.save(outFile, d)



if __name__ == "__main__":

    p = argparse.ArgumentParser()
    p.add_argument('-c', '--config', dest='config',
            default='IC86',
            help='Detector configuration')
    p.add_argument('-o', '--outFile', dest='outFile',
            default=False, action='store_true',
            help='Write to file')
    args = p.parse_args()

    if args.outFile:
        args.outFile = '/home/fmcnally/anisotropy/icesim/%s_hists_test.npy' % \
                args.config

    buildHists(args.config, args.outFile)



