#!/usr/bin/env python

############################################################################
# Makes the bins for the LLH tables
############################################################################

import myGlobals as my
import numpy as np
import argparse
import os


if __name__ == "__main__":

    # Global variables setup for path names
    my.setupShowerLLH(verbose=False)
    resourcedir = my.llh_resource

    p = argparse.ArgumentParser(
        description='Makes binning for ShowerLLH studies')
    p.add_argument('-o', '--outFile', dest='outFile',
                   default=resourcedir + '/ShowerLLH_bins.npy',
                   help='Destination file for bin information')
    args = p.parse_args()

    degree = np.pi / 180.
    binVals = {}

    for bintype in ['standard', 'nozenith', 'logdist']:

        bins = {}
        # Energy bins in Log10(Energy/GeV)
        bins['E'] = np.arange(4, 9.501, 0.05)
        # Distance bins in meters
        bins['D'] = np.append(np.arange(0, 600, 10), np.arange(600, 1051, 50))
        bins['D'] = np.append(bins['D'], np.inf)
        # Zenith Bins in radians (made with equal solid angle bins)
        bins['Z'] = np.linspace(1, np.cos(40 * degree), 4)
        bins['Z'] = np.append(np.arccos(bins['Z']), np.pi / 2)
        # Snow bins in meters
        bins['S'] = np.array([-1, .001, .5, .85])
        bins['S'] = np.append(bins['S'], np.inf)
        # Charge bins in VEM
        bins['C'] = 10**np.linspace(-3, 4.5, 46)[10:38]
        bins['C'] = np.append(0, bins['C'])
        bins['C'] = np.append(bins['C'], np.inf)

        # Custom treatment for alternative bintypes
        order = 'EZSDC'
        if bintype == 'nozenith':
            order = order.replace('Z', '')
        if bintype == 'logdist':
            order = order.replace('Z', '')
            bins['D'] = 10**np.linspace(np.log10(5), np.log10(1050), 149)
            bins['D'] = np.append(0, bins['D'])
            bins['D'] = np.append(bins['D'], np.inf)

        # Bin order determined by order argument
        binVals[bintype] = {}
        for i, param in enumerate(order):
            binVals[bintype][i] = [param, bins[param]]

    #print('binVals = {}'.format(binVals))
    # Save to outFile
    np.save(args.outFile, binVals)
