#!/usr/bin/env python

import healpy as hp
import numpy as np
import os, glob, ROOT, time
from numpy import sqrt, pi
from numpy.polynomial import Legendre
import matplotlib.pyplot as plt
from mapFunctions import getFitParams, outputFit, multifit
import argparse

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(
            description='description: Plotting program for healpix maps.')
    
    parser.add_argument('-f','--filepath', help='Path to healpix fit file to be analyzed')
    parser.add_argument('-m', '--min', dest='min', help='Plot minimum value')
    parser.add_argument('-M', '--max', dest='max', help='Plot maximum value')
    parser.add_argument('-d', '--decmin', dest='decmin', type=float,
            default=-90., help='Minimum declination value (90->-90)')
    parser.add_argument('-D', '--decmax', dest='decmax', type=float,
            default=90., help='Maximum declination value (90->-90)')
    parser.add_argument('-r', '--ramin', dest='ramin', type=float,
            help='Minimum RA value')
    parser.add_argument('-R', '--ramax', dest='ramax', type=float,
            help='Maximum RA value')
    parser.add_argument('--mask', dest='mask', default=False, 
            action='store_true', help='Intelligent masking')
    parser.add_argument('-b', '--batchmode', action='store_true', dest='batch',
            default=False, help='Execute without interaction')
    parser.add_argument('-o', '--output', dest='output', action='store_true',
            default=False, help='Output image file')
    parser.add_argument('-n', '--mapName', dest='mapName',
            help='Map type desired (signal, relint, relint_err, data, bg)')
    parser.add_argument('-s', '--scale', dest='scale', type=float,
            help='Scale the map after input')
    parser.add_argument('-S', '--smooth', dest='smooth', type=float, default=0,
            help='Desired smoothing radius (in degrees)')
    parser.add_argument('--stype', dest='stype', default='tophat',
            help='Option for smoothing type [tophat|gauss]')
    parser.add_argument('--swindow', dest='swindow', type=float, default=3,
            help='Option for smoothing window')
    parser.add_argument('--multi', dest='multi', type=int,
            help='Use Multipole subtraction')
    parser.add_argument('--fix_multi', dest='fix_multi',
            default=False, action='store_true',
            help='Fix multipole subtraction to values from cumulative map')
    parser.add_argument('-x', '--threshold', dest='threshold', type=float,
            help='Threshold value for plotting data')
    parser.add_argument('-c', '--coords', dest='coords', default='C',
            help='C=equatorial, G=galactic, E=ecliptic')
    parser.add_argument('--gplane', dest='gplane',
            default=False, action='store_true',
            help='Show the galactic plane')
    parser.add_argument('--half', action="store_true", dest='half',
            default=False, help='Show only bottom half of the sky')
    parser.add_argument('--title', action="store_true", dest="title",
            default=False, help='Show the title on the plot')
    parser.add_argument('--outDir', dest='outDir',
            default='/home/jbourbeau/public_html/figures/',
            help='Option for changing output directory')
    parser.add_argument('--prelim', action='store_true', dest='prelim',
            default=False, help='Indicate plot is preliminary')
    parser.add_argument('--llabel', dest='llabel',
            default=False, help='Optional left label overlay on map')
    parser.add_argument('--rlabel', dest='rlabel',
            default=False, help='Optional right label overlay on map')
    parser.add_argument('--polar', dest='polar',
            default=False, action='store_true',
            help='Polar gnomonic view of map')
    parser.add_argument('--customOut', dest='customOut',
            help='Option for custom output file name')
    parser.add_argument('--ext', dest='ext',
            default='png', help='Output file extension')
    parser.add_argument('-v', '--verbose', action='store_true', dest='verbose',
            default=False, help='Verbose output')
    parser.add_argument('--params', action='store_true', dest='params',
            default=True, help='Print params')


    args = parser.parse_args()
    opts = vars(args).copy()
   
    alpha = 1/20.
    # Read in (multiple) input files
    #data, bg, local = np.sum([hp.read_map(f, range(3), verbose=False)
    #        for f in inFiles], axis=0)
    data, bg, local = hp.read_map(args.filepath, range(3))
    print('multi = {}'.format(args.multi))
    print('params = {}'.format(args.params))
    p = multifit(args.multi, data, bg, alpha, **opts)
    for keys in p:
        print('p[{}] = {}'.format(keys,p[keys]))
    dp = np.sqrt(3.0/(4.0*np.pi))*np.sqrt(p['Y(1,-1)']**2.0 + p['Y(1,1)']**2.0)
    print('dipole amplitude = {}'.format(dp))








    
