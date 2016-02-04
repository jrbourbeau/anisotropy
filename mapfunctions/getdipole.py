#!/usr/bin/env python

import healpy as hp
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import quad
import os, glob, ROOT, time
import matplotlib.pyplot as plt
from matplotlib import rc
from mapFunctions import getFitParams, outputFit, multifit
import argparse
from energyCuts import getEnergyMaps

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

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
    parser.add_argument('--lmax', dest='lmax', type=int,
            help='Highest l value to be fit. (l=1: dipole, l=2: quadrupole, etc.)')
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
            default='/home/jbourbeau/public_html/figures/dipole/',
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
            default=False, help='Print params')
    parser.add_argument('--usefit', action='store_true', dest='usefit',
            default=False, help='Use spherical harmonic fits instead of raw data')
    parser.add_argument('--plotname', dest='plotname',
            help='Output plot filename')


    args = parser.parse_args()
    opts = vars(args).copy()
   
    alpha = 1/20.
    config = 'IC'
    maps = getEnergyMaps(config)
    dipoles = []
    bins = []
    energy_bin_median = [4.12, 4.38, 4.58, 4.85, 5.12, 5.38, 5.77, 6.13, 6.73]
    energy_err_lower = [0.62, 0.65, 0.68, 0.73, 0.74, 0.75, 0.6, 0.52, 0.46]
    energy_err_upper = [0.5, 0.54, 0.55, 0.64, 0.72, 0.78, 0.83, 0.63, 0.58]
    '''for i in range(len(maps)):
        if len(maps[i]) == 1:
            data, bg, local = hp.read_map(maps[i][0], range(3),verbose=False)
        else:
            # Read in (multiple) input files
            data, bg, local = np.sum([hp.read_map(f, range(3), verbose=False)\
                    for f in maps[i]], axis=0)
        p = multifit(args.lmax, data, bg, alpha, **opts)
        dp = np.sqrt(3.0/(4.0*np.pi))*np.sqrt(p['Y(1,-1)']**2.0 + p['Y(1,1)']**2.0)
        dipoles.append(dp)
        bins.append(i+1)
        print('dipole amplitude = {}'.format(dp))'''
   

    # Right ascension projection of relative intensity map
    data, bg, local = hp.read_map(maps[0][0], range(3),verbose=False)
    nside = hp.get_nside(data)
    dec, ra = hp.pix2ang(nside,[i for i in range(hp.nside2npix(nside))])
    rad2deg = 180./np.pi
    deg2rad = np.pi/180.
    dec = np.array([rad2deg*i for i in dec])
    decmin = 115.
    decmax = 180.
    deccut = (dec >= decmin)*(dec <= decmax)
    ra = np.array([rad2deg*i for i in ra])
    dipole_avg = []
    fig = plt.figure(0)
    num_bands = 24.
    for j in range(len(maps)):
        dI = []
        angle = []
        if len(maps[j]) == 1:
            data, bg, local = hp.read_map(maps[j][0], range(3),verbose=False)
        else:
            # Read in (multiple) input files
            data, bg, local = np.sum([hp.read_map(f, range(3), verbose=False)\
                    for f in maps[j]], axis=0)
        for i in range(int(num_bands)):
            if args.usefit:
                relint = multifit(args.lmax, data, bg, alpha, **opts)
                #print('relint = {}'.format(relint))
            else:
                with np.errstate(invalid='ignore'):
                    relint = (data - bg) / bg
                relint[np.isnan(relint)] = 0.
            
            # Right ascension cut 
            ramin = (360./num_bands)*i
            ramax = (360./num_bands)*(i+1)
            racut = (ra >= ramin)*(ra <= ramax)

            #hp.write_map("racut.fits", [data, bg, local])
            n_pix = (racut*deccut).sum()
            #print('n_pix = {}'.format(n_pix))
            relint[np.logical_not(racut)] = 0.0 
            relint[np.logical_not(deccut)] = 0.0 
            dI.append(relint.sum()/n_pix)
            angle.append(ramin+((360./num_bands)/2.))
        
        anglemin = min(angle)
        anglemax = max(angle)

        def relint_avg(x):
            for i in range(int(num_bands)):
                if rad2deg*x >= (360./num_bands)*i and rad2deg*x <= (360./num_bands)*(i+1):
                    return dI[i]
        
        relint_interp = interp1d([deg2rad*x for x in angle], dI, kind='cubic')
        def integrand1(x):
            return np.cos(x)*relint_interp(x)/np.pi 
        def integrand2(x):
            return np.sin(x)*relint_interp(x)/np.pi
        '''def integrand1(x):
            return np.cos(x)*relint_avg(x)/np.pi 
        def integrand2(x):
            return np.sin(x)*relint_avg(x)/np.pi'''
        eq1 = quad(integrand1,deg2rad*anglemin,deg2rad*anglemax)[0]
        eq2 = quad(integrand2,deg2rad*anglemin,deg2rad*anglemax)[0]
        a0 = np.arctan(eq2/eq1)
        A1 = eq1/np.cos(a0)
        #print(a0)
        #print(A1)
        if A1 < 0.:
            a0 += np.pi
            A1 = eq1/np.cos(a0)
            dipole_avg.append(A1)
        else:
            dipole_avg.append(A1)
        
        if j==0:
            xnew = np.arange(anglemin, anglemax, 1)
            #plt.plot(angle,dI,marker='.',markersize=5,linestyle='None')
            plt.plot(xnew,[relint_interp(i*deg2rad) for i in xnew],marker='None',lw=1,label='IC energy bin {}'.format(j+1))
            plt.plot(xnew,[relint_avg(i*deg2rad) for i in xnew],marker='None',lw=1)
            plt.hlines(0,0.0,360.0,linestyle='-.')
            plt.legend()
            ax = fig.axes[0]
            ax.axis('on')
            tPars = {'fontsize':16}
            ax.set_xlim(0.0,360.)
            ax.set_ylim(-0.0015,0.0010)
            ax.set_xlabel(r'Right Ascension', **tPars)
            ax.set_ylabel(r'Relative Intensity',**tPars)
            ax = plt.gca()
            ax.invert_xaxis()

    #print('dipole_avg = {}'.format(dipole_avg))
    
    if args.output:
        if args.plotname:
            outFile  = args.outDir + args.plotname
        else:
            if args.lmax == 1:
                name = 'relint_dipolefit'
            elif args.lmax == 2:
                name = 'relint_quadrupolefit'
            elif args.lmax == 3:
                name = 'relint_octupolefit'
            else:
                name = 'relint_data'
            outFile  = args.outDir + name
        plt.savefig(outFile+'.'+args.ext, dpi=300, bbox_inches='tight')
    
    fig = plt.figure(1)
    plt.errorbar(energy_bin_median,[1.*i for i in dipole_avg],xerr=[energy_err_lower, energy_err_upper], fmt='o')
    #plt.errorbar(energy_bin_median,dipoles,xerr=[energy_err_upper, energy_err_lower], fmt='o')
    ax = fig.axes[0]
    ax.axis('on')
    tPars = {'fontsize':16}
    ax.set_ylim(0.0,0.0016)
    ax.set_xlabel(r'$\log_{10}(E/\mathrm{GeV})', **tPars)
    ax.set_ylabel(r'Dipole Amplitude',**tPars)
    if args.output:
        if args.plotname:
            outFile  = args.outDir + args.plotname
        else:
            if args.lmax == 1:
                name = 'dipole_dipolefit'
            elif args.lmax == 2:
                name = 'dipole_quadrupolefit'
            elif args.lmax == 3:
                name = 'dipole_octupolefit'
            else:
                name = 'dipole_data'
            outFile  = args.outDir + name
        plt.savefig(outFile+'.'+args.ext, dpi=300, bbox_inches='tight')
    
    '''fig = plt.figure(2)
    plt.errorbar(energy_bin_median,dipoles,xerr=[energy_err_lower, energy_err_upper], fmt='o')
    #plt.errorbar(energy_bin_median,dipoles,xerr=[energy_err_upper, energy_err_lower], fmt='o')
    ax = fig.axes[0]
    ax.axis('on')
    tPars = {'fontsize':16}
    ax.set_ylim(0.0,0.0016)
    ax.set_xlabel(r'$\log_{10}(E/\mathrm{GeV})', **tPars)
    ax.set_ylabel(r'Dipole Amplitude',**tPars)
    if args.output:
        if args.plotname:
            outFile  = args.outDir + args.plotname
        else:
            if args.lmax == 1:
                name = 'dipolefit'
            if args.lmax == 2:
                name = 'quadrupolefit'
            if args.lmax == 3:
                name = 'octupolefit'
            outFile  = args.outDir + name
        plt.savefig(outFile+'.'+args.ext, dpi=300, bbox_inches='tight')'''
    #plt.show()








    
