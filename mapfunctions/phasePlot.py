#!/usr/bin/env python

import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from scipy.optimize import curve_fit
from tabulate import tabulate
import glob, argparse, os

import myGlobals as my
from useful import getMids
from mapFunctions import getMap
from proj1d import returnRI_new as returnRI
from projFunctions import *
from energyCuts import getEbins, getEnergyMaps
from anisotropy.icesim.analysis import readDist as readDist_IC
from anisotropy.topsim.analysis import readDist as readDist_IT

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)


if __name__ == "__main__":

    # Setup global path names
    my.setupAnisotropy(verbose=False)

    p = argparse.ArgumentParser(
            description='')
    p.add_argument('-f', '--files', dest='files',
            default='energy',
            choices=['energy','config','comp','solar'],
            help='Choose from presets for files to run over')
    p.add_argument('--chi2', dest='chi2',
            default=False, action='store_true',
            help='Show chi2 for a variety of "multipole" fits')
    p.add_argument('-l', '--lvalue', dest='l', type=int,
            help="Option to fix l-value for multipole fit")
    p.add_argument('-n', '--nbins', dest='nbins', 
            type=int, default=24,
            help='Number of bins for the 1D fit')
    p.add_argument('--solo', dest='solo', type=int,
            help="Option to plot an individual file's sine fits")
    p.add_argument('--offset', dest='offset', type=int,
            help='Shift phase by some offset value')
    p.add_argument('-b', '--batch', dest='batch',
            default=False, action='store_true',
            help="Execute without interaction")
    p.add_argument('-o', '--out', dest='out',
            default=None,
            help="Option to write to file")
    args = p.parse_args()

    if args.files == 'config':
        files = glob.glob('%s/IC?*_24H_sid.fits' % my.ani_maps)
        configs = np.array([os.path.basename(f).split('_')[0] \
                for f in files])
        csort = configs.argsort()
        configs, files = configs[csort], np.asarray(files)[csort]
        x = range(len(configs))
        smooth = 5
    if args.files == 'energy':
        files = getEnergyMaps('IC')
        ebins = [float(i) for i in getEbins()]
        nbins = len(ebins)-1
        x = [readDist_IC('IC', ebins[i], ebins[i+1])[0] for i in range(nbins)]
        xL = [readDist_IC('IC', ebins[i], ebins[i+1])[1] for i in range(nbins)]
        xR = [readDist_IC('IC', ebins[i], ebins[i+1])[2] for i in range(nbins)]
        #x = getMids(ebins, infvalue=100)
        smooth = 20
        # Custom treatment of IT map
        files += [['{}/IT_24H_sid_STA8.fits'.format(my.ani_maps)]]
        #x = np.append(x, np.log10(2*10**6))
        median, sigL, sigR = readDist_IT('IT', 8, 100)
        x  += [median]
        xL += [sigL]
        xR += [sigR]
    if args.files == 'comp':
        files = glob.glob('%s/IT_24H_sid_STA8_?.fits' % my.ani_maps)
        comps = [f[-6] for f in files]
        #comp2e = {'p':1.,'h':2.,'o':8.,'f':26.}
        comp2e = {'p':1.,'h':4.,'o':16.,'f':56.}
        x = [1./comp2e[comp] for comp in comps]
        smooth = 20
    if args.files == 'solar':
        files = glob.glob('%s/IC_24H_solar.fits' % my.ani_maps)
        smooth = 5
        x = ['solar']

    if args.solo == 0:
        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot(111)

    nfiles = len(files)
    lmax = args.nbins/2 - 1
    if lmax > 15:
        lmax = 15

    opts = {'mask':True, 'ramin':0., 'ramax':360., 'nbins':args.nbins}
    amp, phase, amp_err, phase_err = np.zeros((4, nfiles))
    chi2array = np.zeros((nfiles,lmax-1))

    # Calculate best harmonic function fit value based on chi2
    # Fill chi-squared array
    for i, file in enumerate(files):
        if type(file) != list:
            file = [file]
        # Get relint and relerr maps
        relint = getMap(*file, mapName='relint', **opts)
        relerr = getMap(*file, mapName='relerr', **opts)
        ra, ri, ra_err, ri_err = getRIRAProj(relint, relerr,**opts)
        for l in range(1, lmax):
            popt, perr, chi2 = getHarmonicFitParams(ra, ri, l, ri_err)
            chi2array[i][l-1] = chi2
            if i == args.solo:
                if l==1:
                    ax.errorbar(ra, ri, ri_err, fmt='.',label=r'Data')
                fullra = range(0,360,5)
                fit = [cosFit(j, *popt) for j in fullra]
                ax.plot(fullra, fit, label=r'$l_{} = {}$ Fit'.format('{max}',l))
                ax.set_xlim(0.,360.)
                ax.invert_xaxis()

    if args.chi2:
        chi2table = chi2array.tolist()
        chi2table = [[x[i]] + chi2 for i, chi2 in enumerate(chi2table)]
        print
        print tabulate(chi2table, headers=range(1,lmax), floatfmt='.2f')

    if not args.l:
        chi2array[chi2array < .7] = 100
        args.l = chi2array.sum(axis=0).argmin() + 1

    opts['lmax'] = args.l
    for i, file in enumerate(files):
        # Get relint and relerr maps
        relint = getMap(*file, mapName='relint', **opts)
        relerr = getMap(*file, mapName='relerr', **opts)
        amp[i],amp_err[i],phase[i],phase_err[i] = getProjDipole(relint, relerr, **opts)

    # Deal with negative amplitudes/phases
    c0 = amp < 0
    amp[c0] *= -1.
    phase[c0] += 180
    phase[phase > 360] -= 360
    phase[phase < 0] += 360

    if args.files == 'energy':
        idx0 = len(ebins) - 1
        x_ic, x_it = x[:idx0], x[idx0:]
        xerr_ic = [xL[:idx0], xR[:idx0]]
        xerr_it = [xL[idx0:], xR[idx0:]]
        amp_ic, phase_ic = amp[:idx0], phase[:idx0]
        amp_err_ic, phase_err_ic = amp_err[:idx0], phase_err[:idx0]
        amp_it, phase_it = amp[idx0:], phase[idx0:]
        amp_err_it, phase_err_it = amp_err[idx0:], phase_err[idx0:]

    pltParams = {'fmt':'.', 'lw':2, 'ms':14}

    if args.solo == None:
        fig = plt.figure(figsize=(17,6))
        ax = fig.add_subplot(121)
        ax.set_xlabel(r'$\mathrm{log}_{10}(E/\mathrm{GeV})$', fontsize=14)
        #ax.set_ylabel(r'$\Delta N/\langle N \rangle$', fontsize=14)
        ax.set_ylabel(r'$A_1$', fontsize=14)
        if args.files == 'config':
            ax.set_xlim(-1, x[-1]+1)
            ax.set_xticks(x)
            ax.set_xticklabels(configs)
            ax.set_xlabel('Detector Configuration')
        if args.files == 'comp':
            ax.set_xlabel('Energy per Nucleon')
            ax.set_xlim(0, 1.1)
        if args.files == 'energy':
            ax.errorbar(x_ic, amp_ic, xerr=xerr_ic, yerr=amp_err_ic, 
                    c='b', **pltParams)
            ax.errorbar(x_it, amp_it, xerr=xerr_it, yerr=amp_err_it, 
                    c='r', **pltParams)
        else:
            ax.errorbar(x, phase, yerr=phase_err, **pltParams)

        if args.offset != None:
            if args.files == 'energy':
                pcut = phase_ic > args.offset - 20
                phase_ic[pcut] -= args.offset
                phase_ic[np.logical_not(pcut)] += (360 - args.offset)
                pcut = phase_it > args.offset
                phase_it[pcut] -= args.offset
                phase_it[np.logical_not(pcut)] += (360 - args.offset)
            else:
                pcut = phase > args.offset
                phase[pcut] -= args.offset
                phase[np.logical_not(pcut)] += (360 - args.offset)

        ax = fig.add_subplot(122)
        ax.set_xlabel(r'$\mathrm{log}_{10}(E/\mathrm{GeV})$', fontsize=14)
        #ax.set_ylabel(r'Dipole phase $(\phi/\mathrm{deg})$', fontsize=14)
        #ax.set_ylabel(r'Right Ascension $[^{\circ}]$', fontsize=14)
        ax.set_ylabel(r'$\alpha_1 \ [^{\circ}]$', fontsize=14)
        if args.files == 'config':
            ax.set_xlim(-1, x[-1]+1)
            ax.set_xticks(x)
            ax.set_xticklabels(configs)
            ax.set_xlabel('Detector Configuration')
        if args.files == 'comp':
            ax.set_xlabel('Energy per Nucleon')
            ax.set_xlim(0, 1.1)
        if args.files == 'energy':
            ax.errorbar(x_ic, phase_ic, xerr=xerr_ic, yerr=phase_err_ic, 
                    c='b', **pltParams)
            ax.errorbar(x_it, phase_it, xerr=xerr_it, yerr=phase_err_it, 
                    c='r', **pltParams)
        else:
            ax.errorbar(x, phase, yerr=phase_err, **pltParams)
        ax.set_ylim(0,360)
    else:
        ax.legend(loc='lower right')

    if args.offset != None:
        ax.set_ylim(-25, 370)
        ax.set_yticks(range(0, 361, 45))
        ylabels  = [str(i) for i in range(args.offset, 360-1, 45)]
        ylabels += [str(i) for i in range(0, args.offset+1, 45)]
        ax.set_yticklabels(ylabels)

    if args.out != None:
        plt.savefig(args.out, dpi=300, bbox_inches='tight')

    if not args.batch:
        plt.show()





