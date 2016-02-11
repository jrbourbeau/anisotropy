#!/usr/bin/env python

import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from tabulate import tabulate
import glob, argparse, os

import myGlobals as my
from useful import getMids
from mapFunctions import getMap
from proj1d import returnRI_new as returnRI
from energyCuts import getEbins, getEnergyMaps
from anisotropy.icesim.analysis import readDist as readDist_IC
from anisotropy.topsim.analysis import readDist as readDist_IT


def sinefit(x, *p):
    f0 = 2*np.pi / 360
    #return sum([p[3*i] * np.sin(f0*(i+1) * (x + p[3*i+1])) + p[3*i+2]
    #        for i in range(len(p)/3)])
    return sum([p[2*i] * np.sin(f0*(i+1) * (x + p[2*i+1])) 
            for i in range(len(p)/2)])

def cosinefit(x, *p):
    f0 = 2*np.pi / 360
    return sum([p[2*i+1] * np.cos(f0*(i+1) * (x - p[2*i+2]))
            for i in range(len(p)/2)]) + p[0]
    #return sum([p[2*i] * np.cos(f0*i * (x - p[2*i+1]))
    #        for i in range(len(p)/2)])

def multipoleFit(x, y, l, sigmay, fittype='cos'):
    # Guess at best fit parameters
    amplitude = (3./np.sqrt(2)) * np.std(y)
    phase     = 0
    p0 = [amplitude, phase]*l
    # Reduce amplitude as we go to higher l values
    for i in range(0,len(p0)/2):
        p0[2*i] *= 2.**(-i)
    if fittype == 'cos':
        p0 = [0] + p0
        #p0 = [0, 0] + p0
        fitfunc = cosinefit
    else:
        fitfunc = sinefit
    # Do best fit
    popt, pcov = curve_fit(fitfunc, x, y, p0, sigma=sigmay)
    fitVals = fitfunc(x, *popt)
    ndof  = len(popt)
    chi2 = (1. / (len(y)-ndof)) * sum((y - fitVals)**2 / sigmay**2)
    perr = np.sqrt(np.diag(pcov))
    #perr = perr_scaled / chi2
    return popt, perr, chi2


if __name__ == "__main__":

    # Setup global path names
    my.setupAnisotropy(verbose=False)

    p = argparse.ArgumentParser(
            description='')
    #p.add_argument('-n', '--nbins', dest='nbins', type=int,
    #        default=72,
    #        help='Number of bins to split RA into')
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
        configs = np.array([os.path.basename(f).split('_')[0] for f in files])
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
        files += [['%s/IT_24H_sid_STA8.fits' % my.ani_maps]]
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
    #nbins = 360 / smooth
    lmax = args.nbins/2 - 1
    if lmax > 15:
        lmax = 15

    opts = {'mask':True, 'ramin':0., 'ramax':360., 'nbins':args.nbins}
    #opts = {'decmax':-35.,'decmin':-88.,'ramin':0.,'ramax':360.,'nbins':nbins}
    amp, phase, amp_err, phase_err = np.zeros((4, nfiles))
    chi2array = np.zeros((nfiles,lmax-1))

    def getRelInt(file):
        if type(file) != list:
            file = [file]
        #bgmap = getMap(*file, mapName='bg', **opts)
        #datamap = getMap(*file, mapName='data', **opts)
        #bgmap[bgmap==hp.UNSEEN] = 0
        #datamap[datamap==hp.UNSEEN] = 0
        #ra, ri, sigmax, sigmay = returnRI(bgmap, datamap, **opts)
        relint = getMap(*file, mapName='relint', **opts)
        relerr = getMap(*file, mapName='relerr', **opts)
        ra, ri, sigmax, sigmay = returnRI(relint, relerr, **opts)
        return ra, ri, sigmax, sigmay

    # Calculate best multipole-fit value based on chi2
    for i, file in enumerate(files):
        ra, ri, sigmax, sigmay = getRelInt(file)
        for l in range(1, lmax):
            popt, perr, chi2 = multipoleFit(ra, ri, l, sigmay)
            chi2array[i][l-1] = chi2
            if i == args.solo:
                #fig, ax = plt.subplots()
                ax.errorbar(ra, ri, sigmay, fmt='.')
                fullra = range(360)
                #fit = [sinefit(j, *popt) for j in fullra]
                fit = [cosinefit(j, *popt) for j in fullra]
                ax.plot(fullra, fit, label=l)

    if args.chi2:
        chi2table = chi2array.tolist()
        chi2table = [[x[i]] + chi2 for i, chi2 in enumerate(chi2table)]
        print
        print tabulate(chi2table, headers=range(1,lmax), floatfmt='.2f')

    if not args.l:
        chi2array[chi2array < .7] = 100
        args.l = chi2array.sum(axis=0).argmin() + 1

    for i, file in enumerate(files):
        ra, ri, sigmax, sigmay = getRelInt(file)
        popt, perr, chi2 = multipoleFit(ra, ri, args.l, sigmay)
        #a = np.reshape(popt, (-1,2))
        a = np.reshape(popt[1:], (-1,2))
        amp[i], phase[i] = a[0]
        #e = np.reshape(perr, (-1,2))
        e = np.reshape(perr[1:], (-1,2))
        amp_err[i], phase_err[i] = e[0]

    # Deal with negative amplitudes/phases
    c0 = amp < 0
    amp[c0] *= -1.
    phase[c0] += 180
    phase[phase > 360] -= 360
    phase[phase < 0] += 360

    #print amp
    #print amp_err
    #print phase
    #print phase_err

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
        #ax.set_title('Amplitude', fontsize=16)
        ax.set_xlabel(r'$\mathrm{log}_{10}(E/\mathrm{GeV})$', fontsize=14)
        ax.set_ylabel(r'$\Delta N/\langle N \rangle$', fontsize=14)
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
        #ax.set_title('Phase', fontsize=16)
        ax.set_xlabel(r'$\mathrm{log}_{10}(E/\mathrm{GeV})$', fontsize=14)
        #ax.set_ylabel(r'Dipole phase $(\phi/\mathrm{deg})$', fontsize=14)
        ax.set_ylabel(r'Right Ascension $[^{\circ}]$', fontsize=14)
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





