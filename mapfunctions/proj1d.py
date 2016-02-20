#!/usr/bin/env python

from mapFunctions import getMap

#from pylab import *
import numpy as np
import healpy as hp
import os, optparse

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.container import ErrorbarContainer
from scipy.optimize import curve_fit
from scipy import stats

#rc('xtick', labelsize=18)
#rc('ytick', labelsize=18)

#-------------------------------------------------------

class re_order_errorbarHandler(mpl.legend_handler.HandlerErrorbar):
    def create_artists(self, *args, **kwargs):
        a = mpl.legend_handler.HandlerErrorbar.create_artists(self, 
                *args, **kwargs)
        a = a[-1:] + a[:-1]
        return a


def returnRI(bgmap, datamap, **opts):

    # Setup right-ascension bins
    degree = np.pi / 180
    ramin = opts['ramin'] * degree
    ramax = opts['ramax'] * degree
    rabins = np.linspace(ramin, ramax, opts['nbins']+1)

    # Calculate phi for each pixel
    npix  = len(bgmap)
    nside = hp.npix2nside(npix)
    theta, phi = hp.pix2ang(nside, range(npix))
    # Bin in right ascension
    data = np.histogram(phi, bins=rabins, weights=datamap)[0]
    bg   = np.histogram(phi, bins=rabins, weights=bgmap)[0]

    with np.errstate(invalid='ignore'):
        ri = (data - bg) / bg
        sigmay = np.sqrt(data * (bg + 1./20 * data)/bg**3)
    dx = (ramax - ramin)/(2*opts['nbins'])
    ra = np.linspace(ramin+dx, ramax-dx, opts['nbins']) / degree
    sigmax = dx * np.ones(opts['nbins']) / degree

    return (ra, ri, sigmax, sigmay)


def getRelInt(file, **opts):

        # Get background and data maps
        bgmap = getMap(file, mapName='bg', **opts)
        datamap = getMap(file, mapName='data', **opts)
        # They'll be used as weights, so set unseen values to 0
        bgmap[bgmap==hp.UNSEEN] = 0
        datamap[datamap==hp.UNSEEN] = 0

        vals = returnRI(bgmap, datamap, **opts)
        return vals


def returnRI_new(relint_map, relerr_map, **opts):

    # Setup right-ascension bins
    deg2rad = np.pi / 180
    ramin = opts['ramin'] * deg2rad
    ramax = opts['ramax'] * deg2rad
    rabins = np.linspace(ramin, ramax, opts['nbins']+1)
    #print('rabins = {}'.format(rabins))

    # Calculate phi for each pixel
    npix  = len(relint_map)
    nside = hp.npix2nside(npix)
    theta, phi = hp.pix2ang(nside, range(npix))
    # Bin in right ascension
    phiBins = np.digitize(phi, rabins) - 1
    # UNSEEN cut
    cut = (relint_map != hp.UNSEEN)

    ri, sigmay = np.zeros((2,opts['nbins']))
    for i in range(opts['nbins']):
        phiCut = (phiBins == i)
        c0 = cut * phiCut
        ri[i] = np.mean(relint_map[c0])
        sigmay = np.sqrt(np.sum(relerr_map[c0]**2))/c0.sum()
        #sigmay = np.sqrt(data * (bg + 1./20 * data)/bg**3)
    dx = (ramax - ramin)/(2*opts['nbins'])
    ra = np.linspace(ramin+dx, ramax-dx, opts['nbins']) / deg2rad
    sigmax = dx * np.ones(opts['nbins']) / deg2rad

    return (ra, ri, sigmax, sigmay)


def getRelInt_new(file, **opts):
    # Get relative intensity map
    relint = getMap(file, mapName='relint', **opts)
    relerr = getMap(file, mapName='relerr', **opts)
    vals = returnRI_new(relint, relerr, **opts)
    return vals


def lineFit(x, y, sigmay):
    p = np.polyfit(x, y, 0, w=sigmay)
    fit = np.poly1d(p)
    fitVals = fit(np.asarray(y))
    chi2 = (1. / (len(y)-1)) * sum((y - fitVals)**2 / sigmay**2)
    return fitVals, p, chi2

def multipoleFit(x, y, l, sigmay):
    # Define sine function  (fix frequency)
    f0 = 2*np.pi / 360
    func = lambda x, *p: sum([p[3*i] * np.cos(f0/(i+1)*x + p[3*i+1])
            + p[3*i+2] for i in range(l)])
    # Guess at best fit parameters
    amplitude = (3./np.sqrt(2)) * np.std(y)
    phase     = 0
    offset    = np.mean(y)
    p0 = [amplitude, phase, offset]*l
    # Do best fit
    popt, pcov = curve_fit(func, x, y, p0, sigma=sigmay)
    print popt, np.sqrt(np.diagonal(pcov))
    fitVals = func(x, *popt)
    ndof  = len(popt)
    chi2 = sum((y - fitVals)**2 / sigmay**2)
    return fitVals, popt, chi2


def calcAvgBkg(ax, datamap, basename, **opts):

    nside = hp.get_nside(datamap)
    npix = hp.nside2npix(nside)

    counts = zeros(4*nside, dtype=float)
    norm = zeros(4*nside, dtype=float)

    ringId = hp.pix2ring(nside, np.array(range(npix)))
    for i, ring in enumerate(ringId):
        counts[ring] += datamap[i]
        norm[ring] += 1.

    bgmap = zeros(npix, dtype=float)
    mapmax = max(datamap)
    for i, ring in enumerate(ringId):
        bgmap[i] = counts[ring] / norm[ring]

    ra, ri, sigmax, sigmay = returnRI(bgmap, datamap, **opts)
    ax.errorbar(ra, ri, xerr=0*sigmax, yerr=sigmay, marker='o', fmt='.',\
            capsize=4, label=basename+" (avg)", linewidth=2, markersize=8, mew=0)


#-------------------------------------------------------
#-------------------------------------------------------

if __name__ == "__main__":

    # Set up command line options
    usage = "usage: %prog [options] INPUT.fits"
    parser = optparse.OptionParser(usage)
    parser.add_option("-r", "--ramin", dest="ramin", type=float,
            default=0, help="minimum RA")
    parser.add_option("-R", "--ramax", dest="ramax", type=float,
            default=360, help="maximum RA")
    parser.add_option("-D", "--decmax", dest="decmax", type=float,
            help="maximum Dec")
    parser.add_option("-d", "--decmin", dest="decmin", type=float,
            help="minimum Dec")
    parser.add_option("-n", "--nbins", dest="nbins", type=int,
            default=24, help="number of bins")
    parser.add_option("-z","--zero", action="store_true", dest="zeroline",
            default=False, help="Draw zero line")
    parser.add_option("-f","--flipra", action="store_true", dest="flipra",
            default=False, help="Flips RA in x axis")
    parser.add_option("--avg", action="store_true", dest="avgbg",
            default=False, help="Bkg from avg data")
    parser.add_option("-o", "--output", dest="output", default=None,
            help="Output image file")
    parser.add_option("-b","--batchmode", action="store_true", dest="batchMode",
            default=False, help="Execute without interaction")
    parser.add_option("-s","--syserr", action="store_true", dest="syserr",
            default=False, help="Draw systematic error bars")
    parser.add_option("--multi", dest='multi', type=int,
            default=None, help='Use multipole subtraction')
    parser.add_option('--fit', dest='fit',
            default=False, action='store_true',
            help='Show best-fit dipole on plot')
    parser.add_option("--labels", dest='labels',
            help='Custom label options built-in [configs, method]')
    parser.add_option("-v","--verbose", action="store_true", dest='verbose',
            default=False, help='Optional additional output')

    options, args = parser.parse_args()
    #print('args = {}'.format(args))
    #print('options = {}'.format(options))
    opts = vars(options).copy()
    #print('opts = {}'.format(opts))

    # Default masking behavior
    if not options.decmax:
        opts['mask'] = True

    if options.verbose:
        for key in sorted(opts.keys()):
            print ' --%s %s' % (key, opts[key])

    # Setup dictionaries for default values
    labels = {'sid':'sidereal', 'solar':'solar', 'anti':'anti-sidereal'}
    labels['ext'] = 'extended-sidereal'
    #print('labels = {}'.format(labels))
    errDict = {'sid':'anti', 'solar':'ext'}
    my_handler_map = {ErrorbarContainer: re_order_errorbarHandler(numpoints=1)}

    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111)
    p = {}

    for i, file in enumerate(args):

        # Labeling and possible error file
        basename = os.path.basename(file)[:-5]
        #print('basename = {}'.format(basename))
        label = basename
        if options.labels == 'configs':
            label = basename.split('_')[0]
        if options.labels == 'method':
            label = [labels[k] for k in labels if k in basename.split('_')][0]
        errfile = ''
        for key in errDict.keys():
            if key in file:
                errfile = file.replace(key, errDict[key])

        # Get RA, relint, and errors
        ra, ri, sigmax, sigmay = getRelInt_new(file, **opts)

        # Plot 1-D relative intensity projection
        l = ax.errorbar(ra, ri, xerr=0*sigmax, yerr=sigmay, marker='.', fmt='.',
                    capsize=4, label=label, linewidth=2, markersize=8, mew=0)
        if options.verbose:
            print "sigma: ", basename, np.std(ri)
            print "max: ", basename, np.max(np.abs(ri))

        # Optional best fit(s)
        #fitVals, popt, chi2 = lineFit(ra, ri, sigmay)
        fitVals, popt, chi2 = multipoleFit(ra, ri, 1, sigmay)
        if options.fit:
            ax.plot(ra, fitVals)

        # Print chi2 and p-value
        ndof = len(ri) - 1
        pvalue = 1 - stats.chi2.cdf(chi2, ndof)
        print '%.5f : %.3e' % (chi2, 1-pvalue)

        # Additional plotting options
        if options.avgbg:
            calcAvgBkg(ax, datamap, basename, **opts)
        if options.syserr and os.path.isfile(errfile):
            ra_err, ri_err, sigx, sigy = getRelInt_new(errfile, **opts)
            # Old method - fit sine function and take amplitude
            #fitVals, popt, chi2  = multipoleFit(ra_err, ri_err, 1, sigy)
            #amp = popt[0]
            # New method - just take max value
            syserr = ri_err.max()
            patches = [mpl.patches.Rectangle([ra[j]-5, ri[j]-syserr], 10, 
                    2*syserr, ec="none") for j in range(len(ra))]
            cln = PatchCollection(patches, cmap=mpl.cm.jet, 
                    alpha=0.5, color=l[0].get_color())
            ax.add_collection(cln)

    ax.set_xlabel(r"Right Ascension $[^{\circ}]$",fontsize=14)
    ax.set_ylabel(r"$\Delta N/\langle N \rangle$",fontsize=14)
    #ax.grid()

    if options.zeroline:
        xzero = np.arange(0, 360, 1)
        yzero = 0 * xzero
        ax.plot(xzero,yzero,linewidth=1.5,linestyle='--',color='black')

    ax.set_xlim(options.ramax, options.ramin)
    if options.flipra:
        ax.set_xlim(options.ramin, options.ramax)

    leg = ax.legend(loc='lower right', handler_map=my_handler_map)
    plt.draw()

    if options.output:
        plt.savefig(options.output, dpi=300, bbox_inches='tight')
    if not options.batchMode:
        plt.show()


