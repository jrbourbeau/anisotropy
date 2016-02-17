#!/usr/bin/env python

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import glob, argparse
from scipy import stats
from scipy.special import gamma, gammaln

import myGlobals as my
from useful import getMids

from anisotropy.mapFunctions import mapFunctions as mf
#from anisotropy.mapFunctions.plotFITS import SetupAbsThresholdColormap
from anisotropy.mapFunctions.plotFITS import SetupColorBar


def getSmallRange(x):
    maxbin = x.argmax()
    badbins = np.array([i for i, j in enumerate(x) if j < 5])
    try: redge = badbins[badbins>maxbin].min()
    except ValueError:
        redge = len(x)
    try: ledge = badbins[badbins<maxbin].max()
    except ValueError:
        ledge = 0
    return ledge, redge


# Compare shape of flux for each bin to ensemble average
def bayesFactor(x, xerr, args):

    output = np.zeros(len(x), dtype=np.double)
    # Calculate totals
    xtot = np.sum(x, axis=0)
    ntot = np.sum(xtot)
    for i, x_i in enumerate(x):
        n_i = np.sum(x_i)
        #if n_i == 0:
        #    continue
        output[i] = (sum(gammaln(x_i+1) + gammaln(xtot+1) - gammaln(x_i+xtot+2))
                + gammaln(n_i+ntot+2) - gammaln(n_i+1) - gammaln(ntot+1))

    return output

# Compare shape of flux for each bin to ensemble average
def bayes2(x1, x2):

    # Calculate totals
    n1 = np.sum(x1)
    n2 = np.sum(x2)
    output = (sum(gammaln(x1+1) + gammaln(x2+1) - gammaln(x1+x2+2)) 
            + gammaln(n1+n2+2) - gammaln(n1+1) - gammaln(n2+1))

    return output


def bayes3(x, xerr, args):

    output = np.zeros(len(x), dtype=np.double)
    xtot = np.sum(x, axis=0)
    ntot = np.sum(xtot)
    n_i = np.sum(x, axis=1)
    c0 = (n_i == 0)
    n_i = n_i[c0]
    output[c0] = np.sum(gammaln(x[c0]+1) + gammaln(xtot+1)
            - gammaln(x[c0]+xtot+2) + gammaln(n_i+ntot+2)
            - gammaln(n_i+1) - gammaln(ntot+1), axis=1)

    return output


def getTest(x, xerr, args):

    if args.name in ['llh','llhcut']:
        bins = np.linspace(-20, 20, 151)
    if args.name == 'energy':
        bins = np.arange(6.25, 9.501, 0.05)
    mids = getMids(bins)
    with np.errstate(invalid='ignore'):
        x = np.sum((mids * x) / np.sum(x, axis=1)[:,np.newaxis], axis=1)
    x[x!=x] = 0.
    ave_x = x.mean()
    return x - ave_x


def getChi2(x, err, args):

    # Get average distribution
    ave_x = np.sum(x, axis=0)
    ave_err = np.sqrt(np.sum(err, axis=0))
    # Option to limit range
    if args.small:
        ledge, redge = getSmallRange(ave_x)
        ave_x = ave_x[ledge:redge]
        ave_err = ave_err[ledge:redge]

    # Normalize
    ave_xtot = float(ave_x.sum())
    xtot = np.sum(x, axis=1).astype(float)
    if args.small:
        xtot = np.sum(x[:,ledge:redge], axis=1).astype(float)

    with np.errstate(invalid='ignore'):
        ave_x /= ave_xtot
        ave_err /= ave_xtot
        temp_x = x/xtot[:,np.newaxis]
        temp_err = np.sqrt(err)/xtot[:,np.newaxis]

    if args.small:
        temp_x = temp_x[:,ledge:redge]
        temp_err = temp_err[:,ledge:redge]

    # Bins w/ 0 events -> inf err
    nancut = (ave_err == 0)
    ave_err[nancut] = np.inf
    temp_err[:,nancut] = np.inf

    with np.errstate(invalid='ignore'):
        chi2 = np.sum((temp_x - ave_x)**2 / (temp_err**2 + ave_err**2), axis=1)
    ndof = temp_x.shape[1] - 1

    return chi2, ndof


def smoothMap(map, smooth=5.):

    if smooth == 0:
        return map

    npix = len(map)
    nside = hp.npix2nside(npix)
    smooth_rad = smooth * np.pi/180.
    smooth_map = np.zeros(map.shape)

    vec = np.transpose(hp.pix2vec(nside, np.arange(npix)))
    for i in range(npix):
        neighbors = hp.query_disc(nside, vec[i], smooth_rad)
        smooth_map[i] += np.sum(map[neighbors], axis=0)

    return smooth_map


def loadHists(config=None):

    my.setupShowerLLH(verbose=False)
    histFiles = glob.glob('%s/*_data/IT*sky.npy' % my.llh_data)
    if config != None:
        histFiles = glob.glob('%s/%s_data/IT*sky.npy' % (my.llh_data, config))

    h = {}
    for hfile in histFiles:
        htemp = np.load(hfile)
        htemp = htemp.item()
        for key in htemp.keys():
            h[key] = htemp[key]

    return h


def histReader(h, name):

    keyList = [k for k in h.keys() if name in k and 'err' not in k]
    errList = [k for k in h.keys() if name in k and 'err' in k]
    if '_w' not in name:
        keyList = [k for k in keyList if '_w' not in k]
        errList = [k for k in errList if '_w' not in k]
    x = np.sum([h[k] for k in keyList], axis=0)
    xerr = np.sum([h[k] for k in errList], axis=0)

    stbin = 25
    if name == 'energy':
        x = x[:,stbin:]
        xerr = xerr[:,stbin:]

    return x, xerr


def singlePlot(x, xerr, args):

    fig, ax = plt.subplots()
    if args.name in ['llh','llhcut']:
        bins = np.linspace(-20, 20, 151)
        xyDict = {'Lighter':(33,97), 'Heavier':(33, 269)}
    if args.name == 'energy':
        bins = np.arange(6.25, 9.501, 0.05)
        xyDict = {'Excess':(150, 185)}
    mids = getMids(bins)
    degree = np.pi / 180.
    npix = len(x)
    nside = hp.npix2nside(npix)

    for i, key in enumerate(xyDict.keys()):
        # Plot distribution for pixel
        thetaPix, phiPix = xyDict[key][0]*degree, xyDict[key][1]*degree
        pix = hp.ang2pix(nside, thetaPix, phiPix)
        normFactor = float(x[pix].sum())
        ax.errorbar(mids, x[pix]/normFactor, 
                yerr=np.sqrt(xerr[pix])/normFactor, label=key)
        theta, phi = hp.pix2ang(nside, range(npix))
        thetaRing, phiRing = hp.pix2ang(nside, pix)
        thetaCut = (theta == thetaRing)
        print key, ':', bayes2(x[pix], np.sum(x[thetaCut], axis=0))

    # Plot distribution for rest of declination band
    normFactor = float(x[thetaCut].sum())
    ax.errorbar(mids, np.sum(x[thetaCut], axis=0)/normFactor, 
            yerr=np.sqrt(np.sum(xerr[thetaCut], axis=0))/normFactor, 
            label='Average')
    #ledge, redge = getSmallRange(x[pix])
    #ax.errorbar(mids[ledge:redge], x[pix][ledge:redge], 
    #        yerr=np.sqrt(xerr[pix][ledge:redge]), label=key)

    ax.set_xlabel(r'Energy ($\log_{10}(E/\mathrm{GeV})$)')
    ax.set_ylabel(r'Normalized Counts')
    ax.set_yscale('log')

    plt.legend()
    if args.out:
        plt.savefig(args.out, dpi=300, bbox_inches='tight')
    if not args.batch:
        plt.show()


def makeTitle(args):

    title = '%s %s %ddeg' % (args.name, args.maptype, args.smooth)
    if args.config != None:
        title = '%s %s' % (args.config, title)
    if args.decbands:
        title += ' decbands'

    return title


def checker(nguess):

    from scipy.optimize import fsolve
    x = 10**np.arange(6.2, 9.501, 0.05)

    def myspec(ntot, dg, g0=2.7):
        g1 = g0 + dg
        if ntot < 0:
            return -10**5
        y0 = np.random.poisson(ntot * (x**(-g0) / (x**-g0).sum()))
        y1 = np.random.poisson(ntot * (x**(-g1) / (x**-g1).sum()))
        bf = bayes2(y0, y1)
        return bf

    fig, ax = plt.subplots()
    dgs = np.arange(0.2, 1.501, 0.05)
    for bf in [100, 1000, 10000]:
        n = []
        for dg in dgs:
            myfunc = lambda ntot: np.log(bf) - myspec(ntot, dg)
            #nguess = 10**4
            n += [fsolve(myfunc, nguess)]
            #ave_tot = bayesFactor(np.array([y1,y2]), [], [])
        ax.plot(dgs, n, label=bf)

    ax.set_yscale('log')
    ax.legend()
    plt.show()


def myPlot(x, xerr, args):

    # General setup
    npix = len(x)
    nside = hp.npix2nside(npix)
    theta, phi = hp.pix2ang(nside, range(npix))
    degree = np.pi / 180.
    fDict = {'chi2':getChi2, 'pvals':getChi2, 'bayes':bayesFactor, 'test':getTest}
    f = fDict[args.maptype]

    # Compare each pixel to ensemble distribution for given declination
    map = np.zeros(npix, dtype=np.double)
    if args.decbands:
        thetaRings = np.array(sorted(list(set(theta))))
        thetaRings = thetaRings[thetaRings > (180-args.thetaMax)*degree]
        for i, th in enumerate(thetaRings):
            c0 = (theta == th)
            map[c0] += f(x[c0], xerr[c0], args)
    else:
        map = f(x, xerr, args)

    # Option to calculate p-values
    #if args.maptype == 'pvals':
    #    pvals = np.zeros(len(chi2), dtype=np.double)
    #    for i, chi in enumerate(chi2):
    #        if chi == 0:
    #            continue
    #        pvals[i] = stats.chisqprob(chi2[i], ndof[i])
    #    map = np.log10(pvals)

    # Mask
    map = mf.maskMap(map, args.thetaMin-90, args.thetaMax-90)
    if args.scale:
        map[map!=hp.UNSEEN] *= (10**args.scale)

    temp = {'Excess':(150,185)}

    title = makeTitle(args)

    ### SETUP FROM PLOTFITS ###
    unmasked = np.array([i for i in map if (i!=hp.UNSEEN and i!=np.inf)])
    min = float(args.min) if args.min else unmasked.min()
    max = float(args.max) if args.max else unmasked.max()
    if not args.min:
        args.min = '%.2f' % min
    if not args.max:
        args.max = '%.2f' % max

    # Setup colormap with option for threshold
    colormap = plt.get_cmap('jet')
    colormap.set_under('white')
    colormap.set_bad('gray')

    pltParams = {'fig':1, 'rot':[0, -90, 180], 'title':'', 'min':min, \
            'max':max, 'cbar':False, 'notext':True, 'coord':'C', \
            'cmap':colormap}
    ###

    if args.polar:
        #hp.orthview(map, rot=(0,90,180), title=title, half_sky=True)
        hp.orthview(map, half_sky=True, **pltParams)
        for key in temp.keys():
            hp.projtext(temp[key][0]*degree, temp[key][1]*degree, key)
    else:
        #hp.mollview(map, rot=(0,0,180), title=title)
        #hp.mollview(map, rot=(0,0,180), title='')
        hp.mollview(map, rot=180, title='')
    #for key in temp.keys():
    #    hp.projtext(temp[key][0]*degree, temp[key][1]*degree, key)
    hp.graticule(verbose=False)

    ### MORE SETUP FROM PLOTFITS ###
    labelDict = {'bayes':r'Bayes Factor [$\ln(B_{21})$]'}
    labelDict['test'] = r'Difference from Mean [$\overline{D_i} - \overline{D_{tot}}$]'
    label = labelDict[args.maptype]
    if args.scale:
        label += ' [x 10$^{-%d}$]' % args.scale
    cbarParams = {'label':label, 'min':args.min, 'max':args.max, \
            'coord':'C', 'fontsize':'small', \
            'projaxis':hp.projaxes.HpxOrthographicAxes}
    SetupColorBar(**cbarParams)


    if args.out != None:
        plt.savefig(args.out, dpi=300, bbox_inches='tight')

    if not args.batch:
        plt.show()


def newPlot(x, xerr, args):

    # General setup
    npix = len(x)
    nside = hp.npix2nside(npix)
    theta, phi = hp.pix2ang(nside, range(npix))
    degree = np.pi / 180.

    cuts = [phi < np.pi, phi >= np.pi]
    thcut = (theta > (180-args.thetaMax)*degree) * (theta >= 20*degree)
    xtemp = []
    for cut in cuts:
        c0 = cut * thcut
        xtemp += [np.sum(x[c0], axis=0)]
    print bayes2(xtemp[0], xtemp[1])


if __name__ == "__main__":

    p = argparse.ArgumentParser()
    p.add_argument('-c', '--config', dest='config',
            help='Detector configuration')
    p.add_argument('-n', '--name', dest='name',
            default='llhcut', choices=['energy','dist','llh','llhcut'],
            help='Desired quantity to plot')
    p.add_argument('--maptype', dest='maptype',
            choices=['bayes','chi2','single','pvals','test'],
            help='Desired plot type')
    p.add_argument('-s', '--smooth', dest='smooth', type=float,
            default=5.,
            help='Top-hat smoothing angle')
    p.add_argument('-m', '--min', dest='min')
    p.add_argument('-M', '--max', dest='max')
    p.add_argument('--scale', dest='scale', type=float)
    p.add_argument('--decbands', dest='decbands',
            default=False, action='store_true',
            help='Option to calculate map by declination bands')
    p.add_argument('--small', dest='small',
            default=False, action='store_true',
            help='Limit range of distributions')
    p.add_argument('--polar', dest='polar',
            default=False, action='store_true',
            help='Plot in polar coords')
    p.add_argument('--thetaMax', dest='thetaMax', type=float,
            default=37.9,
            help='Maximum angle (degrees)')
    p.add_argument('--thetaMin', dest='thetaMin', type=float,
            default=0,
            help='Minimum angle (degrees)')
    p.add_argument('--checker',
            default=False, action='store_true',
            help='')
    p.add_argument('--new',
            default=False, action='store_true',
            help='')
    p.add_argument('-b', '--batch', dest='batch',
            default=False, action='store_true',
            help='')
    p.add_argument('-o', '--out', dest='out',
            help='')
    args = p.parse_args()

    #histSpectrum('IT81')

    if args.checker:
        checker()

    # Load data
    h = loadHists(args.config)
    x, xerr = histReader(h, args.name)

    mask = (x != 0)

    # Smooth
    x = smoothMap(x, args.smooth)
    xerr = smoothMap(xerr, args.smooth)

    print x[mask].min(), x[mask].max(), np.percentile(x[mask], 50)
    raise

    # Option to plot a single bin
    if args.maptype == 'single':
        singlePlot(x, xerr, args)

    # Option to plot chi-square
    if args.maptype in ['chi2', 'pvals', 'bayes', 'test']:
        myPlot(x, xerr, args)

    if args.new:
        newPlot(x, xerr, args)

