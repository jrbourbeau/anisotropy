#!/usr/bin/env python

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import matplotlib as mpl

import myGlobals as my
import simFunctions_IC as simFunctions
from plotFunctions import cmap_discretize
from useful import getMedian, getMids
degree = np.pi / 180.

import os, time, glob

##============================================================================
## Basic functions

""" Get the simulation basename for anisotropy studies """
def getSimBase(config):

    # Setup global path names
    my.setupAnisotropy(verbose=False)
    # Get base for simulation path name
    sim = simFunctions.cfg2sim(config)
    simBase = '%s/%s_%s' % (my.ani_sim, config, sim)
    return simBase


def getSimFiles(config):
    simBase = getSimBase(config)
    simFiles = ['%s.npy' % simBase]
    if not os.path.isfile(simFiles[0]):
        simFiles = sorted(glob.glob('%s_Part*.npy' % simBase))
    return simFiles


""" Standard energy bins """
def getEbins():
    ebins = [4, 4.25, 4.5, 4.75, 5, 5.25, 5.5, 6, 6.5]
    return ebins

def test(s):
    ebins = getEbins() + [100]
    typeList = sorted(list(set(s['type'])))
    typeBreaks = [1000020040, 1000060120, 1000130270, 1000220480]
    types = ['p','he','cno','al','fe']
    typeBins = np.digitize(s['type'], typeBreaks)
    t = {}
    for i, type in enumerate(types):
        t[type] = typeBins==i

    bins = np.linspace(3, 9.5, 1000)
    for i in range(len(ebins)-1):
        fig, ax = plt.subplots()
        emin = ebins[i]
        emax = ebins[i+1]
        ax.set_title('%s to %s' % (emin, emax))
        c0 = (s['energy_sp'] >= emin) * (s['energy_sp'] < emax)
        for key in t:
            c1 = c0 * t[key]
            ax.hist(np.log10(s['energy'][c1]), bins=bins, 
                    histtype='step', log=True, label=key, normed=True)
        ax.legend()

    plt.show()


""" Load simulation fro a given configuration """
def load_sim(config, spline=True, part=None):

    # Get simulation file names
    simFiles = getSimFiles(config)

    # Option to limit to single part file
    if part != None:
        simFiles = [simFiles[part]]

    # Load simulation
    for i, simFile in enumerate(simFiles):
        print 'Working on %s' % simFile
        if i == 0:
            s = np.load(simFile)
            s = s.item()
        else:
            temp = np.load(simFile)
            temp = temp.item()
            for key in temp.keys():
                old = len(s[key])
                new = len(temp[key])
                s[key].resize(old+new)
                s[key][old:] = temp[key]

    # Calculate reconstructed energies from spline table
    if 'reco2' in s.keys():
        theta, phi = hp.pix2ang(1024, s['reco2'].astype('int'))
        s['dstZenith']  = theta
        s['dstAzimuth'] = phi
    x = np.cos(s['dstZenith'])
    with np.errstate(divide='ignore'):
        y = np.log10(s['nchannel'])
    if spline:
        s['energy_sp'] = getSplineEnergies(config, x, y)

    return s


""" Calculate the spline energies for a given detector configuration. """
def getSplineEnergies(config, x, y):

    # Load spline table
    from icecube import photospline
    simBase = getSimBase(config)
    splineFile = '%s_spline.fits' % simBase
    s = photospline.I3SplineTable(splineFile)

    t0 = time.time()
    # Calculate median energies according to spline table
    def spline_energy(x, y):
        if (x < 0.3):
            return 0
        return s.eval([x, y])
    esplines = np.array(map(spline_energy, x, y))
    print 'Time to get splines:', time.time()-t0

    return esplines


""" Read in energy distribution information from file """
def readDist(config, emin, emax, old=False, t=False):

    # Read in distribution information
    infoFile = '/data/user/fmcnally/anisotropy/sim/eDist_IC.txt'
    if old:
        infoFile = '/data/user/fmcnally/anisotropy/sim/eDist_IC_orig.txt'
    with open(infoFile, 'r') as f:
        lines = f.readlines()
    table = [line.strip().split(' ') for line in lines]

    # Limit to desired data point
    data = [line[1:] for line in table if line[0]==config]
    if config == 'IC':
        data = [line[1:] for line in table if line[0][:2]==config]
    emin, emax = str(float(emin)), str(float(emax))
    data = np.array([l for l in data if l[:2]==[emin,emax]], dtype='float')

    median = data[:,2]
    sigL = data[:,2] - data[:,3]
    sigR = data[:,4] - data[:,2]

    # Create combined distribution (start unweighted)
    ave_median = np.average(median)
    ave_sigL = np.average(sigL)
    ave_sigR = np.average(sigR)

    if t:
        weights = 1./(sigL + sigR)
        ave_median = np.average(median, weights=weights)

    return ave_median, ave_sigL, ave_sigR


##============================================================================
## Plotting functions

""" Creates 2D hist binned in zenith and nchannel for energy reconstruction """
def recoHist(s, nx=100, ny=50, minCount=1, ztype='reco'):

    xbins = np.linspace(0.3, 1, nx+1)
    ybins = np.linspace(1, 3.2, ny+1)
    h = [[[] for j in range(ny+1)] for i in range(nx+1)]

    x = np.cos(s['zen']) if ztype=='true' else np.cos(s['dstZenith'])
    with np.errstate(divide='ignore'):
        y = np.log10(s['nchannel'])
        e = np.log10(s['energy'])

    # Fill histogram
    xbin = np.digitize(x, xbins) - 1
    ybin = np.digitize(y, ybins) - 1
    for i, energy in enumerate(e):
        h[xbin[i]][ybin[i]] += [energy]

    # Eliminate under/overflow bins
    h = h[:-1]
    h = [i[:-1] for i in h]

    for i in range(len(h)):
        for j in range(len(h[i])):
            if len(h[i][j]) < minCount:
                h[i][j] = [0]

    return h


""" Returns table of spline fits to histogram """
def getSplineTable(h, pens=1e-3):

    from icecube.photospline import spglam as glam
    # Bin setup
    nx, ny = len(h), len(h[0])
    xbins = np.linspace(0.3, 1, nx+1)
    ybins = np.linspace(1, 3.2, ny+1)

    # Extract desired information from histogram
    meds = [[np.percentile(h[i][j], 50) for j in range(ny)] for i in range(nx)]
    vars = [[np.var(h[i][j]) for j in range(ny)] for i in range(nx)]
    meds, vars = np.array(meds), np.array(vars)
    vars[vars==0] = np.inf
    w = 1/vars

    # Spline parameter setup
    axes, knots = [],[]
    binList = [xbins, ybins]
    nKnots = 30
    step_scale = 2/5.
    for bins in binList:
        mids = (bins[1:] + bins[:-1]) / 2.
        axes += [mids]
        step = (bins.max() - bins.min()) * step_scale
        knots += [np.linspace(bins.min()-step, bins.max()+step, nKnots)]

    # Calculate actual spline table
    tab = glam.fit(meds, w, axes, knots, order=(4), penalties={2:pens})

    return tab

""" Calculate energy resolution from full simulation info.
    NOTE: slow - recalculates histogram and spline tables each time
"""
def eres(s, nx=60, ny=40, minCount=1, ztype='reco', pens=1e-3):

    from icecube import photospline
    from icecube.photospline import splinefitstable

    fig, ax = plt.subplots()

    # Get histogram and spline information
    h = recoHist(s, nx, ny, minCount, ztype)
    tab = getSplineTable(h, pens=pens)

    # Binning setup
    xbins = np.linspace(0.3, 1, nx+1)
    ybins = np.linspace(1, 3.2, ny+1)
    x = np.cos(s['dstZenith'])
    with np.errstate(divide='ignore'):
        y = np.log10(s['nchannel'])
        eTrue = np.log10(s['energy'])

    # Calculate spline values
    splinefitstable.write(tab, 'testSpline.fits')
    tab2 = photospline.I3SplineTable('testSpline.fits')
    os.remove('testSpline.fits')
    c0 = (x >= 0.3)
    x1, y1, eTrueSpl = x[c0], y[c0], eTrue[c0]
    def spline_energy(x, y):
        return tab2.eval([x, y])
    print 'Calculating spline energies...'
    eSpline = np.array(map(spline_energy, x1, y1))

    emin, emax = eTrueSpl.min(), eTrueSpl.max()
    ebins = np.array([emin] + getEbins() + [emax])
    emids = (ebins[1:] + ebins[:-1]) / 2.

    # Bin and extract info
    x = eSpline
    values = eSpline - eTrueSpl
    medians, sigL, sigR, vars = getMedian(x, values, ebins)
    ax.errorbar(emids, medians, yerr=(sigL, sigR), fmt='.')

    ax.set_title('Energy Resolution vs True Energy')
    ax.set_xlabel('True Energy (log10(E/GeV))')
    ax.set_ylabel('Median Reconstructed Energy (log10(E/GeV))')

    plt.legend()
    plt.show()


""" Plot binned histogram and corresponding spline fit """
def splinePlot(s, nx=100, ny=50, minCount=1, ztype='reco', pens=1e-3):

    from icecube.photospline import spglam as glam
    # Get spline information
    h = recoHist(s, nx, ny, minCount, ztype)
    tab = getSplineTable(h, pens=pens)
    # Extract desired information from histogram
    meds = [[np.percentile(h[i][j], 50) for j in range(ny)] for i in range(nx)]
    meds = np.array(meds)

    # Energy splits for colorbar
    emin, emax = meds[meds!=0].min(), meds.max()
    ebins = [emin] + getEbins + [emax]
    xbins = np.linspace(0.3, 1, nx+1)
    ybins = np.linspace(1, 3.2, ny+1)
    xmids = (xbins[1:] + xbins[:-1]) / 2.
    ymids = (ybins[1:] + ybins[:-1]) / 2.

    # Plot with custom colormap
    fig = plt.figure(figsize=(17,6))
    mpl.rc('font', family='serif')
    cmap = plt.cm.jet
    cmap = cmap_discretize(cmap, ebins)
    cmap.set_under('white')

    # Subplot with binned medians
    ax = fig.add_subplot(121)
    X, Y = np.meshgrid(xmids, ymids)
    p = ax.pcolor(X, Y, meds.T, cmap=cmap, vmin=emin, vmax=emax)
    ax.set_title('Median Energy')
    ax.set_xlabel('cos(zenith)')
    ax.set_ylabel('log10(Nchannel)')
    ax.set_xlim(xbins.min(), xbins.max())
    ax.set_ylim(ybins.min(), ybins.max())

    fitaxes = [[],[]]
    fitaxes[0] = np.linspace(xbins.min(), xbins.max(), 3*nx+1)
    fitaxes[1] = np.linspace(ybins.min(), ybins.max(), 3*ny+1)
    fitX, fitY = np.meshgrid(fitaxes[0], fitaxes[1])
    fit = glam.grideval(tab, fitaxes)

    # Subplot with spline fit
    ax = fig.add_subplot(122)
    p = ax.pcolor(fitX, fitY, fit.T, cmap=cmap, vmin=emin, vmax=emax)
    cb = fig.colorbar(p, ax=ax)
    ax.set_title('Median Energy (Spline Fit)')
    ax.set_xlabel('cos(zenith)')
    ax.set_ylabel('log10(Nchannel)')
    ax.set_xlim(xbins.min(), xbins.max())
    ax.set_ylim(ybins.min(), ybins.max())

    plt.show()


def ang_hist(s, nbins=40):

    # Get angular information
    r_theta = s['dstZenith']
    r_phi   = s['dstAzimuth']
    t_theta, t_phi = s['zen'], s['azi']

    # Conversion to cartesian
    r_x = np.sin(r_theta) * np.cos(r_phi)
    r_y = np.sin(r_theta) * np.sin(r_phi)
    r_z = np.cos(r_theta)
    t_x = np.sin(t_theta) * np.cos(t_phi)
    t_y = np.sin(t_theta) * np.sin(t_phi)
    t_z = np.cos(t_theta)
    dphi = np.arccos(r_x*t_x + r_y*t_y + r_z*t_z)
    dphi /= degree

    qcut = (r_theta / degree <= 65)

    fig, ax = plt.subplots()
    ax.hist(dphi[qcut], bins=nbins,
            log=True, histtype='step', normed=True)
    plt.show()


def ang_hist2(s):

    # Get angular information
    r_theta = s['dstZenith']
    r_phi   = s['dstAzimuth']
    t_theta, t_phi = s['zen'], s['azi']

    dtheta = (r_theta - t_theta) / degree
    dphi = (r_phi - t_phi) / degree
    thetabins = np.arange(-50,50,1)
    phibins = np.arange(-50,50,1)

    qcut = (r_theta / degree <= 65)
    h, xedges, yedges = np.histogram2d(dtheta[qcut], dphi[qcut], 
            bins=[thetabins, phibins])
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

    fig, ax = plt.subplots()
    im = ax.imshow(np.log10(h.T), extent=extent, 
            origin='lower', interpolation='none')
    fig.colorbar(im)
    #ax.hist2d(dtheta[qcut], dphi[qcut], bins=[thetabins,phibins], log=True)
    plt.show()


""" Plot angular resolution """
def ang_res(s, xaxis='energy', nbins=10, rt='t', out=False, logPlot=False):

    # Extract info from simulation file
    if rt == 't':
        energy = np.log10(s['energy'])
        #energy = np.log10(s['nchannel'])
        zenith = np.cos(s['zen'])
    if rt == 'r':
        energy = s['energy_sp']
        zenith = np.cos(s['dstZenith'])

    # Get angular information
    r_theta = s['dstZenith']
    r_phi   = s['dstAzimuth']
    t_theta, t_phi = s['zen'], s['azi']

    # Conversion to cartesian
    r_x = np.sin(r_theta) * np.cos(r_phi)
    r_y = np.sin(r_theta) * np.sin(r_phi)
    r_z = np.cos(r_theta)
    t_x = np.sin(t_theta) * np.cos(t_phi)
    t_y = np.sin(t_theta) * np.sin(t_phi)
    t_z = np.cos(t_theta)

    # Quality cut
    qcut = (r_theta <= 65*degree)

    # Setup bins
    if xaxis == 'energy':
        bins = np.linspace(energy[qcut].min(), energy[qcut].max(), nbins+1)
    if xaxis == 'zenith':
        thetamax = 65. * degree
        bins = np.linspace(1, np.cos(thetamax), nbins+1)[::-1]

    # Store median and standard deviation
    x = getMids(bins)
    x1 = {'energy':energy, 'zenith':zenith}
    x1 = x1[xaxis]
    y = np.arccos(r_x*t_x + r_y*t_y + r_z*t_z)
    medians, sigL, sigR, vars = getMedian(x1[qcut], y[qcut], bins)

    # Setup plot
    titleDict = {'energy':'Energy', 'zenith':'Zenith Angle'}
    xDict = {'energy':'Log10(Energy/GeV)', 'zenith':'Zenith (degrees)'}
    rtDict = {'r':'Reconstructed', 't':'True'}

    fig, ax = plt.subplots()
    ax.set_title('Angular Resolution vs %s %s' % (rtDict[rt], titleDict[xaxis]),
            fontsize=18)
    ax.set_xlabel(xDict[xaxis], fontsize=16)
    ax.set_ylabel(r'$\Delta \psi$', fontsize=16)

    # Plot parameters
    lw = 2
    ms = 7*lw

    medians /= degree
    sigL /= degree
    sigR /= degree

    ax.errorbar(x, medians, yerr=[sigL,sigR], fmt='.', ms=ms)
    print medians, sigL, sigR
    if logPlot:
        ax.set_yscale('log')
    if out:
        plt.savefig(out)
    plt.show()


""" Plot energy reconstruction histogram """
def plotHist(s, nx=100, ny=50, minCount=1, ztype='reco'):

    xbins = np.linspace(0.3, 1, nx+1)
    ybins = np.linspace(1, 3.2, ny+1)
    xmids = (xbins[1:] + xbins[:-1]) / 2.
    ymids = (ybins[1:] + ybins[:-1]) / 2.

    h = recoHist(s, nx, ny, minCount, ztype)
    medians = [[np.percentile(h[i][j], 50) for j in range(ny)] for i in range(nx)]
    medians = np.array(medians)

    emin, emax = medians[medians!=0].min(), medians[medians!=0].max()
    ebins = [emin] + getEbins() + [emax]

    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111)
    X, Y = np.meshgrid(xbins, ybins)
    cmap = plt.cm.jet
    cmap = cmap_discretize(cmap, ebins)
    cmap.set_under('white')

    p = ax.pcolor(X, Y, medians.T, cmap=cmap, vmin=emin, vmax=emax)
    cb = fig.colorbar(p, ax=ax, ticks=ebins)
    cb.ax.set_yticklabels(['%.2f' % ebin for ebin in ebins])
    cb.set_label(r'$log_{10}(E/GeV)$', rotation=270, labelpad=20)
    ax.set_title('Energy vs Zenith and Nchannel')
    ax.set_xlabel(r'$cos(zenith)$')
    ax.set_ylabel(r'$log_{10}(Nchannel)$')
    ax.set_xlim(xbins.min(), xbins.max())
    ax.set_ylim(ybins.min(), ybins.max())

    plt.show()
    return h


def eres2(s, minCount=1, ztype='reco', zero=False, out=False):

    fig, ax = plt.subplots()
    mpl.rc('font', family='serif')
    ax.set_title('Energy Distributions of Energy Cuts')
    ax.set_xlabel(r'Reconstructed Energy Range ($log_{10}(E/GeV)$)')
    ax.set_ylabel(r'True Energy Distribution ($log_{10}(E/GeV)$)')

    eSpline = s['energy_sp']
    with np.errstate(divide='ignore'):
        eTrue = np.log10(s['energy'])
    c0 = (eSpline != 0)
    eSpline, eTrue = eSpline[c0], eTrue[c0]

    emin, emax = eSpline.min(), eSpline.max()
    ebins = np.array(getEbins() + [emax])
    emids = (ebins[1:] + ebins[:-1]) / 2.
    xL = emids - ebins[:-1]
    xR = ebins[1:] - emids
    # Setup for arrow
    xR[-1] = 0

    lw = 2
    ms = 7*lw

     # Bin and extract info
    x = eSpline
    values = eTrue
    if zero:
        values = eSpline - eTrue
    medians, sigL, sigR, vars = getMedian(x, values, ebins)
    ax.errorbar(emids, medians, xerr=(xL, xR), yerr=(sigL, sigR), \
            fmt='b.', lw=lw, ms=ms)
    ax.arrow(emids[-1], medians[-1], 0.5, 0.0, fc='b', ec='b',\
            head_width=.1, head_length=.2)

    ax.set_xlim(3.5, 8.0)

    if out:
        outPrefix = '/home/fmcnally/papers/anisotropy/figures'
        outFile = '%s/IC86_Energy_Distribution' % (outPrefix)
        #plt.savefig(outFile+'.eps', dpi=300, bbox_inches='tight')
        plt.savefig(outFile+'.png', dpi=300, bbox_inches='tight')

    #plt.legend(loc='lower right')
    plt.show()
