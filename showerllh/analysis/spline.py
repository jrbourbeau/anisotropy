#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from icecube.photospline import spglam as glam

from llhtools import *
from prob import getProbs

##===========================================================================##
## Attempt to spline probability tables

## Look at our probability space ##
def pplot(s, name, spl=True, nk=2, npoly=3, clean=1):

    r = np.log10(s['ML_energy'])
    cut = s['cuts']['llh']
    Ebins, Emids = getEbins(), getEmids()

    f = plt.figure()
    ax1 = f.add_subplot(1,1,1)
    ax1.set_title(name)
    ax1.set_xlabel('true')
    ax1.set_ylabel('reco')

    p = getProbs(s, cut, clean=clean)

    st = len(Emids) - len(p[name][0])
    x, y = Emids, Emids[st:]
    X, Y = np.meshgrid(x, y)
    Z0 = np.log10(p[name])
    Z0[Z0==-np.inf] = -6
    if spl:
        Z1 = spline(s, p[name], nk=nk, npoly=npoly, plot=True)
    im = ax1.imshow(Z1-Z0, interpolation='bilinear',
            origin='lower', extent=[Ebins[st], 9.5, Ebins[0], 9.5])
    f.colorbar(im)
    plt.show()


## Return knots for a given set of points ##
def getKnots(pts, k, n):

    dist = float(pts[1] - pts[0])
    step = dist/k
    start = pts[0] - (n+1)*step
    end = pts[-1] + n*step
    num = np.around((end - start) / step) + 1
    knots = np.linspace(start, end, num=num)
    knots += step/10.
    return knots


## Attempt to spline probability tables ##
def spline(s, p0, nk=2, npoly=3, plot=False):

    Ebins, Emids = getEbins(), getEmids
    st = len(Emids) - len(p0[0])
    axes = [Emids, Emids[st:]]
    knots = [getKnots(axes[0], nk, npoly), getKnots(axes[1], nk, npoly)]

    # Increment all 0's by a fraction of the minimum value
    min_value = p0[p0!=0].min()
    min_value /= 100.
    min_mat = min_value * (p0==0)
    p0 += min_mat
    # Move to log space for smoother fit
    p0 = np.log10(p0)

    # Weight bins
    wts = np.ones(p0.shape)
    x = np.log10(s['MC_energy'][s['cuts']['llh']])
    y = np.log10(s['ML_energy'][s['cuts']['llh']])
    h, xedges, yedges = histogram2d(y, x, bins=(Ebins, Ebins[st:]))
    sigma = np.sqrt(np.var(h, axis=0))
    ntrue = np.sum(h, axis=0)

    # Deal with points with no hits
    weight = -1/p0
    minwt = weight.min()
    # Set the corners
    #weight[0][len(p0[0])-1] *= 10**(-4)
    #weight[len(p0)-1][0] *= 10**(-4)
    #weight[weight==minwt] = 0
    #weight[weight==minwt] *= 10**(-3)
    #weight /= sigma
    weight *= np.log10(ntrue)
    wts *= weight
    pen  = 10**(-4) * weight.min()

    spl = glam.fit(p0, wts, axes, knots, order=(npoly,npoly), periods=(0,0), 
            penalties={2:[pen,pen]})

    if plot:
        f = plt.figure()
        ax1 = f.add_subplot(1,1,1)
        ax1.set_title('')
        ax1.set_xlabel('true')
        ax1.set_ylabel('reco')

        x, y = [Emids, Emids[st:]]
        t0 = glam.grideval(spl, [x, y])
        X, Y = np.meshgrid(x, y)
        im = ax1.imshow(t0, interpolation='bilinear', vmin=-4.3, vmax=0,
                origin='lower', extent=[Ebins[st], 9.5, Ebins[0], 9.5])
        f.colorbar(im)
        plt.show()

    test = glam.grideval(spl, axes)

    return test


""" Check properties of splined matrices """
def spltest(s, bin=40, nk=2, npoly=3, clean=1):

    Emids = getEmids()
    cut = s['cuts']['llh']
    p = getProbs(s, cut, clean=clean)
    p0 = p['Rp|Tf']
    st = len(Emids) - len(p0[0])
    p1 = p['Rp|Tf'][:,bin]
    oldsum = p1.sum()

    axes = Emids
    knots = getKnots(axes, nk, npoly)

    # Increment all 0's by a fraction of the minimum value
    min_value = p1[p1!=0].min()
    min_value /= 100.
    min_mat = min_value * (p1==0)
    p1 += min_mat
    # Move to log space for smoother fit
    p1 = np.log10(p1)

    # Weight bins by probability
    wts = np.ones(p1.shape)
    weight = -1./p1
    minwt = weight.min()
    weight[weight==minwt] = 0       # bins with no data get 0 weight

    # Others weighted by variance
    #cut = s['cuts']['llh']
    #x = np.log10(s['MC_energy'][cut])
    #y = np.log10(s['ML_energy'][cut])
    #h, xedges, yedges = histogram2d(y, x, bins=(Ebins, Ebins[st:]))
    #sigma = np.sqrt(np.var(h, axis=0))[bin]
    #weight[weight!=0] /= sigma

    wts *= weight
    pen  = 10**(-8) * minwt

    spl = glam.fit(p1, wts, [axes], [knots], order=[npoly], 
            penalties={2:[pen]})

    test = glam.grideval(spl, [axes])
    newsum = (10**test).sum()

    print 'original', oldsum
    print 'splined', newsum

    f = plt.figure()
    ax1 = f.add_subplot(1,1,1)
    ax1.set_title('')
    ax1.set_xlabel('reco')
    ax1.set_ylabel('prob')

    plt.plot(Emids, test, '.', label='spline')
    plt.plot(Emids, p1, 'x', label='orig')

    #plt.legend(loc='upper left')
    plt.ylim(-7, 0)
    plt.show()

    #return test












