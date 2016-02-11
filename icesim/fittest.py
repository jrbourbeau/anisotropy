#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import argparse

from analysis import load_sim
from useful import getMids
degree = np.pi / 180.


def ang_hist(s, nbins=50000):

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
    #dphi = r_x*t_x + r_y*t_y + r_z*t_z
    dphi = np.arccos(r_x*t_x + r_y*t_y + r_z*t_z)
    dphi /= degree

    qcut = (r_theta / degree <= 65)

    bins = np.arccos(np.linspace(1, -1, nbins+1)) / degree
    h = np.histogram(dphi[qcut], bins=bins)

    #fig, ax = plt.subplots()
    #ax.hist(dphi[c0], bins=bins, histtype='step')
    #plt.show()

    return h


def gauss(x, p):
    return p[0] * np.exp(-x**2 / (2*p[1]**2))
def kent(x, p):
    sig_rad = p[1] * degree
    x_rad = x * degree
    k = 1. / sig_rad**2
    #c = k / (4*np.pi * np.sinh(k))
    #return c * np.exp(k * np.cos(x_rad))
    # Small angle approximation
    #return k/(2*np.pi) * np.exp(k * (np.cos(x_rad)-1))
    # Floating amplitude
    return p[0] * np.exp(k * (np.cos(x_rad)-1))


def myfunc(x, *p):

    if len(p) == 2:
        return kent(x, p)
    if len(p) == 4:
        g1 = gauss(x, p[:2])
        g2 = gauss(x, p[2:])
        return g1 + g2
        #k1 = kent(x, p[:2])
        #k2 = kent(x, p[2:])
        #return k1 + k2
    if len(p) == 6:
        return p[0] * np.exp(-x**2 / (2*p[1]**2)) + \
               p[2] * np.exp(-x**2 / (2*p[3]**2)) + \
               p[4] * (x+1)**(p[5])
        #return p[0] * np.exp(-x**2 / (2*p[1]**2)) + \
        #       p[2] * np.exp(-x**2 / (2*p[3]**2)) + \
        #       p[4] * np.exp(-x**2 / (2*p[5]**2))

def fitter(x, y, yerr, kg='gauss'):

    p = []
    if kg == 'gauss':
        p.append(y.max())       # amplitude 1
        p.append(0.5)           # sigma 1
        p.append(p[0] / 100.)   # amplitude 2
        p.append(p[1] * 5)      # sigma 2
        #p.append(p[0] / 1000.)  # amplitude of power law fit
        #p.append(-1)            # power law index
    if kg == 'kent':
        p.append(y.max())
        p.append(0.5)
        p.append(p[0] / 100.)   # amplitude 2
        p.append(p[1] * 5)      # sigma 2
    sigma = 1. / yerr
    popt, pcov = curve_fit(myfunc, x, y, p0=p, sigma=sigma)

    return popt, pcov


def test(config, log=False, kg='gauss'):

    if config in ['IC59','IC79']:
        nfiles = 5 if config=='IC79' else 14
        for i in range(nfiles):
            s = load_sim(config, spline=False, part=i)
            if i == 0:
                y, bins = ang_hist(s)
                x = getMids(bins)
            else:
                h = ang_hist(s)
                y += h[0]

    if config in ['IC86','IC86-II','IC86-III']:
        s = load_sim(config, spline=False)
        y, bins = ang_hist(s)
        x = getMids(bins)

    # Normalize values
    y = y.astype('float')
    ycut = (y!=0)
    x = x[ycut]
    y = y[ycut]
    dy = 1/np.sqrt(y)
    y /= y.sum()
    yerr = y * dy

    popt, pcov = fitter(x, y, yerr, kg=kg)
    print config, popt, np.sqrt(np.diagonal(pcov))
    f0 = np.array([gauss(i, popt[:2]) for i in x])
    f1 = np.array([gauss(i, popt[2:]) for i in x])
    fit = np.array([myfunc(i, *popt) for i in x])

    # Test fit for gaussian with median angular resolution as sigma
    #sig = 3.65
    #amp = 2 * 1 / (sig * np.sqrt(2*np.pi))
    #f2 = np.array([gauss(i, [amp, sig]) for i in x])

    # Get area of each gaussian
    a0, a1 = f0.sum(), f1.sum()
    atot = float(a0 + a1)
    #print a0
    #print a1
    print 'Fraction1: %.2f' % (a0 / atot)
    #print 'Unaccounted: %.2f' % (y.sum() - atot)


    fig, ax = plt.subplots()
    ax.errorbar(x, y, yerr, fmt='.')
    ax.plot(x, fit)
    #ax.plot(x, f0)
    #ax.plot(x, f1)
    #ax.plot(x, f2)
    ax.set_xlim(0, 20)
    if log:
        ax.set_yscale('log')
        ax.set_ylim(1e-6, 1)
    plt.show()

if __name__ == "__main__":

    p = argparse.ArgumentParser()
    p.add_argument('-c', '--config', dest='config',
            default='IC86',
            help='Detector configuration to run over')
    p.add_argument('-f', '--ftype', dest='ftype',
            default='gauss',
            help='Fit type [gauss | kent]')
    p.add_argument('--log', dest='log',
            default=False, action='store_true',
            help='Option to show the plot in log scale (y-axis)')
    args = p.parse_args()

    test(args.config, log=args.log, kg=args.ftype)

