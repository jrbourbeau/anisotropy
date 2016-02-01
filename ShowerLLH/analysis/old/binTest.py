#!/usr/bin/env python

from numpy import *
import matplotlib.pyplot as plt

def plotter(a, abins, log=False):

    temp = a.flatten()
    h, bins = histogram(temp, bins=abins)
    amids = bin2mid(abins)
    print len(h[h==0]), 'zeros in bins'
    plt.plot(amids, h, '.')
    plt.ylabel('log10(counts)')
    plt.xlabel('log10(charge(VEM))')
    #plt.xscale('log')
    if log:
        plt.yscale('log')
    plt.show()

def bin2mid(bins):
    mids = (bins[:-1] + bins[1:])/2.
    if inf in bins:
        mids[-1] = mids[-2] + (mids[-2]-mids[-3])
    return mids

def binFinder(a, nbins):

    temp = a.flatten()
    temp.sort()
    temp = temp[temp!=0]
    num = len(temp) / nbins
    bins = [temp[i*num] for i in range(nbins)]
    return asarray(bins)


def binTest(cs, Cbins, nbins=False):

    # Non-zero plot for charge
    if not nbins:
        nbins = len(Cbins)
    cflat = cs.flatten()
    h1, b1 = histogram(cflat[cflat!=0], bins=Cbins[1:])
    Cbin2 =  binFinder(cflat, nbins+1)
    h2, b2 = histogram(cflat[cflat!=0], bins=Cbin2)

    print 'Original bins:'
    print Cbins
    print 'New bins:'
    print Cbin2
    Cmids = bin2mid(Cbins[1:])
    Cmid2 = bin2mid(Cbin2)

    plt.plot(Cmids, h1, 'gx', label='old')
    #plt.plot(Cmid2, h2, 'rx', label='new')
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig('chargebins.png')
    plt.legend()

    plt.show()


def eLook(en, cs, Cbins):

    ecuts = [5, 6, 7, 8]
    index = digitize(en, ecuts) - 1
    Cbins2 =  binFinder(cs.flatten(), 30)
    for i in range(len(ecuts)-1):
        cut = (index==i)
        cflat = (cs[cut]).flatten()
        h1, b1 = histogram(cflat[cflat!=0], bins=Cbins[1:])
        print h1
        h2, b2 = histogram(cflat[cflat!=0], bins=Cbins2)

        plt.plot(log10(Cbins[1:-1]), log10(h1), '.', label='old'+str(i))
        plt.plot(log10(Cbins2[1:]), log10(h2), 'x', label='new'+str(i))

    plt.legend()
    plt.show()

