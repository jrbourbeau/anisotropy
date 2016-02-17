#!/usr/bin/env python

import argparse, glob, os
import numpy as np
import matplotlib.pyplot as plt

import myGlobals as my
from useful import getMids
from llhtools import getEbins
from goodrunlist.goodRunReader import getRunTime

def getDate(file):
    return os.path.basename(file).split('_')[1]

def dateCheck(file, date):
    fdate = getDate(file)
    if fdate[:len(date)] == date:
        return True
    return False

def getFiles(config, month=None):
    my.setupShowerLLH(verbose=False)
    histFiles = glob.glob('%s/%s_data/hists/*.npy' % (my.llh_data, config))
    if month != None:
        histFiles = [f for f in histFiles if dateCheck(f, month)]
    return sorted(histFiles)

def getValues(h):
    ykeys = [k for k in h.keys() if 'energy_w_z' in k and 'err' not in k]
    errkeys = [k for k in h.keys() if 'energy_w_z' in k and 'err' in k]
    y = np.sum([h[k] for k in ykeys], axis=0)
    yerr = np.sum([h[k] for k in errkeys], axis=0)
    return y, yerr


def monthCheck(config):

    histFiles = getFiles(config)
    dateList = [getDate(f) for f in histFiles]
    monthList = sorted(list(set([d[:6] for d in dateList])))

    # Build total histograms
    x = getMids(getEbins(reco=True))
    c0 = (x >= 6.2)
    ycfg, yerrcfg = np.zeros((2, len(x)))
    tcfg = 0.
    for f in histFiles:
        h = np.load(f)
        h = h.item()
        y, yerr = getValues(h)
        date = getDate(f)
        t = float(getRunTime(config, date=date))
        ycfg += y
        yerrcfg += yerr
        tcfg += t

    for month in monthList:
        tempFiles = [f for f in histFiles if getDate(f)[:6]==month]
        ymonth, yerrmonth = np.zeros((2, len(x)))
        tmonth = 0.
        for f in tempFiles:
            h = np.load(f)
            h = h.item()
            y, yerr = getValues(h)
            date = getDate(f)
            t = float(getRunTime(config, date=date))
            ymonth += y
            yerrmonth += yerr
            tmonth += t

        fmonth = (ymonth/tmonth)[c0]
        ferrmonth = (yerrmonth/tmonth)[c0]
        ferrmonth[ferrmonth==0] = np.inf
        fcfg = (ycfg/tcfg)[c0]
        ferrcfg = (yerrcfg/tcfg)[c0]
        ferrcfg[ferrcfg==0] = np.inf
        chi2 = np.sum((fmonth - fcfg)**2 / (ferrmonth**2 + ferrcfg**2))
        print '%s : %s' % (month, chi2)


def myplot2():

    histFiles = getFiles('IT81', '201204')
    dates = ['20120417','20120418']
    histFiles = [f for f in histFiles if any([d in f for d in dates])]

    fig, ax = plt.subplots()
    x = getMids(getEbins(reco=True))
    ytot, yerrtot = np.zeros((2, len(x)))
    ttot = 0.
    c0 = (x >= 6.2)

    for f in histFiles:

        h = np.load(f)
        h = h.item()
        y, yerr = getValues(h)
        date = getDate(f)
        t = float(getRunTime('IT81', date=date))
        ax.errorbar(x, y/t, np.sqrt(yerr)/t, label=date)

        # Increment total values
        ytot += y
        yerrtot += yerr
        ttot += t

    #ax.errorbar(x, ytot/ttot, np.sqrt(yerrtot)/ttot, label='All', fmt='k')

    ax.set_xlim((6, 9.5))
    ax.set_yscale('log')
    ax.legend()

    plt.show()        

def myplot(config, month, st, fin):

    histFiles = getFiles(config, month)

    fig, ax = plt.subplots()
    x = getMids(getEbins(reco=True))
    ytot, yerrtot = np.zeros((2, len(x)))
    ttot = 0.
    c0 = (x >= 6.2)

    for f in histFiles[st:fin]:

        h = np.load(f)
        h = h.item()
        y, yerr = getValues(h)
        date = getDate(f)
        t = float(getRunTime(config, date=date))
        ax.errorbar(x, y/t, np.sqrt(yerr)/t, label=date)

        # Increment total values
        ytot += y
        yerrtot += yerr
        ttot += t

    ax.errorbar(x, ytot/ttot, np.sqrt(yerrtot)/ttot, label='All', fmt='k')

    ax.set_xlim((6, 9.5))
    ax.set_yscale('log')
    ax.legend()

    plt.show()        


if __name__ == "__main__":

    p = argparse.ArgumentParser()
    p.add_argument('-c', '--config', dest='config',
            default='IT81',
            help='Detector configuration')
    p.add_argument('-d', '--date', dest='date',
            default='201105',
            help='Month to run over (yyyymm)')
    p.add_argument('--month', dest='month',
            default=False, action='store_true')
    p.add_argument('--plot', dest='plot',
            default=False, action='store_true')
    #p.add_argument('-r', '--range', dest='range', nargs='*', type=int)
    args = p.parse_args()

    if args.month:
        monthCheck(args.config)
    if args.plot:
        #myplot(args.config, args.date, args.range[0], args.range[1])
        myplot2()
