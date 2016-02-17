#!/usr/bin/env python

import glob, os, load_data
import matplotlib.pyplot as plt
from numpy import *
from llhtools import *

def makeCountHist(config, yyyymm):

    d = load_data.load_data(config, yyyymm)
    r = log10(d['ML_energy'])
    w = d['weights']
    cut = d['cuts']['llh']
    Emids = getEmids()
    eList = ['p', 'f']

    # Get the counts and error
    h = {}
    h['N'], h['err'] = {},{}
    h['N']['All']   = Nfinder(r, cut, w=w)
    h['err']['All'] = Efinder(r, cut, w=w)
    for e in eList:
        h['N'][e]   = Nfinder(r, cut*d[e], w=w)
        h['err'][e] = Efinder(r, cut*d[e], w=w)

    # Save to file
    outFile = 'hists/CountHists_%s_%s.npy' % (config, yyyymm)
    save(outFile, h)

def plotCountHist():

    fileList = glob.glob('hists/CountHists_????_??????.npy')
    fileList.sort()

    for i in range(len(fileList)):
        temp = load(fileList[i])
        temp = temp.item()
        if i == 0:
            joint = temp
        else:
            for key in joint.keys():
                for key2 in joint[key].keys():
                    joint[key][key2] += temp[key][key2]

    Emids = getEmids()
    plt.errorbar(Emids, joint['N']['f'], yerr=joint['err']['f'], fmt='r.')
    plt.errorbar(Emids, joint['N']['p'], yerr=joint['err']['p'], fmt='b.')
    plt.errorbar(Emids, joint['N']['All'], yerr=joint['err']['All'], fmt='k.')
    plt.yscale('log')
    plt.show()


if __name__ == "__main__":

    configList = ['IT73','IT81']
    for config in configList:

        prefix = '/net/user/fmcnally/ShowerLLH/'+config+'_data/'
        fileList = glob.glob(prefix + 'DataPlot_20*.npy')
        ymList = [os.path.basename(f)[9:15] for f in fileList]
        ymList = list(set(ymList))      # Remove duplicates
        ymList.sort()

        for yyyymm in ymList:
            print 'Working on %s %s...' % (config, yyyymm)
            makeCountHist(config, yyyymm)


