#!/usr/bin/env python

from numpy import *
import os, glob
import healpy as hp
import matplotlib.pyplot as plt

def savePlot():

    prefix = '/net/user/fmcnally/ShowerLLH/maps/smoothed/'
    fileList = glob.glob(prefix + 'IT??_*_24H*nocuts*20deg_relint.fits')
    fileList.sort()
    for f in fileList:
        print 'Working on', f
        test = hp.read_map(f)
        outFile = os.path.basename(f).replace('.fits','.png')
        title = outFile[:-4]
        hp.mollview(test, rot=-180, min=-.002, max=.002, title=title)
        plt.savefig(outFile)
        plt.close()

if __name__ == "__main__":

    savePlot()
