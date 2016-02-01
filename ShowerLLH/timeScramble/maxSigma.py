#!/usr/bin/env python

from numpy import *
import matplotlib.pyplot as plt
import healpy as hp
import glob

def plotSigmas(params):

    prefix = '/net/user/fmcnally/ShowerLLH/maps/smoothed/'
    masterList = glob.glob(prefix + '*signal.fits')
    masterList.sort()
    fileList = [f for f in masterList if all(x in f for x in params)]
    degList, minList, maxList = [],[],[]

    for file in fileList:

        # Get min and max
        map = hp.read_map(file)
        minList.append(map.min())
        maxList.append(map.max())

        # Get degree of smoothing
        end = file.find('deg')
        st  = file.rfind('_', 0, end) + 1
        deg = int(file[st:end])
        degList.append(deg)

    plt.plot(degList, minList, 'r.')
    plt.plot(degList, maxList, 'b.')
    plt.show()


if __name__ == "__main__":

    params = ['IT73', '24H']
    plotSigmas(params)






