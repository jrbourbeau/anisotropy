#!/usr/bin/env python

from numpy import *
import healpy as hp
import matplotlib.pyplot as plt
import glob

if __name__ == "__main__":

    mapDir = '/net/user/fmcnally/ShowerLLH/maps/smoothed/'
    #masterFile = mapDir + 'IT73_24H_nocuts_20deg_relint.fits'
    #masterFile = mapDir + 'IT73_24H_nocuts_20deg_signal.fits'
    masterFile = mapDir + 'IT_24H_nocuts_20deg_signal.fits'

    m = hp.read_map(masterFile)
    minPix = m.argmin()
    stationList, minList, errList = [],[],[]

    #fileList = glob.glob(mapDir + 'IT73_24H_NStations_*20deg_relint.fits')
    fileList = glob.glob(mapDir + 'IT_24H_NStations_*20deg_relint.fits')
    fileList.sort()

    for file in fileList:

        # Extract the minimum number of stations
        st = file.find('NStations_') + 10
        nstations = file[st:st+2]
        stationList += [nstations]

        # Get the relative intensity at the minimum pixel
        m = hp.read_map(file)
        relint = m[minPix]
        minList += [relint]

        # Get the error in the relative intensity at the minimum pixel
        errFile = file.replace('relint', 'relint_err')
        m = hp.read_map(errFile)
        relint_err = m[minPix]
        errList += [relint_err]

    plt.errorbar(stationList, minList, yerr=errList, fmt='o')
    plt.title('Relative Intensity vs NStations')
    plt.xlabel('NStations')
    plt.ylabel('Relative Intensity')
    plt.xlim(-1,17)
    #plt.plot(stationList, minList, 'o')
    plt.savefig('IT_RelInt_vs_NStations.png')
    plt.show()

