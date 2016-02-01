#!/usr/bin/env python

############################################################################
# Saves basic MC information for finding efficiency
############################################################################

import tables, pickle, time, sys
from numpy import *

if __name__ == "__main__":

    if len(sys.argv) != 2:
        print 'Incorrect usage: Save_MC.py [IT59|IT73]'
        raise 
    config = sys.argv[1]

    # Input/Output directories
    prefix = '/net/user/fmcnally/ShowerLLH/%s_sim/' % config
    resourcedir = '/net/user/fmcnally/ShowerLLH/resources/'
    if config == 'IT73':
        simDict = {}
        simDict['P'] = ['7351','7006','7579']
        simDict['H'] = ['7483','7241','7263','7791']
        simDict['O'] = ['7486','7242','7262','7851']
        simDict['F'] = ['7394','7007','7784']
        eDict = {}
        eDict['low']  = ['7351','7394','7483','7486']
        eDict['med']  = ['7006','7007','7241','7242','7262','7263']
        eDict['high'] = ['7579','7784','7791','7851']

    # Setup dictionary for storage
    d = {}
    for key in simDict.keys():
        d[key] = {}

    # Import Ebins
    t = tables.openFile(resourcedir + 'ShowerLLH_bins.hdf5')
    Ebins = t.root.bins.col('Ebins')[0]
    t.close()

    # Create empty hists
    for e in simDict.keys():
        for erange in eDict.keys():
            d[e][erange] = zeros(len(Ebins)-1, dtype=int)
            for i in range(3,5):
                for j in range(i):
                    k = str(i) + str(j)
                    d[e][erange+'_'+k] = zeros(len(Ebins)-1, dtype=int)

    # Import MC information
    for e in simDict.keys():

        simList = simDict[e]        
        print 'Working on ' + e

        for sim in simList:
            # Load the energies
            t = tables.openFile(prefix + 'files/SimLLH_'+sim+'_MC.hdf5')
            energies = log10(t.root.MCPrimary.col('energy'))
            cosz = cos(t.root.MCPrimary.col('zenith'))
            t.close()

            # Get the energy range for the simulation
            for key in eDict.keys():
                if sim in eDict[key]:
                    erange = key

            # Bin in energy bins
            d[e][erange] += histogram(energies, bins=Ebins)[0]
            # Additionally bin in zenith bins
            zcuts = {}
            zcuts['40'] = (cosz >= 0.95)
            zcuts['41'] = (cosz >= 0.90) * (cosz < 0.95)
            zcuts['42'] = (cosz >= 0.85) * (cosz < 0.90)
            zcuts['43'] = (cosz >= 0.80) * (cosz < 0.85)
            zcuts['30'] = (cosz >= 0.93)
            zcuts['31'] = (cosz >= 0.87) * (cosz < 0.93)
            zcuts['32'] = (cosz >= 0.80) * (cosz < 0.87)
            for zbin in zcuts.keys():
                c0 = zcuts[zbin]
                d[e][erange+'_'+zbin] += histogram(energies[c0], bins=Ebins)[0]

    fl = open(prefix + '%s_MC.pkl' % config, 'wb')
    pickle.dump(d, fl)
    fl.close()

