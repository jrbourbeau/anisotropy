#!/usr/bin/env python

############################################################################
# Saves basic MC information for finding efficiency
############################################################################

import tables, time, glob
from numpy import *

if __name__ == "__main__":

    d = {}
    config = 'IT73'
    prefix = '/net/user/fmcnally/ShowerLLH/%s_sim/' % config
    resourcedir = '/net/user/fmcnally/ShowerLLH/resources/'
    fileList = glob.glob(prefix + 'files/SimLLH_????_MC.hdf5')
    fileList.sort()

    # Input/Output directories
    typeDict = {0:'none',14:'P',402:'He',1407:'N',1608:'O',2713:'Al',5626:'Fe'}
    if config == 'IT73':
        eDict = {}
        eDict['low']  = ['7351','7394','7483','7486']
        eDict['med']  = ['7006','7007','7241','7242','7262','7263']
        eDict['high'] = ['7579','7784','7791','7851']
    #tList = ['P', 'He', 'O', 'Fe']
    tList = ['P', 'Fe']

    # Setup dictionary for storage
    for t in tList:
        d[t] = {}

    # Import Ebins
    binDict = load(resourcedir + 'ShowerLLH_bins.npy')
    binDict = binDict.item()
    Ebins = binDict['Ebins']

    # Create empty hists
    for t in tList:
        for erange in eDict.keys():
            d[t][erange] = zeros(len(Ebins)-1, dtype=int)
            for i in range(3,5):
                for j in range(i):
                    k = str(i) + str(j)
                    d[t][erange+'_'+k] = zeros(len(Ebins)-1, dtype=int)

    # Import MC information
    for fl in fileList:

        print 'Working on ' + fl
        st  = fl.find('SimLLH_') + 7
        sim = fl[st:st+4]

        # Load the energies and composition (assumes one composition per file)
        t = tables.openFile(fl)
        energies = log10(t.root.MCPrimary.col('energy'))
        cosz = cos(t.root.MCPrimary.col('zenith'))
        comp = typeDict[t.root.MCPrimary.col('type')[0]]
        t.close()

        # Get the energy range for the simulation
        for key in eDict.keys():
            if sim in eDict[key]:
                erange = key

        # Bin in energy bins
        d[comp][erange] += histogram(energies, bins=Ebins)[0]

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
            d[comp][erange+'_'+zbin] += histogram(energies[c0], bins=Ebins)[0]

    outFile = prefix + config+'_MC.npy'
    save(outFile, d)

