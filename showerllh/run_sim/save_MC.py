#!/usr/bin/env python

############################################################################
# Saves basic MC information for finding efficiency
############################################################################

import tables, time, glob, argparse, os
import numpy as np

import myGlobals as my
import simFunctions_IT as simFunctions
from showerllh.analysis.llhtools import getEbins

if __name__ == "__main__":

    # Import global path names
    my.setupShowerLLH(verbose=False)

    p = argparse.ArgumentParser(description='Merge MC files')
    p.add_argument('-c', '--config', dest='config',
            default='IT73',
            help='Detector configuration')
    p.add_argument('--overwrite', dest='overwrite',
            default=False, action='store_true',
            help='Option to overwrite existing output file')
    args = p.parse_args()

    d = {}
    prefix = '%s/%s_sim' % (my.llh_data, args.config)
    outFile = '%s/SimPlot_MC.npy' % prefix
    if os.path.isfile(outFile) and not args.overwrite:
        raise SystemExit('Outfile %s exists.\nNothing to do...' % outFile)

    fileList = glob.glob('%s/files/SimLLH_*_MC.hdf5' % prefix)
    fileList = [f for f in fileList if '_part' not in f]
    fileList.sort()

    # Basic setup
    eDict = simFunctions.getErange()
    inv_dict = simFunctions.getErange(inverted=True)
    compList = ['P', 'He', 'O', 'Fe']
    Ebins = getEbins()

    # Create empty hists
    for comp in compList:
        d[comp] = {}
        for erange in eDict.keys():
            d[comp][erange] = np.zeros(len(Ebins)-1, dtype=int)

    # Import MC information
    for f in fileList:

        print 'Working on ' + f
        st  = f.find('SimLLH_') + 7
        sim = f[st:st+4]

        # Load the energies and composition (assumes one composition per file)
        t = tables.openFile(f)
        energies = np.log10(t.root.MCPrimary.col('energy'))
        t.close()

        # Bin in energy bins
        comp = simFunctions.sim2comp(sim)
        erange = inv_dict[sim]
        d[comp][erange] += np.histogram(energies, bins=Ebins)[0]

    np.save(outFile, d)


