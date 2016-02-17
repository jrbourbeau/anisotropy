#!/usr/bin/env python

#############################################################################
# Converts individual count tables into a single normalized LLH table
#############################################################################

import myGlobals as my
import numpy as np
import glob, re, os
import simFunctions_IT as simFunctions

def converter(fileList, outFile):

    orig, norm = {},{}

    # Sum all CountTables
    for i, file in enumerate(fileList):

        print 'Loading', file
        q = np.load(file)
        q = q.item()

        # Record and check binning
        if i == 0:
            bins = q['bins']

        # Get composition from simulation
        sim = re.split('_|\.', os.path.basename(file))[1]
        comp = simFunctions.sim2comp(sim, full=True)
        # Add to cumulative table
        if comp not in orig.keys():
            orig[comp] = np.zeros(q['counts'].shape, dtype=float)
        orig[comp] += q['counts']

    # Normalize and log tables (log10(Hist/Hist.sum()))
    print 'Normalizing tables...'
    for comp in orig.keys():
        orig[comp] += .1    # baseline so zeros don't give errors
        norm[comp] = np.zeros(orig[comp].shape, dtype=float)
        for idx in np.ndindex(orig[comp].shape[:-1]):
            norm[comp][idx] = np.log10(orig[comp][idx] / sum(orig[comp][idx]))

    # Write to file
    print 'Writing to %s' % outFile
    np.save(outFile, {'bins':bins, 'llhtables':norm})


if __name__ == "__main__":

    
    # Global variables setup for path names
    my.setupShowerLLH(verbose=False)

    masterList = glob.glob('%s/CountTables/CountTable_*.npy' % my.llh_resource)
    masterList = [f for f in masterList if '_Part' not in f]

    binList = sorted(list(set([re.split('_|\.', f)[-2] for f in masterList])))

    for bintype in binList:

        fileList = [f for f in masterList if bintype in f]
        fileList.sort()
        fileList = [f for f in fileList if '_7394_' not in f]
        fileList = [f for f in fileList if '_7483_' not in f]
        fileList = [f for f in fileList if '_7486_' not in f]
        fileList = [f for f in fileList if '_7351_' not in f]
        outFile = '%s/LLHTables_%s.npy' % (my.llh_resource, bintype)

        converter(fileList, outFile)


