#!/usr/bin/env python

###########################################################################
# Takes a list of processed hdf5 files and returns desired information    #
# as a binary dictionary for rapid reading for use with plotting          #
# Usage: save_sim.py [config]                                             #
###########################################################################

import sys, time, tables, glob, re, os
import numpy as np

def saver(config):

    d = {}
    prefix = '/data/user/fmcnally/anisotropy/sim'
    fileList = glob.glob('%s/%s_????.hdf5' % (prefix, config))
    fileList.sort()

    #typeDict = {'none':0,'P':14,'He':402,'N':1407,'O':1608,'Al':2713,'Fe':5626}
    # Type values appear to change between simulations
    typeDict = {}
    typeDict['P']  = [14, 2212]
    typeDict['He'] = [1000020040]
    typeDict['O']  = [1000080160]
    typeDict['Fe'] = [5626, 1000260560]
    # Invert dictionary
    typeDict = dict((v,k) for k in typeDict for v in typeDict[k])

    ## Keys to import ##
    flt_values  = []
    for key in ['zenith','azimuth']:
        flt_values += ['reco_'+key]
    for key in ['energy','zenith','azimuth']:
        flt_values += ['MC_'+key]

    int_values = ['nstation','MC_type','sim']
    str_values = ['comp']
    values = flt_values + int_values + str_values

    # Initialize arrays
    for value in flt_values:
        d[value] = np.array([], dtype='float')
    for value in int_values:
        d[value] = np.array([], dtype='int')
    for value in str_values:
        d[value] = np.array([], dtype='|S2')

    t0 = time.time()
    for file in fileList:

        print "Working on", file
        t = tables.openFile(file)
        q = {}

        # Get reconstruction info
        for value in ['zenith', 'azimuth']:
            q['reco_'+value] = t.root.ShowerPlane.col(value)

        # MCPrimary information
        for value in ['energy', 'zenith', 'azimuth', 'type']:
            q['MC_'+value] = t.root.MCPrimary.col(value)
        q['comp'] = np.array([typeDict[type] for type in q['MC_type']])

        # Other info
        q['nstation'] = t.root.NStations.col('value')
        sim = int(re.split('_|\.', os.path.basename(file))[1])
        q['sim'] = sim

        # append to existing arrays
        new = len(q['nstation'])
        old = len(d['nstation'])
        for value in values:
            d[value].resize(old+new)
            d[value][old:] = q[value]

        t.close()

    print "Time taken:", time.time()-t0
    print "Average time per run:", (time.time()-t0)/len(fileList)

    np.save('%s/%s_sim.npy' % (prefix, config), d)



if __name__ == "__main__":

    config = sys.argv[1]
    saver(config)
