#!/usr/bin/env python

###########################################################################
# Takes a list of processed hdf5 files and returns desired information    #
# as a binary dictionary for rapid reading for use with plotting          #
# Usage: save_data.py                                                     #
###########################################################################

import numpy as np
import argparse, os, time, tables, glob
from copy import deepcopy
from collections import defaultdict

from icecube import astro

import myGlobals as my
import dataFunctions as df

def saver(config, outFile, fileList):

    ##=======================================================================
    ## Starting parameters

    d = defaultdict(list)
    rDict = {'proton':'p','helium':'h','oxygen':'o','iron':'f'}
    t1 = astro.Time()

    t0 = time.time()
    for file in fileList:

        print 'Working on', file
        t = tables.openFile(file)
        q = {}

        # Get reconstructed compositions from list of children in file
        children = []
        for node in t.walk_nodes('/'):
            try: children += [node.name]
            except tables.NoSuchNodeError:
                continue
        children = list(set(children))
        compList = [n.split('_')[-1] for n in children if 'ShowerLLH_' in n]

        # Get ShowerLLH cuts and info
        rrc = t.root.ShowerLLH_proton.col('exists').astype('bool')
        for value in ['zenith', 'azimuth']:
            q[value] = t.root.ShowerLLH_proton.col(value)
        for comp in compList:
            r = rDict[comp]
            for value in ['x','y','energy']:
                q[r+'ML_'+value] = t.getNode('/ShowerLLH_'+comp).col(value)
            q[r+'LLH'] = t.getNode('/ShowerLLHParams_'+comp).col('maxLLH')

        # Timing
        mjd_day = t.root.I3EventHeader.col('time_start_mjd_day')
        mjd_sec = t.root.I3EventHeader.col('time_start_mjd_sec')
        mjd_ns  = t.root.I3EventHeader.col('time_start_mjd_ns')
        q['mjd'] = np.zeros(len(mjd_day), dtype=np.float64)
        for i in range(len(mjd_day)):
            day = int(mjd_day[i])
            sec = int(mjd_sec[i])
            ns  = int(mjd_ns[i])
            t1.SetTime(day, sec, ns)
            q['mjd'][i] = t1.GetMJD()

        # Event ID
        run = t.root.I3EventHeader.col('Run')
        event = t.root.I3EventHeader.col('Event')
        subevent = t.root.I3EventHeader.col('SubEvent')
        eventIDs = []
        for i in range(len(run)):
            eventIDs += ['%s_%s_%s' % (run[i], event[i], subevent[i])]
        q['eventIDs'] = np.asarray(eventIDs)

        # Condition and prescale passed (filter[condition, prescale])
        # For notes on weights see bottom of file
        filtermask  = df.filter_mask(config)
        filternames = df.filter_names(config)
        f = {}
        for fname in filternames:
            f[fname] = t.getNode('/'+filtermask).col(fname)
            f[fname] = f[fname][:,0].astype(float)
        filterArray = np.array([f[fname] * df.it_weights(fname) 
                for fname in f.keys()])
        filterArray[filterArray == 0] = 100.
        q['weights'] = np.amin(filterArray, axis=0)

        # Other reconstruction info
        q['NStations'] = t.root.NStations.col('value')
        q['LoudestOnEdge'] = t.root.LoudestOnEdge.col('value')
        q['Q1'] = t.root.Q1.col('value')

        # Append to existing arrays (only keep events where ShowerLLH ran)
        for key in q.keys():
            #d[key] += q[key][rrc].tolist()
            d[key] += q[key].tolist()

        t.close()

    # Laputop values
    #for key in ['x','y','zenith','azimuth','s125','e_proton','e_iron','beta']:
    #    arrays += ['lap_'+key]

    # Get Laputop info
    #for value in ['x', 'y', 'zenith', 'azimuth']:
    #    q['lap_'+value] = t.root.Laputop.col(value)
    #for value in ['s125', 'e_proton', 'e_iron', 'beta']:
    #    q['lap_'+value] = t.root.LaputopParams.col(value)

    # Convert value lists to arrays
    for key in d.keys():
        d[key] = np.asarray(d[key])

    # Get most likely composition
    rList = [rDict[comp] for comp in compList]
    full_llhs = np.array([d[r+'LLH'] for r in rList])
    max_llh = np.amax(full_llhs, axis=0)
    d['llh_comp'] = np.array(['' for i in range(len(d['pLLH']))])
    for r in rList:
        d['llh_comp'][d[r+'LLH'] == max_llh] = r
    for key in ['x', 'y', 'energy']:
        d['ML_'+key] = np.array([d[r+'ML_'+key][i]
                for i, r in enumerate(d['llh_comp'])])

    # Check for multiple most-likely compositions (mark as bad)
    badVals = np.sum(full_llhs == max_llh, axis=0)
    badVals = (badVals-1).astype('bool')
    d['llh_comp'][badVals] = ''
    for key in ['x','y','energy']:
        d['ML_'+key][badVals] = np.nan

    print "Time taken:", time.time()-t0
    print "Average time per run:", (time.time()-t0)/len(fileList)

    np.save(outFile, d)



if __name__ == "__main__":

    p = argparse.ArgumentParser(description='Converts hdf5 files to npy dict')
    p.add_argument('-c', '--config', dest='config',
            default='IT81',
            help='Detector configuration')
    p.add_argument('-d', '--date', dest='date',
            help='Month to work on [yyyymm] (optional)')
    p.add_argument('-b', '--bintype', dest='bintype',
            default='logdist',
            help='ShowerLLH binning to run over [standard|nozenith|logdist]')
    p.add_argument('--overwrite', dest='overwrite',
            default=False, action='store_true',
            help='Option to overwrite existing merged files')
    args = p.parse_args()

    my.setupShowerLLH(verbose=False)
    prefix = '%s/%s_data' % (my.llh_data, args.config)
    masterList = glob.glob('%s/files/DataLLH_????????_%s.hdf5' % 
            (prefix, args.bintype))
    masterList.sort()
    print len(masterList)

    months = [os.path.basename(f).split('_')[1][:-2] for f in masterList]
    months = sorted(list(set(months)))
    if args.date:
        months = [m for m in months if args.date in m]

    for month in months:
        fileList = glob.glob('%s/files/DataLLH_%s*_%s.hdf5' %
                (prefix, month, args.bintype))
        fileList.sort()
        print len(fileList)
        outFile = '%s/DataPlot_%s_%s.npy' % (prefix, month, args.bintype)
        if os.path.isfile(outFile) and not args.overwrite:
            print 'Outfile %s already exists. Skipping...' % outFile
            continue
        saver(args.config, outFile, fileList)


###############################################################################
## Notes on weights

"""
 - Events that pass STA8 condition have a prescale and weight of 1.
 - Events that pass STA3ii condition have a 1/2 chance to pass the STA3ii
prescale. Those that fail have a 1/3 chance to pass the STA3 prescale. So, the 
total chance of an event passing is 1/2+(1/3*1/2) = 2/3. Weight = 1.5
 - Events that pass STA3 condition but not STA3ii condition have a prescale and
weight of 3
"""
