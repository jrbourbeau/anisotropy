#!/usr/bin/env python

import argparse, ROOT, tables, glob, os, re
import numpy as np
import root_numpy

import myGlobals as my
import dataFunctions as df

if __name__ == "__main__":

    # Setup global paths
    my.setupAnisotropy(verbose=False)
    my.setupShowerLLH(verbose=False)

    p = argparse.ArgumentParser(description='Creates root files from hdf5 \
            files with ShowerLLH info for use with anisotropy scripts')
    p.add_argument('-c', '--config', dest='config',
            default='IT73',
            help='Detector configuration')
    p.add_argument('-b', '--bintype', dest='bintype',
            default='logdist',
            help='ShowerLLH bin type')
    p.add_argument('--overwrite', dest='overwrite',
            default=False, action='store_true',
            help='Overwrite existing files')
    #p.add_argument('-d', '--date', dest='date',
    #        default='20100601',
    #        help='Date to run over yyyymmdd (Optional)')
    args = p.parse_args()

    hdf5Prefix = '%s/%s_data/files' % (my.llh_data, args.config)
    cfg = df.cfgconvert(args.config)
    rootPrefix = '/data/ana/CosmicRay/Anisotropy/IceTop/%s' % cfg
    outDir = '%s/%s_data/ani_files' % (my.llh_data, args.config)
    hdf5Files = glob.glob('%s/DataLLH_*_%s.hdf5' % (hdf5Prefix, args.bintype))
    hdf5Files.sort()

    for hdf5File in hdf5Files:

        # Load hdf5 file
        date = os.path.basename(hdf5File).split('_')[1]
        yy, mm, dd = date[:4], date[4:6], date[6:]

        # Check for existing outfile
        rootFile = '%s/%s-%s-%s.root' % (rootPrefix, yy, mm, dd)
        outFile = '%s/%s' % (outDir, os.path.basename(rootFile))
        if os.path.isfile(outFile) and not args.overwrite:
            print 'File %s already exists...' % outFile
            continue
        if not os.path.isfile(rootFile):
            print 'ROOT file %s does not exist. Skipping...' % rootFile
            continue

        # Extract information from hdf5 file
        print 'Working on %s...' % hdf5File
        t = tables.openFile(hdf5File)
        compDict = {'proton':'p','helium':'h','oxygen':'o','iron':'f'}
        d = {}
        # Reconstruction
        for comp in ['proton','helium','oxygen','iron']:
            print 'Working on %s...' % comp
            r = compDict[comp]
            d[r+'ML_energy'] = t.getNode('/ShowerLLH_'+comp).col('energy')
            d[r+'LLH'] = t.getNode('/ShowerLLHParams_'+comp).col('maxLLH')
        # Event header
        for key in ['Run','Event','SubEvent']:
            d[key] = t.root.ShowerLLH_proton.col(key)
        eventID = np.array(['%s_%s_%s' % 
                (d['Run'][i], d['Event'][i], d['SubEvent'][i])
                for i in range(len(d['Run']))])
        t.close()

        # Get most likely composition
        rList = ['p','h','o','f']
        full_llhs = np.array([d[r+'LLH'] for r in rList])
        max_llh = np.amax(full_llhs, axis=0)
        d['llh_comp'] = np.array(['' for i in range(len(d['pLLH']))])
        for r in rList:
            d['llh_comp'][d[r+'LLH'] == max_llh] = r
        d['ML_energy'] = np.array([d[r+'ML_energy'][i]
                for i, r in enumerate(d['llh_comp'])])

        # Check for multiple most-likely compositions (mark as bad)
        badVals = np.sum(full_llhs == max_llh, axis=0)
        badVals = (badVals-1).astype('bool')
        d['llh_comp'][badVals] = ''
        d['ML_energy'][badVals] = np.nan

        # Get cut from ROOT file
        eventID = set(eventID)
        f = ROOT.TFile(rootFile)
        t = f.NStations
        nentries = t.GetEntriesFast()
        cut = np.zeros(nentries, dtype='bool')
        for i in range(nentries):
            t.GetEntry(i)
            rootID = '%s_%s_%s' % (t.Run, t.Event, t.SubEvent)
            if rootID in eventID:
                cut[i] = True
        f.Close()

        ## WRITE TO FILE ##
        # Most likely composition
        try:
            values = np.zeros(nentries, dtype=[('comp','S1')])
            values[cut] = d['llh_comp'][:]
            root_numpy.array2root(values, outFile, 'llh_comp', 'recreate')
        except ValueError:
            print 'Length mismatch. Skipping...'
            continue

        # Most likely energy
        values = np.zeros(nentries, dtype=[('energy',float)])
        values[cut] = d['ML_energy'][:]
        root_numpy.array2root(values, outFile, 'ML_energy')

        # Likelihoods
        keys = ['pLLH','hLLH','oLLH','fLLH']
        for key in keys:
            values = np.zeros(nentries, dtype=[('llh', float)])
            values[cut] = d[key][:]
            root_numpy.array2root(values, outFile, key)




