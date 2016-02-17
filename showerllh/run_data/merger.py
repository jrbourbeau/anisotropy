#!/usr/bin/env python

import glob, os, argparse

import myGlobals as my
from npx4.submit_npx4 import py_submit

if __name__ == "__main__":

    p = argparse.ArgumentParser(description='Merge showerllh hdf5 files')
    p.add_argument('-c', '--config', dest='config',
            default='IT81',
            help='Detector configuration')
    p.add_argument('-d', '--date', dest='date',
            help='Date to run over (optional)')
    p.add_argument('-b', '--bintype', dest='bintype',
            default='logdist',
            help='ShowerLLH binning used [standard|nozenith|logdist]')
    p.add_argument('--test', dest='test',
            default=False, action='store_true',
            help='Run off cluster to test')
    p.add_argument('--overwrite', dest='overwrite',
            default=False, action='store_true',
            help='Overwrite existing merged files')
    args = p.parse_args()

    my.setupGlobals(verbose=False)
    my.setupShowerLLH(verbose=False)

    prefix = '%s/%s_data/files' % (my.llh_data, args.config)
    hdfMerge = '%s/build/hdfwriter/resources/scripts/merge.py' % my.offline
    tempDir = '${_CONDOR_SCRATCH_DIR}/data-llh'

    masterList = glob.glob('%s/DataLLH_*_%s_0*.hdf5' % (prefix, args.bintype))
    masterList.sort()

    dates = [os.path.basename(f).split('_')[1] for f in masterList]
    dates = sorted(list(set(dates)))
    if args.date:
        dates = [ymd for ymd in dates if args.date in ymd]

    exList, jobIDs = [],[]
    for date in dates:

        # Build list of files and destination
        files = glob.glob('%s/DataLLH_%s_%s_0*.hdf5' % \
                (prefix, date, args.bintype))
        files.sort()
        tempFiles = ['%s/%s' % (tempDir, os.path.basename(f)) for f in files]
        files = ' '.join(files)
        tempFiles = ' '.join(tempFiles)
        tempOut = '%s/DataLLH_%s_%s.hdf5' % (tempDir, date, args.bintype)
        outFile = '%s/DataLLH_%s_%s.hdf5' % (prefix, date, args.bintype)
        if args.test:
            tempOut = outFile
            tempFiles = files

        # Check if file exists
        if os.path.isfile(outFile) and not args.overwrite:
            continue
        if os.path.isfile(outFile) and args.overwrite:
            print 'Outfile %s already exists. Deleting...' % outFile
            os.remove(outFile)

        cmd = 'python %s -o %s %s' % (hdfMerge, tempOut, tempFiles)
        if not args.test:
            ex = [
                'mkdir -p %s' % tempDir,
                'cp %s %s' % (files, tempDir),
                '%s %s' % (my.env, cmd),
                'rm -f %s' % tempFiles,
                'mv %s %s' % (tempOut, prefix)
            ]
        else:
            ex = [cmd]

        exList += [ex]
        jobIDs += ['llhmerge_%s_%s' % (args.config, date)]

    if args.test:
        exList = exList[:2]
    njobs = len(exList)
    if njobs > 500:
        yn = raw_input('Submit %i jobs? [y|n]: ' % njobs)
        if yn != 'y':
            raise SystemExit('Aborting...')

    print 'Submitting %i batches...' % njobs
    for ex, jobID in zip(exList, jobIDs):
        py_submit(ex, test=args.test, jobID=jobID)


