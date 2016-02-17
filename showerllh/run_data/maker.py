#!/usr/bin/env python

import os, argparse, glob

import myGlobals as my
import dataFunctions as df
from goodrunlist.goodRunReader import getGoodDates
from npx4.submit_npx4 import py_submit

if __name__ == "__main__":

    # Setup global path names
    my.setupGlobals(verbose=False)
    my.setupShowerLLH(verbose=False)

    p = argparse.ArgumentParser(
            description='Mass runs MakeShowerLLH.py on cluster')
    p.add_argument('-c', '--config', dest='config',
            help='Detector configuration to run over')
    p.add_argument('-d', '--date', dest='date',
            help='Date to run over (yyyy[mmdd])')
    p.add_argument('-n', '--n', dest='n', type=int,
            help='Number of files to run per batch')
    p.add_argument('--old', dest='old',
            default=False, action='store_true',
            help='Run in old Python implementation')
    p.add_argument('-b', '--bintype', dest='bintype',
            default='logdist',
            #choices=['standard','nozenith','logdist','lap','extras'],
            help='Choose bin config for llh tables OR other desired info')
    p.add_argument('--missing', dest='missing',
            default=False, action='store_true',
            help='Option to only submit files with missing file names')
    p.add_argument('--test', dest='test', action='store_true',
            default=False,
            help='Option for running test off cluster')
    args = p.parse_args()

    # Default parameters
    it_geo = df.it_geo(args.config)
    llhFile = '%s/LLHTables_%s.npy' % (my.llh_resource, args.bintype)
    gridFile = '%s/%s_grid.npy' % (my.llh_resource, it_geo)
    outDir = '%s/%s_data/files' % (my.llh_data, args.config)
    outp = "${_CONDOR_SCRATCH_DIR}/simple-dst"
    if not args.n:
        if args.config == 'IT59':
            args.n = 1
        if args.config in ['IT73','IT81']:
            args.n = 20
        if args.config in ['IT81-II']:  # some jobs held, lower this
            args.n = 80
        if args.config in ['IT81-III','IT81-IV']:
            args.n = 40
        if args.test and args.config != 'IT59':
            args.n = 2

    goodDateList = getGoodDates(args.config)
    goodDates = [i for i in goodDateList if i[:len(args.date)]==args.date]
    if args.test:
        goodDates = goodDates[:2]

    cwd = os.getcwd()
    exList, jobIDs = [],[]

    for yyyymmdd in goodDates:

        # Get list of files
        files, gcdFiles = df.getDataFiles(args.config, yyyymmdd)
        ## TEMPORARY ##
        #badRuns = ['00119990','00119991']
        #files = [f for f in files if df.getRun(f) not in badRuns]
        ## END TEMPORARY ##
        runList = list(set([df.getRun(f) for f in files]))
        gcdFiles = [gcd for gcd in gcdFiles if df.getRun(gcd) in runList]
        files.sort()

        # List of existing files to possibly check against
        existingFiles = glob.glob('%s/DataLLH_%s_%s_*.hdf5' % \
                (outDir, yyyymmdd, args.bintype))
        existingFiles.sort()

        # Split into batches
        batches = [files[i:i+args.n] for i in range(0, len(files), args.n)]
        if args.test:
            batches = batches[:2]

        for bi, batch in enumerate(batches):

            # Place GCD files
            gcds, idxs = [],[]
            run_prev = '0'
            for i, file in enumerate(batch):
                run = df.getRun(file)
                if run != run_prev:
                    gcds += [gcd for gcd in gcdFiles if run in gcd]
                    idxs += [i]
                    run_prev = run
            for j, idx in enumerate(idxs[::-1]):
                batch.insert(idx, gcds[::-1][j])

            # Name outfile
            out = '%s/DataLLH_%s_%s_new' % (outp, yyyymmdd, args.bintype)
            if args.test:
                out = '%s/DataLLH_%s_%s' % (outDir, yyyymmdd, args.bintype)
            for file in [batch[1], batch[-1]]:
                run = df.getRun(file)
                part = df.getSubRun(file)
                out += '_%s_%s' % (run, part)
            out += '.hdf5'

            # Check only for missing files
            if args.missing:
                outtest = '%s/%s' % (outDir, os.path.basename(out))
                if outtest in existingFiles:
                    continue

            cmd = 'python %s/MakeShowerLLH.py' % cwd
            cmd += ' -c %s -o %s' % (args.config, out)
            cmd += ' --gridFile %s --llhFile %s' % (gridFile, llhFile)
            if not args.test:
                lfiles = ['%s/%s' % (outp, os.path.basename(f)) for f in batch]
                lfiles = ' '.join(lfiles)
                batch = ' '.join(batch)
                cmd += ' --files %s' % (lfiles)
                ex = [
                    'mkdir -p %s' % outp,           # Make directory for output
                    'cp %s %s' % (batch, outp),      # Copy files for local I/O
                    '%s %s' % (my.env, cmd),        # Process files
                    'rm -f %s' % lfiles,            # Remove copied i3 files
                    'mv %s %s' % (out, outDir)      # Move output file
                ]
            else:
                batch = ' '.join(batch)
                cmd += ' --files %s' % batch
                ex = [cmd]

            exList += [ex]
            jobIDs += ['showerllh_%s_%s_%i' % (args.config, yyyymmdd, bi)]

    if args.test:
        exList = exList[:2]
    njobs = len(exList)
    if njobs > 500:
        yn = raw_input('Submit %i jobs? [y|n]: ' % njobs)
        if yn != 'y':
            raise SystemExit('Aborting...')

    print 'Submitting %i batches...' % njobs
    for ex, jobID in zip(exList, jobIDs)[:10]:
        py_submit(ex, test=args.test, jobID=jobID)


