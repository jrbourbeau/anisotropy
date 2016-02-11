#!/usr/bin/env python

import glob, os, sys, argparse
import simFunctions_IC as simFunctions
import myGlobals as my


if __name__ == "__main__":

    # Setup global path names
    my.setupGlobals(verbose=False)
    my.setupAnisotropy(verbose=False)

    p = argparse.ArgumentParser(
            description='Recreates DST info from simulation.')
    p.add_argument('-c', '--config', dest='config',
            help='Detector configuration [IC59 - IC86-III]')
    p.add_argument('-n', '--nBatch', dest='nBatch', type=int,
            help='Number of simulation files per submission batch')
    p.add_argument('--check', dest='check',
            default=False, action='store_true',
            help='Run over files checking for stream errors')
    p.add_argument('--test', dest='test',
            default=False, action='store_true',
            help='Option for running off cluster to test')
    p.add_argument('--b0', dest='b0', type=int,
            default=0,
            help='Starting batch number')
    args = p.parse_args()

    # Default batch sizes
    if not args.nBatch:
        if args.config == 'IC59':
            args.nBatch = 300   # ~ 1.5hr for IC59
        if args.config == 'IC79':
            args.nBatch = 30    # < 1hr for IC79
        if args.config in ['IC86','IC86-II','IC86-III']:
            # < 1000 files with ~6s per file
            args.nBatch = 50

    # Get list of files
    files   = simFunctions.getSimFiles(args.config)
    gcd     = simFunctions.getGCD(args.config)

    # Setup script
    cmdDict = {}
    cmdDict['IC59'] = 'MakeSimFiles_IC59.py'
    cmdDict['IC79'] = 'MakeSimFiles_IC79.py'
    for cfg in ['IC86','IC86-II','IC86-III']:
        cmdDict[cfg] = 'MakeSimFiles_IC86.py'
    cmdScript = 'dataCheck.py' if args.check else cmdDict[args.config]

    # Split into batches for submission
    batches = [files[i:i+args.nBatch] for i in range(0,len(files),args.nBatch)]
    if args.test:
        batches = batches[args.b0:args.b0+1]

    if len(batches) > 500:
        print '%i jobs is too many to submit. Aborting...' % len(batches)
        raise

    # Bad runs from file
    badRunFile = '%s/badFiles.txt' % my.ani_sim
    with open(badRunFile, 'r') as f:
        badFiles = f.readlines()
        badFiles = [f.strip() for f in badFiles]

    cwd = os.getcwd()
    os.chdir(my.npx4)

    print 'Submitting %i batches...' % len(batches)
    for batch in batches:

        # Necessary cleaner (temporary)
        if not args.check:
            batch = [f for f in batch if f not in badFiles]

        # Name outfile
        sim = simFunctions.cfg2sim(args.config)
        simBase = '%s/%s_%s' % (my.ani_sim, args.config, sim)
        start = batch[0].split('.')[-3]
        end   = batch[-1].split('.')[-3]
        out = '%s_%s_%s.npy' % (simBase, start, end)
        if args.check:
            out = out.replace('.npy', '_badFiles.npy')
        print out

        batch.insert(0, gcd)
        batch = ' '.join(batch)

        cmd  = 'python %s/%s' % (cwd, cmdScript)
        ex = ' '.join([cmd, args.config, out, batch])
        if not args.test:
            ex = ' '.join(['./submit_npx4.sh', my.env, ex])
        os.system(ex)

    os.chdir(cwd)

