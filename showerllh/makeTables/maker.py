#!/usr/bin/env python

import argparse, glob, os, re
import numpy as np

import myGlobals as my
import simFunctions_IT as simFunctions
from npx4.submit_npx4 import py_submit

if __name__ == "__main__":

    # Global variables setup for path names
    my.setupGlobals(verbose=False)
    my.setupShowerLLH(verbose=False)
    simOutput = simFunctions.getSimOutput()

    p = argparse.ArgumentParser(
            description='Makes binned histograms for use with ShowerLLH',
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog=simOutput)
    p.add_argument('-s', '--sim', dest='sim', nargs='*',
            help='Simulation dataset to run over')
    p.add_argument('-n', '--n', dest='n', type=int,
            default=100,
            help='Number for files to run per submission batch')
    p.add_argument('-b', '--bintype', dest='bintype',
            default='standard',
            choices=['standard','nozenith','logdist'],
            help='Option for a variety of preset bin values')
    p.add_argument('--test', dest='test', action='store_true',
            default=False,
            help='Option for running off cluster to test')
    args = p.parse_args()

    # Default parameters
    binFile = '%s/ShowerLLH_bins.npy' % my.llh_resource
    outDir = '%s/CountTables' % my.llh_resource
    if args.test and args.n==100:
        args.n = 1

    cwd = os.getcwd()
    exList = []

    for sim in args.sim:

        # Build fileList
        config = simFunctions.sim2cfg(sim)
        files  = simFunctions.getSimFiles(sim)
        gcd = simFunctions.getGCD(config)
        # Split into batches
        batches = [files[i:i+args.n] for i in range(0, len(files), args.n)]
        if args.test:
            batches = batches[:2]

        for batch in batches:

            start = re.split('\.', batch[0])[-3]
            end   = re.split('\.', batch[-1])[-3]
            outFile  = '%s/CountTable_%s_%s' % (outDir, sim, args.bintype)
            outFile += '_Part'+start+'-'+end+'.npy'

            batch.insert(0, gcd)
            batch = ' '.join(batch)

            cmd = 'python %s/MakeHist.py' % cwd
            ex  = '%s -f %s -b %s -o %s' % (cmd, batch, args.bintype, outFile)
            if not args.test:
                ex = ' '.join([my.env, ex])
            exList += [[ex]]

    # Moderate number of submitted jobs
    njobs = len(exList)
    if njobs > 500:
        yn = raw_input('About to submit %i jobs. You sure? [y|n]: ' % njobs)
        if yn != 'y':
            raise SystemExit('Aborting...' % njobs)

    # Submit jobs
    print 'Submitting %i batches...' % njobs
    for ex in exList:
        py_submit(ex, test=args.test)






