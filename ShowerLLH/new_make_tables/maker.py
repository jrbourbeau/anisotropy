#!/usr/bin/env python

import argparse, glob, os, re
import simFunctions_IT as simFunctions
import myGlobals as my
import numpy as np

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
    os.chdir(my.npx4)

    for sim in args.sim:

        # Build fileList
        config = simFunctions.sim2cfg(sim)
        files  = simFunctions.getSimFiles(sim)
        gcd = simFunctions.getGCD(config)
        # Split into batches
        batches = [files[i:i+args.n] for i in range(0, len(files), args.n)]
        if args.test:
            batches = batches[:2]
        nbatch = len(batches)
        if nbatch > 500:
            raise SystemExit('%i jobs is too many. Aborting...' % nbatch)

        print 'Submitting %i batches...' % nbatch
        for batch in batches:

            start = re.split('\.', batch[0])[-3]
            end   = re.split('\.', batch[-1])[-3]
            outFile  = '%s/CountTable_%s_%s' % (outDir, sim, args.bintype)
            outFile += '_Part'+start+'-'+end+'.npy'

            batch.insert(0, gcd)
            batch = ' '.join(batch)

            ex   = 'python %s/MakeHist.py' % cwd
            ex  += ' -f %s -b %s -o %s' % (batch, args.bintype, outFile)
            if not args.test:
                ex = ' '.join(['./submit_npx4.sh', my.env, ex])
            os.system(ex)

    os.chdir(cwd)






