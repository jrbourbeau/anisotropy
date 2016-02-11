#!/usr/bin/env python

import glob, os, re, sys, argparse
import simFunctions_IT as simFunctions
import myGlobals as my

def cleanBadRuns(config, fileList):

    badRunFile = 'badFiles.npy'
    badFiles = np.load(badRunFile)
    badFiles = badFiles.item()

    if config not in badFiles.keys():
        return fileList
    fileList = [f for f in fileList if f not in badFiles[config]]
    return fileList


if __name__ == "__main__":

    # Setup global path names
    my.setupGlobals(verbose=False)
    simOutput = simFunctions.getSimOutput()

    p = argparse.ArgumentParser(
            description='Creates npy dictionary files from sim for analysis',
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog=simOutput)
    p.add_argument('-s', '--sim', dest='sim',
            help='Simulation to run over')
    p.add_argument('-n', '--nBatch', dest='nBatch', type=int,
            default=200,
            help='Number of simulation files per submission batch')
    p.add_argument('--test', dest='test',
            default=False, action='store_true',
            help='Option for running off cluster to test')
    p.add_argument('--b0', dest='b0', type=int,
            default=0,
            help='Starting batch number')
    p.add_argument('--part0', dest='part0', type=int,
            help='Optional part of simulation to start at')
    p.add_argument('--part1', dest='part1', type=int,
            help='Optional part of simulation to end at')
    args = p.parse_args()

    # Get config and simulation files
    config = simFunctions.sim2cfg(args.sim)
    files  = simFunctions.getSimFiles(args.sim)
    gcd    = simFunctions.getGCD(config)

    # Limit to desired part ranges (optional)
    getPart = lambda f: int(re.split('_|\.', f)[-3])
    if args.part0:
        files = [f for f in files if getPart(f) >= args.part0]
    if args.part1:
        files = [f for f in files if getPart(f) <= args.part1]

    # Eliminate bad files (hopefully temporary)
    with open('badFiles.txt','r') as f:
        badFiles = f.readlines()
        badFiles = [l.strip() for l in badFiles]
    files = [f for f in files if f not in badFiles]

    # Split into batches for submission
    batches = [files[i:i+args.nBatch] for i in range(0,len(files),args.nBatch)]
    if args.test:
        batches = batches[args.b0:args.b0+1]
    if len(batches) > 500:
        print '%i jobs is too many to submit. Aborting...' % len(batches)
        sys.exit()

    cwd = os.getcwd()
    os.chdir(my.npx4)

    print 'Submitting %i batches...' % len(batches)
    for batch in batches:

        # Name outfile
        outPrefix = '/data/user/fmcnally/anisotropy/sim'
        start = batch[0].split('.')[-3]
        end   = batch[-1].split('.')[-3]
        out = '%s/%s_%s_%s_%s.hdf5' % (outPrefix, config, args.sim, start, end)
        print out

        batch.insert(0, gcd)
        batch = ' '.join(batch)

        cmdScript = 'MakeSimFiles.py'
        cmd  = 'python %s/%s' % (cwd, cmdScript)
        ex = ' '.join([cmd, config, out, batch])
        if not args.test:
            ex = ' '.join(['./submit_npx4.sh', my.env, ex])
        os.system(ex)

    os.chdir(cwd)

