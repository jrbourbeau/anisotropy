#!/usr/bin/env python

import glob, re, os, argparse

import myGlobals as my
import simFunctions_IT as simFunctions
from npx4.submit_npx4 import py_submit

if __name__ == "__main__":

    # Setup global path names
    my.setupGlobals(verbose=False)
    my.setupShowerLLH(verbose=False)
    simOutput = simFunctions.getSimOutput()

    p = argparse.ArgumentParser(
            description='Mass runs MakeLLHFiles.py on cluster',
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog=simOutput)
    p.add_argument('-s', '--sim', dest='sim', nargs='*',
            help='Simulation to run over')
    p.add_argument('-n', '--n', dest='n', type=int,
            default=100,
            help='Number of files to run per batch')
    p.add_argument('-b', '--bintype', dest='bintype',
            default='standard',
            help='Choose bin config for llh tables OR other desired info \
            [standard|nozenith|logdist|MC|lap|extras]')
    p.add_argument('--missing', dest='missing',
            default=False, action='store_true',
            help='Option to only submit files with missing file names')
    p.add_argument('--test', dest='test', action='store_true',
            default=False,
            help='Option for running test off cluster')
    p.add_argument('--cpp', dest='cpp',
            default=False, action='store_true',
            help='Run in C++ (static) mode')
    p.add_argument('--old', dest='old',
            default=False, action='store_true',
            help='Run using old tables')
    args = p.parse_args()

    cwd = os.getcwd()
    exList, jobIDs = [],[]

    prefix = '/data/user/fmcnally/showerllh/IT73_sim/files'

    for sim in args.sim:

        # Get config and simulation files
        config = simFunctions.sim2cfg(sim)
        files  = simFunctions.getSimFiles(sim)
        gcd    = simFunctions.getGCD(config)

        # Default parameters            
        resourcedir = my.llh_resource
        if args.old:
            resourcedir=resourcedir.replace('/resources','/original/resources')
        llhFile = '%s/LLHTables_%s.npy' % (resourcedir, args.bintype)
        gridFile = '%s/%s_grid.npy' % (resourcedir, config)
        outDir = '%s/%s_sim/files' % (my.llh_data, config)
        if args.test:
            args.n = 2

        # List of existing files to possibly check against
        existingFiles = glob.glob('%s/SimLLH_%s_%s_*.hdf5' % \
                (outDir, sim, args.bintype))
        existingFiles.sort()

        # Split into batches
        batches = [files[i:i+args.n] for i in range(0, len(files), args.n)]
        if args.test:
            batches = batches[:2]

        for bi, batch in enumerate(batches):

            # Name outfile
            start = re.split('\.', batch[0])[-3]
            end   = re.split('\.', batch[-1])[-3]
            out  = '%s/SimLLH_%s_%s' % (outDir, sim, args.bintype)
            out += '_part%s-%s.hdf5' % (start, end)

            if args.missing:
                if out in existingFiles:
                    continue

            batch.insert(0, gcd)
            batch = ' '.join(batch)

            if args.bintype == 'MC':
                cmd = 'python %s/MakeMC.py' % cwd
                arg = '--files %s -c %s -o %s' % (batch, config, out)
            elif args.bintype == 'lap':
                cmd = 'python %s/MakeLaputop.py' % cwd
                arg = '--files %s -c %s -o %s' % (batch, config, out)
            elif args.bintype == 'extras':
                cmd = 'python %s/MakeExtras.py' % cwd
                arg = '--files %s -c %s -o %s' % (batch, config, out)
            else:
                cmd = 'python %s/MakeShowerLLH.py' % cwd
                if args.cpp:
                    cmd += ' -cpp'
                    out = out.replace('.hdf5', '_cpp.hdf5')
                if args.old:
                    out = out.replace('.hdf5', '_old.hdf5')
                arg  = '--files %s -c %s -o %s' % (batch, config, out)
                arg += ' --gridFile %s --llhFile %s' % (gridFile, llhFile)

            ex  = ' '.join([cmd, arg])
            if not args.test:
                ex = ' '.join([my.env, ex])

            exList += [[ex]]
            jobIDs += ['showerllh_%s_%04i' % (sim, bi)]

    # Moderate number of submitted jobs
    njobs = len(exList)
    if njobs > 500:
        yn = raw_input('About to submit %i jobs. You sure? [y|n]: ' % njobs)
        if yn != 'y':
            raise SystemExit('Aborting...' % njobs)

    # Submit jobs
    print 'Submitting %i batches...' % njobs
    for ex, jobID in zip(exList, jobIDs):
        py_submit(ex, test=args.test, jobID=jobID)




