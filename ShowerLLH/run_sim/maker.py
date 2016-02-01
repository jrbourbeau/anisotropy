#!/usr/bin/env python

import re, os, sys, argparse, simFunctions
import myGlobals as my


if __name__ == "__main__":

    # Setup global path names
    my.setupGlobals(verbose=False)
    my.setupShowerLLH(verbose=False)
    simOutput = simFunctions.getSimOutput()

    p = argparse.ArgumentParser(
            description='Mass runs MakeLLHFiles.py on cluster',
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog=simOutput)
    p.add_argument('-s', '--sim', dest='sim',
            help='Simulation to run over')
    p.add_argument('-n', '--n', dest='n', type=int,
            help='Number of files to run per batch')
    p.add_argument('--gridFile', dest='gridFile',
            help='File containing locations for iterative grid search')
    p.add_argument('--llhFile', dest='llhFile',
            help='File with lh tables for reconstruction')
    p.add_argument('-j', '--job', dest='job', nargs='*',
            help='Run over files from specific jobs that failed')
    p.add_argument('-o', '--out', dest='out',
            help='Custom start for output filename')
    p.add_argument('--MC', dest='MC', action='store_true',
            default=False,
            help='Get MC information instead of running ShowerLLH')
    p.add_argument('--lap', dest='lap', action='store_true',
            default=False,
            help='Get Laputop information instead of running ShowerLLH')
    p.add_argument('--test', dest='test', action='store_true',
            default=False,
            help='Option for running test off cluster')
    args = p.parse_args()

    # Get config and simulation files
    config = simFunctions.sim2cfg(args.sim)
    if args.job:
        files = []
        for job in args.job:
            f      = open('/home/zgriffith/npx4/npx4-execs/npx4-'+job+'.sh')
            if config == 'IT73':
                files += [i for i in f.readlines()[-2].split(' ') if i[-3:] == 'bz2' and not 'Geo' in i]
            if config == 'IT81':
                files += [i for i in f.readlines()[-2].split(' ') if i[-2:] == 'gz' and not 'Geo' in i]
    else:
        files  = simFunctions.getSimFiles(args.sim)
    gcd    = simFunctions.getGCD(config)

    # Default parameters            
    resourcedir = my.llh_resource
    outDir = '%s/%s_sim/files' % (my.llh_data, config)
    '''
    if args.sim == '11644':
        files.remove('Level2_IC86_corsika_icetop.010042.005936.i3.gz')
        files.remove('Level2_IC86_corsika_icetop.010042.006089.i3.gz')
    '''
    if not args.n:
        eRange = simFunctions.getErange(args.sim)
        args.n = 200 if args.sim in eRange['mid'] else 100
        if args.test:
            args.n = 2
    if not args.gridFile:
        args.gridFile = '%s/%s_grid.npy' % (resourcedir, config)
    if not args.llhFile:
        args.llhFile = '%s/LLHTables_standard.npy' % resourcedir

    # Split into batches
    batches = [files[i:i+args.n] for i in range(0, len(files), args.n)]
    if args.test:
        batches = batches[:2]

    cwd = os.getcwd()
    os.chdir(my.npx4)

    for batch in batches:

        # Name outfile
        start = re.split('\.', batch[0])[-3]
        end   = re.split('\.', batch[-1])[-3]
        out = '%s/SimLLH_%s' % (outDir, args.sim)
        if args.out:
            out += '_%s' % args.out
        if args.lap:
            out += '_lap'
        if args.MC:
            out += '_MC'
        out += '_part%s-%s.hdf5' % (start, end)
        #out += '_part%s-%s.i3' % (start, end)
        print out

        batch.insert(0, gcd)
        batch = ' '.join(batch)

        if args.MC:
            cmd = 'python %s/MakeMCFiles.py' % cwd
            arg = '--files %s -o %s' % (batch, out)
        if args.lap:
            cmd = 'python %s/MakeLapFiles.py' % cwd
            arg = '--files %s -c %s -o %s' % (batch, config, out)
        else:
            cmd = 'python %s/MakeLLHFiles.py' % cwd
            arg  = '--files %s -c %s -o %s' % (batch, config, out)
            arg += ' --gridFile %s --llhFile %s' % (args.gridFile, args.llhFile)

        ex  = ' '.join([cmd, arg])
        if not args.test:
            ex = ' '.join(['./submit_npx4.sh', my.env, ex])
        os.system(ex)

    os.chdir(cwd)




