#!/usr/bin/env python

import re, os, sys, argparse, dataFunctions
import glob
from math import ceil
import myGlobals as my
sys.path.append('/home/zgriffith/ShowerLLH/useful/')
from getGoodRuns import fileCleaner

if __name__ == "__main__":

    # Setup global path names
    my.setupGlobals(verbose=False)
    my.setupShowerLLH(verbose=False)

    p = argparse.ArgumentParser(
            description='Mass runs MakeLLHFiles.py on cluster',
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('-d', '--date', dest='date',
            help='Date of data to run over in YYYYMMDD format')
    p.add_argument('-j', '--job', dest='job', nargs='*',
            help='Run over files from specific jobs that failed')
    p.add_argument('-n', '--n', dest='n', type=int,
            help='Number of files to run per batch')
    p.add_argument('--gridFile', dest='gridFile',
            help='File containing locations for iterative grid search')
    p.add_argument('--llhFile', dest='llhFile',
            help='File with lh tables for reconstruction')
    p.add_argument('-o', '--out', dest='out',
            help='Custom start for output filename')
    p.add_argument('--test', dest='test', action='store_true',
            default=False,
            help='Option for running test off cluster')
    args = p.parse_args()

    yy = args.date[:4]
    md = args.date[4:]
    resourcedir = '/data/user/zgriffith/ShowerLLH/resources/'

    # Get list of files
    config = 'IT81-II'
    if config == 'IT81-II':
        prefix = '/data/exp/IceCube/%s/filtered/level2/' % (yy)

        if args.job:
            files = []
            for job in args.job:
                f      = open('/home/zgriffith/npx4/npx4-execs/npx4-'+job+'.sh')
                files += [i for i in f.readlines()[-2].split(' ') if i[-3:] == 'bz2' and not 'Geo' in i]
        else:
            files = glob.glob(prefix + md + '*/Level2_*Subrun*_IT.i3.bz2')
        gcdFiles = glob.glob(prefix + md + '*/Level2_*_GCD.i3.gz')
        goodRunList = resourcedir + 'burn_sample.txt'

    # Sort and clean for good runs
    files.sort()
    files = fileCleaner(goodRunList, files)

    # Default parameters            
    resourcedir = my.llh_resource
    outDir = '%s/%s_data/files' % (my.llh_data, config)
    #njobs = int(ceil(len(files)/float(args.n))) 
    if not args.n:
        args.n = 4
        njobs = int(ceil(len(files)/(float(args.n)*4))) 
        if args.test:
            args.n = 1 
            njobs  = 1 
    if not args.gridFile:
        args.gridFile = '%s/%s_grid.npy' % (resourcedir, 'IT81')
    if not args.llhFile:
        args.llhFile = '%s/LLHTables_standard.npy' % resourcedir

    cwd = os.getcwd()
    os.chdir(my.npx4)

    for i in range(njobs):
        try:
            fList = files[args.n*i:args.n*(i+1)]
        except IndexError:
            fList = files[args.n*i:]
        if args.test:
            fList = fList[:2]
        if 'GCD' not in fList[0]:
            st  = fList[0].find('Run') + 3
            run = fList[0][st:st+8]
            gcd = [gcdFile for gcdFile in gcdFiles if run in gcdFile][0]
            fList.insert(0, gcd)
        if 'GCD' in fList[-1]:
            fList.remove(fList[-1])

        # Name outfile
        outPrefix = '/data/user/zgriffith/ShowerLLH/'+config+'_data/files/burn_sample/pulse_info/'
        out = outPrefix + 'DataLLH_%s' % args.date
        for file in [fList[1], fList[-1]]:
            run_st  = file.find('Run') + 3
            subrun = 'Part' if 'Part' in file else 'Subrun'
            subrun_st = file.find(subrun) + len(subrun)
            run  = file[run_st:run_st+8]
            part = file[subrun_st:subrun_st+8]
            out += '_%s_%s' % (run, part)
        out += '.hdf5'
        print i, out

        fList = ' '.join(fList)

        cmd = 'python %s/slc_MakeLLHFiles.py' % cwd
        arg  = '--files %s -c %s -o %s' % (fList, config, out)
        arg += ' --gridFile %s --llhFile %s' % (args.gridFile, args.llhFile)

        ex  = ' '.join([cmd, arg])
        if not args.test:
            ex = ' '.join(['./submit_npx4.sh', my.env, ex])
        os.system(ex)

    os.chdir(cwd)




