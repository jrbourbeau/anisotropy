#!/usr/bin/env python

##===========================================================================##
## Runs save_hist.py over cluster
## NOTE: Submits daily files because IT81-III is SLOW (need to figure out why)
##===========================================================================##

import argparse, re, glob, os

import myGlobals as my
from npx4.submit_npx4 import py_submit

if __name__ == "__main__":

    p = argparse.ArgumentParser(description='Converts hdf5 files to npy dict')
    p.add_argument('-c', '--configs', dest='configs', nargs='*',
            default=['IT73','IT81','IT81-II','IT81-III','IT81-IV'],
            help='Detector configuration [IT73 --> IT81-IV]')
    p.add_argument('-b', '--bintype', dest='bintype',
            default='logdist',
            help='ShowerLLH binning to run over [standard|nozenith|logdist]')
    p.add_argument('--sky', dest='sky',
            default=False, action='store_true',
            help='Write the sky histograms')
    p.add_argument('--test', dest='test',
            default=False, action='store_true',
            help='Option for running test off cluster')
    p.add_argument('--overwrite', dest='overwrite',
            default=False, action='store_true',
            help='Option to overwrite existing files')
    args = p.parse_args()

    # Setup global paths
    my.setupGlobals(verbose=False)
    my.setupShowerLLH(verbose=False)

    cwd = os.getcwd()
    exList, jobIDs = [],[]

    for config in args.configs:

        hdf5dir = '%s/%s_data/files' % (my.llh_data, config)
        histdir = '%s/%s_data/hists' % (my.llh_data, config)
        files = glob.glob('%s/DataLLH_*_%s.hdf5' % (hdf5dir, args.bintype))
        files.sort()

        for fidx, f in enumerate(files):

            # Check for existing outfile
            base = re.split('/|\.', f)[-2]
            outFile = '%s/%s.npy' % (histdir, base)
            if args.sky:
                outFile = outFile.replace('.npy', '_sky.npy')
            if os.path.isfile(outFile) and not args.overwrite:
                continue

            cmd = 'python %s/save_hist.py' % cwd
            cmd += ' -c %s -f %s -o %s' % (config, f, outFile)
            if args.sky:
                cmd += ' --sky'
            if not args.test:
                cmd = '%s %s' % (my.env, cmd)

            exList += [[cmd]]
            jobIDs += ['histllh_%s_%i' % (config, fidx)]

    if args.test:
        exList = exList[:2]
        jobIDs = jobIDs[:2]

    njobs = len(exList)
    if njobs > 500:
        yn = raw_input('Submit %i jobs? [y|n]: ' % njobs)
        if yn != 'y':
            raise SystemExit('Aborting...')

    print 'Submitting %i batches...' % njobs
    for ex, jobID in zip(exList, jobIDs):
        py_submit(ex, test=args.test, jobID=jobID)

