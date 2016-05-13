#!/usr/bin/env python
################################################################################
# Parse a good/bad run list for IC86, get a list of DST files from the data
# warehouse, and extract an analysis-grade ntuple from each run.
#
# Then set up simple DST processing for each day by splitting the day's subruns
# into groups of files, and then merging the result.  The script will output a
# DAGMan file that can be executed by running
#
# condor_submit_dag -config dagman.config ntmaker.dag
################################################################################
#
# Version: $Id: ntmaker.py 120535 2014-06-09 21:47:04Z sybenzvi $
#
################################################################################

import os, re, subprocess, argparse
import glob
from collections import defaultdict

# Setup global path names
from goodrunlist.goodRunReader import getGoodRunFile, getGoodDates
from npx4.submit_npx4 import py_submit
import dataFunctions as df
import myGlobals as my
my.setupGlobals(verbose=False)


# -----------------------------------------------------------------------------
def mergeChunks(jobID, fileList):
    """Merge the Simple-DST chunks generated from the Level-2 ROOT files.
    """
    # Create output file
    outd = '%s' % (outp)
    outf = '%s/ic86_%s.root' % (outd, jobID)
    dest = '%s/ic86_%s.root' % (prfx, jobID)

    # Build Simple-DST generator shell script
    #fnam = ['%s/%s' % (outd, os.path.basename(f)) for f in fileList]
    fnam = fileList
    argf = ' \\\n'.join(fnam)
    ocmd = ' \\\n'.join([my.env, 'dstmerge', outf, argf])
    exeScript  = [
        '#!/bin/bash\n',
        'date',
        'hostname\n',
        #'# Copy ROOT files to local disk for I/O',
        'mkdir -p %s' % outd,
        #'cp %s %s' % (' \\\n'.join(fileList), outd),
        '\n# Process files',
        ocmd,
        '\n# Move output file to final destination',
        'mv %s %s' % (outf, dest),
        '\n# Remove chunks from remote disk',
        'rm -f %s' % ' \\\n'.join(fileList),
        '\ndate'
    ]
    scriptName = '%s/npx/exe/%s.sh' % (pdir, jobID)
    script = open(scriptName, 'w')
    for line in exeScript:
        script.write('%s\n' % line)
    script.close()
    subprocess.call(['chmod', '755', scriptName])

    # Build a condor submission script to run the shell script
    cdrScript = [
        'Universe = vanilla',
        'getenv = true',
        'Executable = %s' % scriptName,
        'Log = %s/npx/log/%s.log' % (pdir, jobID),
        'Output = %s/npx/out/%s.out' % (pdir, jobID),
        'Error = %s/npx/err/%s.err' % (pdir, jobID),
        'Notification = NEVER',
        'Queue'
    ]
    scriptName = '%s/npx/sub/%s.sub' % (pdir, jobID)
    script = open(scriptName, 'w')
    for line in cdrScript:
        script.write('%s\n' % line)
    script.close()

    return dest, scriptName


# -----------------------------------------------------------------------------
# Set up the command line options
if __name__ == '__main__':

    p = argparse.ArgumentParser(description='Simple DST ntuple maker')
    p.add_argument('-c', '--config', dest='config',
            help='Detector configuration to run over')
    p.add_argument('-n', '--nBatch', dest='n', type=int,
            default=72,     # Average run time ~ 5 hours
            help='Split sub-run processing into chunks of n files/chunk')
    p.add_argument('-d', '--date', dest='YMD',
            default=None,
            help='Force use of (yyyy[mmdd])')
    p.add_argument('--test', dest='test',
            default=False, action='store_true',
            help='Option for running off cluster to test')
    args = p.parse_args()

    #outdir = '/data/user/fmcnally/anisotropy/simple-dst'
    outdir = '/data/ana/CosmicRay/Anisotropy/IceCube/IC86/2014/simple-dst'

    # Collect dates to run over
    goodDateList = getGoodDates(args.config)
    if args.YMD:
        goodDateList = [i for i in goodDateList if i[:len(args.YMD)]==args.YMD]
    if args.test:
        goodDateList = goodDateList[:1]
        args.n = 5      # ~2 min per file, ~2000 files per day

    # Build list of executables
    exs, jobIDs = [],[]
    for yyyymmdd in goodDateList:

        files = df.getDataFiles(args.config, yyyymmdd, gcd=False)
        yy, mm, dd = yyyymmdd[:4], yyyymmdd[4:6], yyyymmdd[6:]

        # Split into batches for submission
        batches = [files[i:i+args.n] for i in range(0, len(files), args.n)]
        if args.test:
            batches = batches[0:3]

        # Process chunks for this yyyy-mm-dd
        for i, batch in enumerate(batches):
            jobID = 'simpleDST_%04s-%02s-%02s_p%02d' % (yy, mm, dd, i+1)
            files = ' '.join(batch)

            # Mimic Segev's approach. Maybe necessary for large batch sizes
            outp = "${_CONDOR_SCRATCH_DIR}/simple-dst"
            out = '%s/ic86_%s.root' % (outp, jobID)
            lfiles = ['%s/%s' % (outp, os.path.basename(f)) for f in batch]
            lfiles = ' '.join(lfiles)
            cmd = ' '.join([my.env, 'dstcut', out, lfiles])
            ex = [
                'mkdir -p %s' % outp,
                'cp %s %s' % (files, outp),  # Copy files for local I/O
                #'%s' % my.env,               # Load environment shell
                '%s' % cmd,                  # Process files
                'rm -f %s' % lfiles,         # Remove copied ROOT files
                'mv %s %s' % (out, outdir)    # Move output file
            ]
            exs += [ex]
            #out = '%s/ic86_%s.root' % (outdir, jobID)
            #exs += [[' '.join([my.env, 'dstcut', out, files])]]
            jobIDs += [jobID]

    print 'Submitting %i jobs...' % len(exs)
    if len(exs) > 500:
        raise SystemExit('Too many jobs. Use larger batches.')

    # Submit
    cwd = os.getcwd()
    os.chdir(my.npx4)
    for ex, jobID in zip(exs, jobIDs):
        py_submit(ex, jobID=jobID, test=args.test)
    os.chdir(cwd)


