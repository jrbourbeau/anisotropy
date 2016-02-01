#!/usr/bin/env python

import os, glob, re, argparse

import myGlobals as my
import dataFunctions as df
from npx4.submit_npx4 import py_submit

# Special case for default ebins
class myAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        ebins=['4','4.25','4.5','4.75','5','5.25','5.5','5.75','6','6.5','100']
        if not values:
            setattr(namespace, self.dest, ebins)
        else:
            setattr(namespace, self.dest, values)

# Special case for default sbins
class myActionS(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        sbins = ['-2','-1','-0.75','-0.5','-0.25','1','4']
        if not values:
            setattr(namespace, self.dest, sbins)
        else:
            setattr(namespace, self.dest, values)

# Special case for default sbins
class myActionC(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        cbins = ['p','h','o','f']
        if not values:
            setattr(namespace, self.dest, cbins)
        else:
            setattr(namespace, self.dest, values)


if __name__ == "__main__":

    # Setup global variables
    my.setupGlobals(verbose=False)
    my.setupAnisotropy(verbose=False)

    p = argparse.ArgumentParser(
            description='Makes daily healpix maps using time-scrambling code')

    # General parameters
    p.add_argument('-c', '--config', dest='config',
            help='Detector configuration [IC59,IT73,IC86-III...]')
    p.add_argument('-n', '--nBatch', dest='nBatch', type=int,
            help='Number of files to run per submission batch')
    p.add_argument('-t', '--integrationtime', dest='integrationtime',
            default='24',
            help='Integration time (in hours)')
    p.add_argument('-m', '--method', dest='method',
            default='sid',
            help='Time transformation method [sid|anti|ext|solar]')
    # IceTop specific
    p.add_argument('-f', '--filter', dest='filter',
            help='IceTop filter to apply [STA3|STA8|NotSTA8]')
    p.add_argument('--sbins', dest='sbins',
            default=False, action=myActionS, nargs='*',
            help='Option to do all s125 bins')
    p.add_argument('--comp', dest='comp',
            default=False, action=myActionC, nargs='*',
            help='Option to do all composition maps')
    p.add_argument('--emin', dest='emin',
            help='Optional minimum reconstructed energy value')
    # Energy cut options
    p.add_argument('--ebins', dest='ebins',
            default=False, action=myAction, nargs='*',
            help='Option to do all energy bins')
    # Additional options
    p.add_argument('-o', '--outdir', dest='outdir',
            default='%s/maps' % my.ani_data,
            help='Destination directory for output files')
    p.add_argument('--test', dest='test',
            default=False, action='store_true',
            help='Option for running off cluster to test')
    p.add_argument('--b0', dest='b0', type=int,
            default=0,
            help='Starting batch number')
    p.add_argument('--batchfile', dest='batchfile',
            help='Name for output file containing file lists')
    p.add_argument('--overwrite', dest='overwrite',
            default=False, action='store_true',
            help='Option to overwrite existing map files')

    args = p.parse_args()
    c_opts = vars(args).copy()
    if not args.config:
        p.error('Detector configuration not given')

    # Defaults for submission batch size
    nBatchDict = {}
    for config in ['IC59','IC79','IC86','IC86-II','IC86-III','IC86-IV']:
        nBatchDict[config] = 1
    for config in ['IT59','IT73','IT81','IT81-II','IT81-III','IT81-IV']:
        nBatchDict[config] = 10
    if args.test and not args.nBatch:
        args.nBatch = 1
    nBatch = nBatchDict[args.config] if not args.nBatch else args.nBatch

    # Default spline files
    prefix = my.ani_sim
    simDict = {'IC59':'4046','IC79':'6451'}
    for config in ['IC86','IC86-II','IC86-III','IC86-IV']:
        simDict[config] = '10649'
    if args.ebins and args.config[:2]=='IC':
        sim = simDict[args.config]
        c_opts['spline'] = '%s/%s_%s_spline.fits' % (prefix, args.config, sim)
        ## NOTE: Simulation not yet available for newest year. Use IC86-III
        if args.config == 'IC86-IV':
            c_opts['spline'] = '%s/IC86-III_%s_spline.fits' % (prefix, sim)

    # Base name for outfiles
    out = '%s_%02dH_%s' % (args.config, int(args.integrationtime), args.method)
    outBase = '%s/%s/%s' % (args.outdir, args.config, out)
    if args.filter:
        outBase += '_%s' % args.filter
    if args.emin:
        outBase += '_emin'
    c_opts['outBase'] = outBase

    # Base name for batch file
    batchPrefix = '%s/timeScramble/tempFiles' % my.ani_home
    if not args.batchfile:
        c_opts['batchfile'] = '%s/%s.txt' % \
                (batchPrefix, os.path.basename(outBase))

    # Remove dates that have already been run
    fileList = []
    masterList = df.getDSTfiles(args.config)
    for rootFile in masterList:
        date = re.split('/|_|\.', rootFile)[-2]
        testFiles = [outBase]
        if args.ebins:
            testFiles = ['%s_%s-%sGeV' % (f, args.ebins[i], args.ebins[i+1]) \
                    for i in range(len(args.ebins)-1) for f in testFiles]
        if args.sbins:
            testFiles = ['%s_%sto%ss125' % (f, args.sbins[i], args.sbins[i+1]) \
                    for i in range(len(args.sbins)-1) for f in testFiles]
        if args.comp:
            testFiles = ['%s_%s' % (f, comp) \
                    for comp in args.comp for f in testFiles]
        testFiles = ['%s_%s.fits' % (f, date) for f in testFiles]
        if all([os.path.isfile(testFile) for testFile in testFiles]):
            if not args.overwrite:
                continue
            for testFile in testFiles:
                os.remove(testFile)
        fileList.append(rootFile)

    # Split into batches for submission
    batchList = [fileList[i:i+nBatch] for i in range(0, len(fileList), nBatch)]
    if args.test:
        batchList = batchList[args.b0:args.b0+1]
    batchList = [' '.join(batch)+'\n' for batch in batchList]
    with open(c_opts['batchfile'], 'w') as f:
        f.writelines(batchList)

    # Setup parameters to feed C++ script
    print 'Parameters for submission:'
    validArgs  = ['config','integrationtime','method','filter','outBase']
    validArgs += ['spline','batchfile','emin']
    c_opts = {key:c_opts[key] for key in c_opts if c_opts[key]!=None and \
            key in validArgs}
    for key in sorted(c_opts.keys()):
        print '  --%s %s' % (key, c_opts[key])
    c_opts = [['--'+key, c_opts[key]] for key in sorted(c_opts.keys())]
    # Python black magic for flattening list of arbitrary depth
    flatList = lambda *n: (e for a in n \
            for e in (flatList(*a) if isinstance(a, (tuple, list)) else (a,)))
    c_opts = list(flatList(c_opts))
    c_opts = ' '.join(c_opts)
    # Ebins/Sbins need to go last (multitoken parameters)
    if args.ebins:
        print '  --ebins ' + ' '.join(args.ebins)
        c_opts += ' --ebins %s' % ' '.join(args.ebins)
    if args.sbins:
        print '  --sbins ' + ' '.join(args.sbins)
        c_opts += ' --sbins %s' % ' '.join(args.sbins)
    if args.comp:
        print '  --comp ' + ' '.join(args.comp)
        c_opts += ' --comp %s' % ' '.join(args.comp)

    # Increase requested memory for jobs that produce multiple maps
    sublines = None
    if args.ebins or args.sbins or args.comp:
        sublines = ["request_memory = 5000"]

    # Submit files
    cwd = os.getcwd()
    os.chdir(my.npx4)

    print 'Submitting %i files...' % len(batchList)
    for i in range(len(batchList)):

        jobID = 'timescramble_%s_%05i' % (args.config, i)
        cmd  = '%s/TimeScramble' % cwd
        ex = [cmd, '--batch_idx', str(i), c_opts]
        if not args.test:
            ex = [my.env] + ex
        ex = [' '.join(ex)]
        py_submit(ex, sublines=sublines, test=args.test, jobID=jobID)

    os.chdir(cwd)


