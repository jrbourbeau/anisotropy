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
			default='Sidereal',
			help='Time transformation method [Sidereal, Anti, Solar, Extended]')
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
			default='{}/maps'.format(my.ani_data),
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
	#print('args = {}'.format(args))
	#print('vars(args) = {}'.format(vars(args)))
	c_opts = vars(args).copy()
	#print('c_opts = {}'.format(c_opts))
	if not args.config:
		p.error('Detector configuration not given')

	# Defaults for submission batch size (IT=10 and IC=1)
	nBatchDict = {}
	for config in ['IC59','IC79','IC86-2011','IC86-2012','IC86-2013','IC86-2014']:
	#for config in ['IC59','IC79','IC86','IC86-II','IC86-III','IC86-IV']:
		nBatchDict[config] = 1
	for config in ['IT59','IT73','IT81-2011','IT81-2012','IT81-2013','IT81-2014']:
	#for config in ['IT59','IT73','IT81','IT81-II','IT81-III','IT81-IV']:
		nBatchDict[config] = 10
	if args.test and not args.nBatch:
		args.nBatch = 1
	if args.config in nBatchDict.keys():
		nBatch = nBatchDict[args.config] if not args.nBatch else args.nBatch
	else:
		print('An invalid detector configuration was entered. Exiting program.')
		raise SystemExit

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
	out = '{}_{:02d}H_{}'.format(args.config, int(args.integrationtime), args.method)
	outbase = '{}/{}/{}'.format(args.outdir, args.config, out)
	if args.filter:
		outbase += '_{}'.format(args.filter)
	if args.emin:
		outbase += '_emin'
	c_opts['outbase'] = outbase
	print('outbase = {}'.format(c_opts['outbase']))

	# Base name for batch file
	batchPrefix = '{}/timescramble/tempfiles'.format(my.ani_home)
	if not args.batchfile:
		c_opts['batchfile'] = '{}/{}.txt'.format(batchPrefix, os.path.basename(outbase))
	print('batchfile = {}'.format(c_opts['batchfile']))

	# Remove dates that have already been run
	fileList = []
	masterList = df.getDSTfiles(args.config)
	#print('masterList = {}'.format(masterList))
	for rootFile in masterList:
		date = re.split('/|_|\.', rootFile)[-2]
		testFiles = [outbase]
		if args.ebins:
			testFiles = ['%s_%s-%sGeV' % (f, args.ebins[i], args.ebins[i+1]) \
					for i in range(len(args.ebins)-1) for f in testFiles]
		if args.sbins:
			testFiles = ['%s_%sto%ss125' % (f, args.sbins[i], args.sbins[i+1]) \
					for i in range(len(args.sbins)-1) for f in testFiles]
		if args.comp:
			testFiles = ['%s_%s' % (f, comp) \
					for comp in args.comp for f in testFiles]
		testFiles = ['{}_{}.fits'.format(f, date) for f in testFiles]
		#print('testFiles = {}'.format(testFiles))
		if all([os.path.isfile(testFile) for testFile in testFiles]):
			if not args.overwrite:
				continue
			for testFile in testFiles:
				os.remove(testFile)
		fileList.append(rootFile)
	print('fileList = {}'.format(fileList))
    
	# Split into batches for submission
	batchList = [fileList[i:i+nBatch] for i in range(0, len(fileList), nBatch)]
	#print('batchList = {}'.format(batchList))
	if args.test:
		batchList = batchList[args.b0:args.b0+1]
	batchList = [' '.join(batch)+'\n' for batch in batchList]
	#print('batchList = {}'.format(batchList))
	with open(c_opts['batchfile'], 'w') as f:
		f.writelines(batchList)
    
	# Setup parameters to feed C++ script
	print('Parameters for submission:')
	validArgs  = ['config','integrationtime','method','filter','outbase']
	validArgs += ['spline','batchfile','emin']
	c_opts = {key:c_opts[key] for key in c_opts if c_opts[key]!=None and key in validArgs}
	for key in sorted(c_opts.keys()):
		print(  '--{} {}'.format(key, c_opts[key]))
	c_opts = [['--'+key, c_opts[key]] for key in sorted(c_opts.keys())]
	#print('c_opts = {}'.format(c_opts))
	# Python black magic for flattening list of arbitrary depth
	flatList = lambda *n: (e for a in n \
			for e in (flatList(*a) if isinstance(a, (tuple, list)) else (a,)))
	c_opts = list(flatList(c_opts))
	#print('c_opts = {}'.format(c_opts))
	c_opts = ' '.join(c_opts)
	print('c_opts = {}'.format(c_opts))
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

	print('Submitting {} files...'.format(len(batchList)))
	for i in range(len(batchList)):

		jobID = 'timescramble_{}_{:05d}'.format(args.config, i)
		cmd  = '{}/timeScramble'.format(cwd)
		ex = [cmd, '--batch_idx', str(i), c_opts]
		if not args.test:
			ex = [my.env] + ex
			#print('ex = {}'.format(ex))
		ex = [' '.join(ex)]
		#print('ex = {}'.format(ex))
		py_submit(ex, sublines=sublines, test=args.test, jobID=jobID)
	
	os.chdir(cwd)


