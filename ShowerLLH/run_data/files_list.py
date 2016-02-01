#!/usr/bin/env python

import glob, os, sys
from math import ceil
sys.path.append('/home/fmcnally/ShowerLLH/useful/')
from getGoodRuns import fileCleaner

def maker(config, yyyymmdd, batchSize=20):

    yy = yyyymmdd[:4]
    md = yyyymmdd[4:]
    resourcedir = '/net/user/fmcnally/ShowerLLH/resources/'

    # Get list of files
    if config == 'IT73':
        prefix = '/data/exp/IceCube/%s/filtered/level2a/' % (yy)
        files = glob.glob(prefix + md + '*/Level2a_*_IT.i3.bz2')
        gcdFiles = glob.glob(prefix + md + '*/Level2a_*_GCD.i3.bz2')
        goodRunList = resourcedir + 'IC79_GoodRunList_v2pt9.txt'
    if config == 'IT81':
        prefix = '/data/exp/IceCube/%s/filtered/level2/' % (yy)
        files = glob.glob(prefix + md + '*/Level2_*Part*_IT.i3.bz2')
        gcdFiles = glob.glob(prefix + md + '*/Level2_*_GCD.i3.gz')
        goodRunList = resourcedir + 'IC86_GoodRunList_v1pt5a.txt'
    if config == 'IT81-II':
        prefix = '/data/exp/IceCube/%s/filtered/level2/' % (yy)
        files = glob.glob(prefix + md + '*/Level2_*Subrun*_IT.i3.bz2')
        gcdFiles = glob.glob(prefix + md + '*/Level2_*_GCD.i3.gz')
        goodRunList = resourcedir + 'IC86_GoodRunList_v1pt5b.txt'

    # Sort and clean for good runs
    files += gcdFiles
    files.sort()
    files = fileCleaner(goodRunList, files)
    #njobs = (len(files)/batchSize)/2    # run over half the files
    njobs = int(ceil(len(files)/float(batchSize)))
    #njobs = 1
    #job_list = [12,13,19,35,36,37,38,39,41,42,44,45,62,83,106,109]
    #job_list = [18,50,114]
    #job_list = [2,45]
    #job_list = [51,101,103,121,122]
    #job_list = [2,76]
    #job_list = [37,65,80]
    #job_list = [63,85]

    cwd = os.getcwd()

    ## Warning: this was range(6,njobs)! Need to check if I'm missing
    ## information from beginning of each month
    ## Still need to check IT73 data
    for i in range(0,njobs):
        try:
            fList = files[batchSize*i:batchSize*(i+1)]
        except IndexError:
            fList = files[batchSize*i:]

        # Name outfile
        outPrefix = '/net/user/zgriffith/ShowerLLH/'+config+'_data/files/'
        out = outPrefix + 'DataLLH_%s' % yyyymmdd
        print i, fList[0][-34:-10], fList[-1][-34:-10], '\n'

        #os.system('./submit_npx4.sh ' + ex)

    os.chdir(cwd)



if __name__ == "__main__":

    if len(sys.argv)!=3:
        print 'Usage: %s [config] [yyyymmdd]' % sys.argv[0]
        sys.exit(1)

    config = sys.argv[1]
    yyyymmdd = sys.argv[2]
    batchSize = 80                          # ~ 2 hours per job

    maker(config, yyyymmdd, batchSize=batchSize)


