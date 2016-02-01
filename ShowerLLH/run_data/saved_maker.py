#!/usr/bin/env python

import glob, os, sys
from math import ceil
sys.path.append('/home/zgriffith/ShowerLLH/useful/')
from getGoodRuns import fileCleaner

def maker(config, yyyymmdd, batchSize=20):

    yy = yyyymmdd[:4]
    md = yyyymmdd[4:]
    resourcedir = '/data/user/zgriffith/ShowerLLH/resources/'

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
        goodRunList = resourcedir + 'burn_sample.txt'

    # Sort and clean for good runs
    files += gcdFiles
    files.sort()
    files = fileCleaner(goodRunList, files)
    njobs = int(ceil(len(files)/float(batchSize)))
    #njobs = 1
    #job_list = [18,38]
    #job_list = [2,45]
    #job_list = [51,101,103,121,122]
    #job_list = [2,76]
    #job_list = [37,65,80]
    #job_list = [63,85]

    cwd = os.getcwd()
    os.chdir('/home/zgriffith/ShowerLLH/npx4/')

    for i in range(njobs):
        try:
            fList = files[batchSize*i:batchSize*(i+1)]
        except IndexError:
            fList = files[batchSize*i:]
        if 'GCD' not in fList[0]:
            st  = fList[0].find('Run') + 3
            run = fList[0][st:st+8]
            gcd = [gcdFile for gcdFile in gcdFiles if run in gcdFile][0]
            fList.insert(0, gcd)
        if 'GCD' in fList[-1]:
            fList.remove(fList[-1])

        # Name outfile
        outPrefix = '/data/user/zgriffith/ShowerLLH/'+config+'_data/files/burn_sample/'
        out = outPrefix + 'DataLLH_%s' % yyyymmdd
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
        #ex  = '/data/user/fmcnally/offline/V04-05-00/build/env-shell.sh'
        ex  = '/data/user/zgriffith/offline/build/env-shell.sh'
        ex += ' python %s/MakeLLHFiles.py %s %s %s' % (cwd, config, out, fList)
        os.system('./submit_npx4.sh ' + ex)
        #os.system(ex)

    os.chdir(cwd)



if __name__ == "__main__":

    if len(sys.argv)!=3:
        print 'Usage: %s [config] [yyyymmdd]' % sys.argv[0]
        sys.exit(1)

    config = sys.argv[1]
    yyyymmdd = sys.argv[2]
    batchSize = 80                         # ~ 2 hours per job

    maker(config, yyyymmdd, batchSize=batchSize)


