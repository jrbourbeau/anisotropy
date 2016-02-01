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

    # Sample bad file
    bad = '/net/user/zgriffith/ShowerLLH/IT73_data/files/DataLLH_201006_00115995_00000045_00115996_00000019.hdf5'
    badFile = os.path.basename(bad)
    badRun0 = badFile[20:28]
    badPart0 = badFile[29:37]
    badRun1
    badPart1

    newList
    for file in files:
        if Run0 >= badRun0 and Part0 >= badPart0 and Run1 <= badRun1 and Part1 <= badPart1:
            newList.append(file)

    newgcdList
    for file in gcdfiles:
        if Run0 >= badRun0 and Run1 <= badRun1:
            newListgcd.append(file)

    # Sort and clean for good runs
    files += gcdFiles
    files.sort()
    files = fileCleaner(goodRunList, files)
    print len(files)
    #njobs = (len(files)/batchSize)/2    # run over half the files
    njobs = int(ceil(len(files)/float(batchSize)))
    #njobs = 1

    cwd = os.getcwd()
    os.chdir('/home/zgriffith/ShowerLLH/npx4/')

    ## Warning: this was range(6,njobs)! Need to check if I'm missing
    ## information from beginning of each month
    ## Still need to check IT73 data
    for badFile in badFileList:


        if 'GCD' not in fList[0]:
            st  = fList[0].find('Run') + 3
            run = fList[0][st:st+8]
            gcd = [i for i in gcdFiles if run in i][0]
            fList.insert(0, gcd)
        if 'GCD' in fList[-1]:
            fList.remove(fList[-1])

        # Name outfile
        outPrefix = '/net/user/zgriffith/ShowerLLH/'+config+'_data/files/'
        out = outPrefix + 'DataLLH_%s' % yyyymmdd
        for file in [fList[1], fList[-1]]:
            run_st  = file.find('Run') + 3
            subrun = 'Part' if 'Part' in file else 'Subrun'
            subrun_st = file.find(subrun) + len(subrun)
            run  = file[run_st:run_st+8]
            part = file[subrun_st:subrun_st+8]
            out += '_%s_%s' % (run, part)
        out += '.hdf5'
        print out

        fList = ' '.join(fList)
        ex  = '/net/user/fmcnally/offline/V04-05-00/build/env-shell.sh'
        ex += ' python %s/MakeLLHFiles.py %s %s %s' % (cwd, config, out, fList)
        #os.system('./submit_npx4.sh ' + ex)
        os.system(ex)

    os.chdir(cwd)



if __name__ == "__main__":

    if len(sys.argv)!=3:
        print 'Usage: %s [config] [yyyymmdd]' % sys.argv[0]
        sys.exit(1)

    config = sys.argv[1]
    yyyymmdd = sys.argv[2]
    batchSize = 80                          # ~ 2 hours per job

    maker(config, yyyymmdd, batchSize=batchSize)


