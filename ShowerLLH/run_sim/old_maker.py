#!/usr/bin/env python

import glob, os, sys
from math import ceil

def maker(config, sim, batchSize=200, filetype='LLH'):

    # Get list of files
    if config == 'IT73':
        gcd = '/data/sim/sim-new/downloads/GCD_31_08_11/GeoCalibDetectorStatus_IC79.55380_L2a.i3.gz'
        prefix = '/data/sim/IceTop/2010/filtered/level2a/CORSIKA-ice-top/'
        files = glob.glob(prefix + sim+'/*/Level2a_*.i3.bz2')
    if config == 'IT81':
        gcd = '/data/sim/sim-new/downloads/GCD/GeoCalibDetectorStatus_IC86.55697_corrected_V2.i3.gz'
        prefix = '/data/sim/IceTop/2011/filtered/CORSIKA-ice-top/level2/'
        files = glob.glob(prefix + sim+'/*/Level2_*.i3.gz')
    
    files.sort()
    njobs = int(ceil(len(files)/float(batchSize)))

    cwd   = os.getcwd()
    os.chdir('/home/zgriffith/ShowerLLH/npx4/')

    for i in range(njobs):
        try:
            fileList = files[batchSize*i:batchSize*(i+1)]
        except IndexError:
            fileList = files[batchSize*i:]

        # Name outfile
        outPrefix = '/data/user/zgriffith/ShowerLLH/'+config+'_sim/files/'

        if config == 'IT73':
            start = fileList[0][-13:-7]
            end   = fileList[-1][-13:-7]
        if config == 'IT81':
            start = fileList[0][-12:-6]
            end   = fileList[-1][-12:-6]

        out = outPrefix + 'SimLLH_%s_%s_%s.hdf5' % (sim, start, end)
        print out

        fileList.insert(0, gcd)
        fileList = ' '.join(fileList)

        ex  = '/data/user/zgriffith/offline/build/env-shell.sh'
        ex += ' python %s/Make%sFiles.py %s %s %s' % (cwd, filetype, config, out, fileList)
        os.system('./submit_npx4.sh ' + ex)
        #os.system(ex)

    os.chdir(cwd)



if __name__ == "__main__":

    if len(sys.argv)!=3:
        print 'Usage: %s [sim]' % sys.argv[0]
        sys.exit(1)

    config    = sys.argv[1]
    sim       = sys.argv[2]
    batchSize = 50                  # 200 for 7006/7007 = 1/3 hr
    maker(config, sim, batchSize=batchSize, filetype='LLH')


