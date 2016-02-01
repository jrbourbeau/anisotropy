#!/usr/bin/env python

import glob, os, sys
from math import ceil

def maker(config, sim, batchSize=200, filetype='LLH'):

    if config == 'IT81':
        gcd = '/data/sim/sim-new/downloads/GCD/GeoCalibDetectorStatus_IC86.55697_corrected_V2.i3.gz'
        prefix = '/data/sim/IceTop/2011/generated/CORSIKA-ice-top/'#/level2/'
        files = glob.glob(prefix + sim+'/*/Detector_*.i3.gz')
         
    files.sort()
    njobs = int(ceil(len(files)/float(batchSize)))

    cwd   = os.getcwd()
    os.chdir('/home/zgriffith/npx4/')

    for i in range(njobs):
        try:
            fileList = files[batchSize*i:batchSize*(i+1)]
        except IndexError:
            fileList = files[batchSize*i:]

        # Name outfile
        outPrefix = '/data/user/zgriffith/ShowerLLH/'+config+'_sim/files/'

        if config == 'IT81':
            start = fileList[0][-12:-6]
            end   = fileList[-1][-12:-6]

        out = outPrefix + 'Raw_%s_%s_%s.hdf5' % (sim, start, end)
        print out

        fileList.insert(0, gcd)
        fileList = ' '.join(fileList)

        #ex  = '/data/user/zgriffith/offline/build/env-shell.sh'
        ex  = '/cvmfs/icecube.opensciencegrid.org/standard/RHEL_6_x86_64/metaprojects/simulation/trunk/env-shell.sh'
        ex += ' python %s/test_gammas.py %s %s %s' % (cwd, config, out, fileList)
        os.system('./submit_npx4.sh ' + ex)
        #os.system(ex)

    os.chdir(cwd)



if __name__ == "__main__":

    if len(sys.argv)!=3:
        print 'Usage: %s [sim]' % sys.argv[0]
        sys.exit(1)

    config    = sys.argv[1]
    sim       = sys.argv[2]
    batchSize = 50                 # 200 for 7006/7007 = 1/3 hr
    maker(config, sim, batchSize=batchSize, filetype='LLH')


