#!/usr/bin/env python

import glob, os, sys
from math import ceil

def maker(sim, batchSize=200, filetype='LLH'):

    # Get list of files
    gcd = '/data/sim/sim-new/downloads/GCD_31_08_11/GeoCalibDetectorStatus_IC79.55380_L2a.i3.gz'
    prefix = '/data/sim/IceTop/2010/filtered/level2a/CORSIKA-ice-top/'
    files = glob.glob(prefix + sim+'/*/Level2a_*.i3.bz2')
    files.sort()

    #njobs = (len(files)/batchSize)/2    # run over half the files
    njobs = int(ceil(len(files)/float(batchSize)))

    cwd = os.getcwd()
    os.chdir('/home/fmcnally/ShowerLLH/npx4/')

    for i in range(njobs):
        try:
            fileList = files[batchSize*i:batchSize*(i+1)]
        except IndexError:
            fileList = files[batchSize*i:]

        # Name outfile
        outPrefix = '/net/user/fmcnally/ShowerLLH/IT73_sim/files/'
        start = fileList[0][-13:-7]
        end   = fileList[-1][-13:-7]
        out = outPrefix + 'SimLLH_%s_%s_%s.hdf5' % (sim, start, end)
        print out

        fileList.insert(0, gcd)
        fileList = ' '.join(fileList)

        ex  = '/net/user/fmcnally/offline/V04-05-00/build/env-shell.sh'
        ex += ' python %s/Make%sFiles.py %s %s' % (cwd, filetype, out, fileList)
        os.system('./submit_npx4.sh ' + ex)
        #os.system(ex)

    os.chdir(cwd)



if __name__ == "__main__":

    if len(sys.argv)!=2:
        print 'Usage: %s [sim]' % sys.argv[0]

    sim = sys.argv[1]
    #batchSize = 200                         # 200 for 7006/7007 = 1/3 hr
    batchSize = 100                        # Smaller batch size for hi/lowE
    maker(sim, batchSize=batchSize, filetype='LLH')


