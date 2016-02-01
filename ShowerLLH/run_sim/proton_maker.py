#!/usr/bin/env python

import glob, os, sys
from math import ceil

def maker(sim, batchSize=200, filetype='LLH'):

    # Get list of files
    gcd = '/data/sim/sim-new/downloads/GCD/GeoCalibDetectorStatus_IC86.55697_corrected_V2.i3.gz'
    prefix = '/data/sim/IceTop/2011/filtered/CORSIKA-ice-top/level2/'
    files = glob.glob(prefix + sim+'/*/Level2_*.i3.bz2')
    #files = ['/home/zgriffith/Level2_IC86_corsika_icetop.010719.000446.i3'] #gamma 10719 test file
    #files = ['/home/zgriffith/Level2_IC86_corsika_icetop.009166.000666.i3.bz2'] #proton 9166 test file
    #files = ['/home/zgriffith/Level2_IC86_corsika_icetop.009165.000327.i3.bz2'] #iron 9165 test file
    files.sort()

    #njobs = (len(files)/batchSize)/20    # run over a quarter of the files
    njobs = int(ceil(len(files)/float(batchSize)))

    cwd = os.getcwd()
    os.chdir('/home/zgriffith/ShowerLLH/npx4/')

    for i in range(njobs):
        try:
            fileList = files[batchSize*i:batchSize*(i+1)]
        except IndexError:
            fileList = files[batchSize*i:]

        # Name outfile
        outPrefix = '/net/user/zgriffith/ShowerLLH/IT81_sim/files/'
        start = fileList[0][-13:-7]
        end   = fileList[-1][-13:-7]
        out = outPrefix + 'SimLLH_%s_%s_%s.hdf5' % (sim, start, end)
        out_i3 = outPrefix + 'SimLLH_%s_%s_%s.i3' % (sim, start, end)
        print out

        fileList.insert(0, gcd)
        fileList = ' '.join(fileList)

        #ex  = '/home/zgriffith/icetray/V13_03/build/env-shell.sh'
        #ex   = '/net/user/zgriffith/offline/build/env-shell.sh'
        ex  = '/net/user/fmcnally/offline/V04-05-00/build/env-shell.sh'
        ex += ' python %s/Make%sFiles.py %s %s %s' % (cwd, filetype, out, out_i3, fileList)
        os.system('./submit_npx4.sh ' + ex)
        #os.system(ex)

    os.chdir(cwd)



if __name__ == "__main__":

    if len(sys.argv)!=2:
        print 'Usage: %s [sim]' % sys.argv[0]

    sim = sys.argv[1]
    batchSize = 200                         # 200 for 7006/7007 = 1/3 hr
    #batchSize = 100                        # Smaller batch size for hi/lowE
    maker(sim, batchSize=batchSize, filetype='LLH')


