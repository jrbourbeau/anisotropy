#!/usr/bin/env python

import glob, os
from math import ceil

def maker(simList, batchSize=1000):

    # Basic setup for IT8122
    prefix = '/data/sim/IceTop/2011/filtered/CORSIKA-ice-top/level2/'
    outPrefix = '/net/user/zgriffith/ShowerLLH/resources/'
    gcd = '/data/sim/sim-new/downloads/GCD/' + \
                'GeoCalibDetectorStatus_IC86.55697_corrected_V2.i3.gz'

    cwd = os.getcwd()
    os.chdir('/home/zgriffith/ShowerLLH/npx4/')

    for sim in simList:

        # Get list of files
        files = glob.glob(prefix + sim+'/*/Level2_*.i3.gz')
        files.sort()

        #njobs = (len(files/batchSize)/2        # run over half the files
        njobs = int(ceil(len(files)/float(batchSize)))

        for i in range(njobs):
            try:
                fileList = files[batchSize*i:batchSize*(i+1)]
            except IndexError:
                fileList = files[batchSize*i:]

            # Name outfile
            start = fileList[0][-12:-6]
            end   = fileList[-1][-12:-6]
            outFile = outPrefix+'IT81/CountTable_'+sim+'_'+start+'_'+end+'.npy'
            print outFile

            fileList.insert(0, gcd)
            fileList = ' '.join(fileList)

            ex  = '/net/user//offline/V04-05-00/build/env-shell.sh'
            ex += ' python %s/MakeHist.py %s %s' % (cwd, outFile, fileList)
            os.system('./submit_npx4.sh ' + ex)
            #os.system(ex)

    os.chdir(cwd)


if __name__ == "__main__":

    comp = raw_input('Gammas, Proton, Helium, Oxygen, or Iron? [G|P|He|O|Fe]: ')
    batchSize = 1000

    # Choose simulation to run over
    simDict = {}
    simDict['G']  = ['10687']
    simDict['P']  = ['7351','7006','7579']
    simDict['He'] = ['7483','7241','7263','7791']
    simDict['O']  = ['7486','7242','7262','7851']
    simDict['Fe'] = ['7394','7007','7784']
    simList = simDict[comp]
    #simList = ['10687']

    maker(simList, batchSize=batchSize)








