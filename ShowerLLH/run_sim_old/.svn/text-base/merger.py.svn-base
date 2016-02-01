#!/usr/bin/env python

import glob, os, sys

def merger(fileStart):

    prefix = '/net/user/fmcnally/ShowerLLH/IT73_sim/files/'
    fileList = glob.glob(prefix + 'SimLLH_'+fileStart+'_0*.hdf5')
    fileList.sort()

    if len(fileList) == 0:
        print 'No subfiles found for '+fileStart+'. Nothing to merge.'
        return

    fileList = ' '.join(fileList)
    outFile = prefix + 'SimLLH_'+fileStart+'.hdf5'

    hdfPrefix = '/net/user/fmcnally/offline/V04-05-00/build/hdfwriter/resources'
    ex = 'python %s/scripts/merge.py -o %s %s' % (hdfPrefix, outFile, fileList)
    os.system(ex)

if __name__ == "__main__":

    sim = sys.argv[1]

    if sim == 'All':
        simList  = ['7351','7006','7579']
        simList += ['7483','7241','7263','7791']
        simList += ['7486','7242','7262','7851']
        simList += ['7394','7007','7784']
        for sim in simList:
            merger(sim)

    else:
        merger(sim)

