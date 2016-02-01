#!/usr/bin/env python

from numpy import *
import glob 

def merger(fileStart):

    prefix = '/net/user/zgriffith/ShowerLLH/resources/IT81/'
    fileList = glob.glob(prefix + 'CountTable_'+fileStart+'_0*.npy')
    fileList.sort()

    if len(fileList) == 0:
        print 'No subfiles found for '+fileStart+'. Nothing to merge.'
        return

    outFile = prefix + 'CountTable_'+fileStart+'.npy'

    for i in range(len(fileList)):
        print 'Loading', fileList[i]
        q = load(fileList[i])
        if i == 0:
            d = zeros(q.shape)
        d += q

    save(outFile, d)


if __name__ == "__main__":

    simList  = ['10687']                            # Gamma
    '''
    simList  = ['7351','7006','7579']               # Proton
    simList += ['7483','7241','7263','7791']        # Helium
    simList += ['7486','7242','7262','7851']        # Oxygen
    simList += ['7394','7007','7784']               # Iron
    '''
    for sim in simList:
        merger(sim)


