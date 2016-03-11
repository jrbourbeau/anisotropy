#!/usr/bin/env python

import os, glob

def maker():

    cwd = os.getcwd()
    prefix = '/data/user/fmcnally/anisotropy/maps/merged'
    outPrefix = '/data/user/fmcnally/anisotropy/maps/dq-fit'

    #fileList = glob.glob(prefix + '/IC?*_24H_sid.fits')
    fileList = glob.glob(prefix + '/IC_24H_sid.fits')
    fileList.sort()

    output = "true"

    for inFile in fileList[:1]:

        print inFile
        outFile = outPrefix + '/' + os.path.basename(inFile)
        outFile = outFile.replace('.fits', '_dqfit.fits')
        if os.path.isfile(outFile):
            os.remove(outFile)

        cmd = '%s/dq-fit' % cwd
        args = ' '.join([inFile, outFile, output])
        ex   = ' '.join([cmd, args])

        os.system(ex)


if __name__ == "__main__":

    maker()
