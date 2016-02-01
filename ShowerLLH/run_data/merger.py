#!/usr/bin/env python

import glob, os, sys

if __name__ == "__main__":

    if len(sys.argv)!=3:
        print 'Usage: python %s [config] [yyyymmdd]' % sys.argv[0]
        sys.exit(2)

    config   = sys.argv[1]
    yyyymmdd = sys.argv[2]

    prefix = '/data/user/zgriffith/ShowerLLH/'+config+'_data/files/burn_sample/'
    outFile = prefix + 'DataLLH_'+yyyymmdd+'.hdf5'
    #outFile = prefix + 'DataLLH_burn_sample.hdf5'

    #fileList = glob.glob(prefix + 'DataLLH_'+yyyymmdd+'???.hdf5')
    fileList = glob.glob(prefix + 'DataLLH_'+yyyymmdd+'*.hdf5')
    fileList.sort()
    fileList = ' '.join(fileList)

    hdfPrefix = '/data/user/fmcnally/offline/V04-08-00/build/hdfwriter/resources'
    ex = 'python %s/scripts/merge.py -o %s %s' % (hdfPrefix, outFile, fileList)
    os.system(ex)


## Original bash code for running merge and mix ##
#python "${I3}/hdfwriter/resources/scripts/merge.py" "-o" "$outfile" $fileList
#python "${I3}/hdfwriter/resources/scripts/mix.py" "${prefix}/DataLLH_${md}.hdf5" "${prefix}/DataLLH_${md}.hdf5"






