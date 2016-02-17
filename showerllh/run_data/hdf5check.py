#!/usr/bin/env python

import argparse, glob, tables
import myGlobals as my

if __name__ == "__main__":

    # Import global path names
    myGlobals.setupShowerLLH()

    p = argparse.ArgumentParser(
            description='Checks hdf5 files to make sure there are no errors')
    p.add_argument('-c', '--config', dest='config',
            default='IT73',
            help='Detector configuration')
    p.add_argument('-s', '--start', dest='start',
            help='Option for start of filenames')
    args = p.parse_args()

    prefix = '%s/%s_data/' % (my.llh_data, args.config)
    fileList = glob.glob(prefix + 'files/' + args.start + '*.hdf5')
    fileList.sort()

    for f in fileList:

        print 'Working on', f
        t = tables.openFile(f)
        try:
            test = t.root.ShowerLLH_proton.col('energy')
        except tables.exceptions.HDF5ExtError:
            print f, 'is bad'
        t.close()
