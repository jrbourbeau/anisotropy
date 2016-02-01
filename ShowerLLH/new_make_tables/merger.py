#!/usr/bin/env python

import glob, argparse, re
import myGlobals as my
import numpy as np


if __name__ == "__main__":

    # Global variables setup for path names
    my.setupShowerLLH(verbose=False)
    resourcedir = my.llh_resource

    p = argparse.ArgumentParser(
            description='Merges all count tables in a given directory')
    p.add_argument('-p', '--prefix', dest='prefix',
            default=resourcedir+'/CountTables/',
            help='Location of CountTables to merge')
    args = p.parse_args()

    # Build list of simulations
    masterList = glob.glob(args.prefix + 'CountTable_*_Part??????-??????.npy')
    simList = ['_'.join(f.split('_')[:-1]) for f in masterList]
    simList = sorted(list(set(simList)))

    for sim in simList:

        outFile = '%s.npy' % sim
        fileList = glob.glob(sim + '_Part??????-??????.npy')
        fileList.sort()

        if len(fileList) == 0:
            raise SystemError('No subfiles found for '+sim+'.')

        d = {}
        for i, file in enumerate(fileList):
            print 'Loading', file
            q = np.load(file)
            q = q.item()
            if i == 0:
                d['bins']   = q['bins']
                d['counts'] = np.zeros(q['counts'].shape)
            #if cmp(q['bins'], d['bins']) != 0:
            #    print '%s bins do not match %s bins' % (file, fileList[0])
            #    print 'Skipping...'
            #    continue
            d['counts'] += q['counts']

        np.save(outFile, d)

