#!/usr/bin/env python

import glob, os, re
import numpy as np

if __name__ == "__main__":

    hdfPrefix='/data/user/fmcnally/offline/V04-08-00/build/hdfwriter/resources'
    prefix = '/data/user/fmcnally/anisotropy/sim'
    fileList = glob.glob('%s/*.hdf5' % prefix)
    fileList.sort()

    simList = []
    for file in fileList:
        sim = os.path.basename(file)
        sim = '_'.join(re.split('_|\.', sim)[:2])
        if sim not in simList:
            simList += [sim]

    # Files must be of the form [basename_start_end.suffix]
    def recursiveMerge(depth, files):

        if len(files) == 1:
            print 'Nothing to be done for', files[0]
            return

        fileSize = sum(np.array([os.stat(f).st_size for f in files]))
        if fileSize > 1e9:
            depth += 1
            tempLists = [files[:len(files)/2], files[len(files)/2:]]
            for tempList in tempLists:
                recursiveMerge(depth, tempList)

        base  = '_'.join(files[0].split('_')[:-2])
        start = files[0].split('_')[-2]
        end   = files[-1].split('_')[-1]
        outFile = '_'.join([base, start, end])
        if depth == 0:
            outFile = base + '.hdf5'
        if os.path.isfile(outFile):
            print 'Outfile %s already exists.' % outFile
            print 'Skipping...'
            return

        files = ' '.join(files)
        ex  = 'python %s/scripts/merge.py' % hdfPrefix
        ex += ' -o %s %s' % (outFile, files)
        os.system(ex)


    for sim in simList:

        fileList = glob.glob('%s/%s_*.hdf5' % (prefix, sim))
        fileList.sort()

        if len(fileList) == 0:
            print 'No subfiles found for '+sim+'. Nothing to merge.'
            continue

        recursiveMerge(0, fileList)


