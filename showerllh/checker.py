#!/usr/bin/env python

import numpy as np
import glob, re

import myGlobals as my

if __name__ == "__main__":

    my.setupShowerLLH(verbose=False)

    prefix = '%s/IT73_sim/files' % my.llh_data
    fileList = glob.glob('%s/*_part*' % prefix)
    fileList.sort()

    bases = ['_'.join(f.split('_')[:-1]) for f in fileList]
    bases = sorted(list(set(bases)))

    for base in bases:

        print 'Working on %s...' % base
        files = glob.glob('%s_part*' % base)
        files.sort()

        completed = []
        for file in files:
            start, end = re.split('_|\.|-|part',file)[-3:-1]
            completed += range(int(start), int(end)+1)

        print 'Found %i files, ranging from %i to %i' % \
                (len(files), completed[0], completed[-1])
        fullrange = range(completed[0], completed[-1]+1)
        missing = [i for i in fullrange if i not in completed]
        missing.sort()

        if len(missing) == 0:
            print 'No missing parts found'
            continue

        # Split according to continuous ranges
        starts, ends = [],[]
        starts += [missing[0]]
        for i in range(1,len(missing)):
            if missing[i] != missing[i-1]+1:
                starts += [missing[i]]
                ends += [missing[i-1]]
        ends += [missing[-1]]
        for start, end in zip(starts, ends):
            print 'Missing %i to %i' % (start, end)
    

