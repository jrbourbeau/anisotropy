#!/usr/bin/env python

import numpy as np
import glob, os


if __name__ == "__main__":

    prefix = '/data/user/fmcnally/showerllh'
    configs = ['IT73','IT81','IT81-II','IT81-III','IT81-IV']

    for cfg in configs:

        tempPrefix = '%s/%s_data/hists' % (prefix, cfg)
        files = glob.glob('%s/*_sky.npy' % tempPrefix)
        files.sort()

        bases = [os.path.basename(f) for f in files]
        bases = [b.split('_')[1][:6] for b in bases]
        monthList = sorted(list(set(bases)))

        for month in monthList:
            monthFiles = [f for i, f in enumerate(files) if bases[i]==month]
            for j, m in enumerate(monthFiles):

                print 'Working on %s...' % m
                temp = np.load(m)
                temp = temp.item()
                if j == 0:
                    d = temp
                else:
                    for key in temp:
                        d[key] += temp[key]

            print 'Writing to file...'
            outFile = '%s/DataLLH_%s_logdist_sky.npy' % (tempPrefix, month)
            print outFile
            np.save(outFile, d)

