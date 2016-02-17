#!/usr/bin/env python

import tables, glob
import numpy as np

if __name__ == "__main__":

    fileList = glob.glob('/data/user/fmcnally/showerllh/IT81-II_data/files/*')
    fileList = [f for f in fileList if '_00' not in f]
    fileList.sort()

    nlist = np.zeros(len(fileList))

    for i, f in enumerate(fileList):
        t = tables.openFile(f)
        n = t.root.ShowerLLH_proton.nrows
        nlist[i] = n
        print
        print f
        print n
        t.close()

    print
    print nlist.min(), nlist.mean(), nlist.max()
    
