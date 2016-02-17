#!/usr/bin/env python

import glob
import numpy as np

if __name__ == "__main__":

    prefix = '/data/user/fmcnally/ShowerLLH/resources/CountTables/'
    fileList = glob.glob(prefix + '*.npy')
    for file in fileList:
        d = np.load(file)
        d = d.item()
        bins = d['bins']
        print bins

        binOrder = {0:'E',1:'Z',2:'S',3:'D',4:'C'}
        new_bins = {i:[binOrder[i], bins[binOrder[i]]] for i in binOrder}

        d['bins'] = new_bins

        #np.save(file, d)
