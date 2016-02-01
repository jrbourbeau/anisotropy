#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import pickle

def plotLLH(c_llhs, m_llhs, f_llhs, cmf):

    fl = open('/net/user/fmcnally/ShowerLLH/resources/BinSpots.pkl', 'rb')
    grids = pickle.load(fl)
    fl.close()

    if cmf == 'coarse':
        x, y = np.transpose(grids['IT73']['coarse'])
        z = c_llhs
    if cmf == 'middle':
        c_index = c_llhs.argmax()
        x, y = np.transpose(grids['IT73']['middle'][c_index])
        z = m_llhs
    if cmf == 'fine':
        c_index = c_llhs.argmax()
        m_index = m_llhs.argmax()
        x, y = np.transpose(grids['IT73']['fine'][c_index][m_index])
        z = f_llhs

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(x, y, c=z, marker='H', s=200)
    plt.show()

if __name__ == "__main__":

    import tables
    t = tables.openFile('/net/user/fmcnally/ShowerLLH/IT73_sim/SimLLH_7007_00001_00010.hdf5')

    c_llhs = t.root.ShowerLLHParams.col('pCoarse_llhs')[0]
    m_llhs = t.root.ShowerLLHParams.col('pMed_llhs')[0]
    f_llhs = t.root.ShowerLLHParams.col('pFine_llhs')[0]

    plotLLH(c_llhs, m_llhs, f_llhs, 'fine')



