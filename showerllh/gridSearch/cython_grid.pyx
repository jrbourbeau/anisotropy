import numpy as np
import time
cimport numpy as np

""" Calculate array of LLH values for different energy bins
    Uses weave (C) to calculate closest approach distances and llh
    Fastest version of code
 """
def getLLH_C2(np.ndarray[double, ndim=2] grid, 
              np.ndarray[double, ndim=5] llhTable,
              np.ndarray[double, ndim=2] tankxyz,
              double theta, double phi,
              int Zbin,
              np.ndarray[int] Sbin, np.ndarray[int] Cbin,
              np.ndarray[double] Dbins):

    cdef double t0
    # General setup
    cdef int nG = grid.shape[0]
    cdef int ntanks = tankxyz.shape[0]
    cdef np.ndarray[double] RecoMaxLLHs = np.zeros(nG, dtype=np.double)
    cdef np.ndarray[int] argmax = np.zeros(nG, dtype=np.int32)
    cdef int i
    for i in range(nG):
        RecoMaxLLHs[i] = -np.inf
    cdef np.ndarray[double] tankx, tanky, tankz
    tankx, tanky, tankz = np.transpose(tankxyz)

    # Careful treatment of binning
    cdef nE = llhTable.shape[0]
    cdef np.ndarray[int] Dbin = np.zeros(ntanks, dtype=np.int32)

    cdef double Z = 1947
    cdef double ex = -np.sin(theta) * np.cos(phi)
    cdef double ey = -np.sin(theta) * np.sin(phi)
    cdef double ez = -np.cos(theta)

    cdef int j, bin, e, t
    cdef double X, Y, hx, hy, hz, s, x1, y1, z1, dist, llh

    for i in range(nG):

        # Calculate and bin closest approach distances
        X = grid[i][0]
        Y = grid[i][1]
        for j in range(ntanks):
            hx = X - tankx[j]
            hy = Y - tanky[j]
            hz = Z - tankz[j]
            s = ex*hx + ey*hy + ez*hz
            x1 = tankx[j] + s*ex
            y1 = tanky[j] + s*ey
            z1 = tankz[j] + s*ez
            dist = np.sqrt(pow(x1-X,2) + pow(y1-Y,2) + pow(z1-Z,2))
            #bin = np.digitize([dist], Dbins)[0] - 1
            bin = 0
            while dist > Dbins[bin+1]:
                bin += 1
            Dbin[j] = bin

        # Calculate llh for every energy bin, return most likely
        for e in range(nE):
            llh = 0.0
            for t in range(ntanks):
                llh += llhTable[e, Zbin, Sbin[t], Dbin[t], Cbin[t]]
            if llh > RecoMaxLLHs[i]:
                RecoMaxLLHs[i] = llh
                argmax[i] = e

    return RecoMaxLLHs, argmax

