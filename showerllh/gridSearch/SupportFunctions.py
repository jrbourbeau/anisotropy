#!/usr/bin/env python

import numpy as np
from scipy import weave
import tables, time

""" Converts a bin to average energy in GeV """
def b2e(Ebins, bin):
    step = (Ebins[1]-Ebins[0]) / 2.0
    center = Ebins[bin] + step
    return 10**center


""" Calculate the closest approach distance between a shower and tank """
def getDist(x0, y0, z0, theta, phi, x, y, z):

    # Get the unit vector of particle's direction
    ex = -np.sin(theta) * np.cos(phi)
    ey = -np.sin(theta) * np.sin(phi)
    ez = -np.cos(theta)
    # Vector between particle position and tank
    hx = x-x0
    hy = y-y0
    hz = z-z0
    # Distance between particle position and closest approach
    s = ex*hx + ey*hy + ez*hz
    # Closest approach position
    x1 = x0 + s*ex
    y1 = y0 + s*ey
    z1 = z0 + s*ez
    # Closest approach distance
    dist = np.sqrt((x1-x)**2 + (y1-y)**2 + (z1-z)**2)
    return dist.transpose()


""" Calculate array of LLH values for different energy bins.
    Pure python implementation (slow!!)
"""
def getLLH(grid, llhTable, bins, binnedVals, theta, phi, tankxyz):

    t0 = time.time()

    # General setup
    nG = int(len(grid))
    ntanks = int(len(tankxyz))
    RecoMaxLLHs = np.array([-np.inf for i in range(nG)], dtype=np.double)
    argmax = np.zeros(nG, dtype=np.int32)
    tankx, tanky, tankz = tankxyz.T

    # Careful treatment of binning
    # Every bin setup has to have energy, distance, and charge bins
    nE = int(len(bins['E']) - 1)
    idx = [0, binnedVals['Z'], binnedVals['S'], 0, binnedVals['C']]

    z = 1947.
    for i, (x,y) in enumerate(grid):

        # Calculate and bin closest approach distances
        closest_dists = getDist(x, y, z, theta, phi, tankx, tanky, tankz)
        idx[3] = np.digitize(closest_dists, bins['D']) - 1
        # Find and record
        for e in range(nE):
            idx[0] = e
            llh = 0.
            for p in range(ntanks):
                llh += llhTable[idx[0], idx[1], idx[2][p], idx[3][p], idx[4][p]]
            if llh > RecoMaxLLHs[i]:
                RecoMaxLLHs[i] = llh
                argmax[i] = e

    print 'Python:', time.time()-t0

    return RecoMaxLLHs, argmax


""" Calculate array of LLH values for different energy bins
    Uses weave (C) to calculate closest approach distances and llh
    Fastest version of code
 """
def getLLH_C3(grid, llhTable, bins, binnedVals, theta, phi, tankxyz):

    nDims = int(len(bins))
    binDims = np.array([len(bins[i][1])-1 for i in bins], dtype=np.int32)
    nG = int(len(grid))
    nE = int(binDims[0])
    ntanks = int(len(tankxyz))

    # Get indexes for energy and distance bins (special)
    binNames = {bins[k][0]:k for k in bins}
    e_idx = int(binNames['E'])
    d_idx = int(binNames['D'])
    Dbins = np.asarray(bins[d_idx][1], dtype=np.double)

    RecoMaxLLHs = np.array([-np.inf for i in range(nG)],dtype=np.double)
    argmax = np.zeros(nG, dtype=np.int32)

    argnames  = ['grid','llhTable','binnedVals','theta','phi','tankxyz']
    argnames += ['binDims','nDims','nG','nE','ntanks','e_idx','d_idx','Dbins']
    argnames += ['RecoMaxLLHs','argmax']

    loopcode = """

        // Starting parameters
        double Z = 1947;
        double ex = -sin(theta) * cos(phi);
        double ey = -sin(theta) * sin(phi);
        double ez = -cos(theta);

        // Loop over every grid position
        for (int i=0; i<nG; i++) {

            // Calculate and bin closest approach distances
            double X = grid(i,0);
            double Y = grid(i,1);

            for (int j=0; j<ntanks; j++) {
                double hx = X - tankxyz(j,0);
                double hy = Y - tankxyz(j,1);
                double hz = Z - tankxyz(j,2);
                double s = ex*hx + ey*hy + ez*hz;
                double x1 = tankxyz(j,0) + s*ex;
                double y1 = tankxyz(j,1) + s*ey;
                double z1 = tankxyz(j,2) + s*ez;
                double dist = sqrt(pow(x1-X,2) + pow(y1-Y,2) + pow(z1-Z,2));
                int bin = 0;
                while (dist > Dbins(bin+1))
                    bin += 1;
                binnedVals(j*nDims + d_idx) = bin;
            }

            // Calculate llh for every energy bin, return most likely
            int idx, c;
            for (int e=0; e<nE; e++) {
                double llh = 0.0;
                for (int p=0; p<ntanks; p++) {
                    binnedVals(p*nDims + e_idx) = e;
                    idx = 0;
                    c = 1;
                    for (int k=nDims; k-- > 0; ) {
                        idx += c * binnedVals(p*nDims + k);
                        c *= binDims(k);
                    }
                    llh += llhTable(idx);
                }
                if (llh > RecoMaxLLHs(i)) {
                    RecoMaxLLHs(i) = llh;
                    argmax(i) = e;
                }
            }
        }
        return_val = 0;
        """

    weave.inline(loopcode, argnames, type_converters=weave.converters.blitz)

    return RecoMaxLLHs, argmax


""" Calculate array of LLH values for different energy bins
    Uses weave (C) to calculate closest approach distances and llh
    Fastest version of code
 """
def getLLH_C(grid, llhTable, Zbin, Sbin, Cbin, theta, phi, tankxyz, Dbins):

    # 5-D hists:  ( E x Z x S x D x C )
    nG = int(len(grid))
    nE = int(len(llhTable))
    ntanks = int(len(tankxyz))
    RecoMaxLLHs = np.array([-np.inf for i in range(nG)], dtype=np.double)
    argmax = np.zeros(nG, dtype=np.int32)
    Dbin = np.zeros(ntanks, dtype=np.int32)
    tankx, tanky, tankz = tankxyz.T

    loopcode = """

        // Starting parameters
        double Z = 1947;
        double ex = -sin(theta) * cos(phi);
        double ey = -sin(theta) * sin(phi);
        double ez = -cos(theta);

        // Loop over every grid position
        for (int i=0; i<nG; i++) {

            // Calculate and bin closest approach distances
            double X = grid(i,0);
            double Y = grid(i,1);

            for (int j=0; j<ntanks; j++) {
                double hx = X - tankx(j);
                double hy = Y - tanky(j);
                double hz = Z - tankz(j);
                double s = ex*hx + ey*hy + ez*hz;
                double x1 = tankx(j) + s*ex;
                double y1 = tanky(j) + s*ey;
                double z1 = tankz(j) + s*ez;
                double dist = sqrt(pow(x1-X,2) + pow(y1-Y,2) + pow(z1-Z,2));
                int bin = 0;
                while (dist > Dbins(bin+1)) {
                    bin += 1;
                }
                Dbin(j) = bin;
            }

            // Calculate llh for every energy bin, return most likely
            for (int e=0; e<nE; e++) {
                double llh = 0.0;
                for (int p=0; p<ntanks; p++) {
                    llh += llhTable(e, Zbin, Sbin(p), Dbin(p), Cbin(p));
                }
                if (llh > RecoMaxLLHs(i)) {
                    RecoMaxLLHs(i) = llh;
                    argmax(i) = e;
                }
            }
        }
        return_val = 0;
        """

    weave.inline(loopcode, ['nG','grid','theta','phi','ntanks','tankx','tanky','tankz','Dbins','nE','llhTable','RecoMaxLLHs','Zbin','Sbin','Dbin','Cbin','argmax'], type_converters=weave.converters.blitz)

    return RecoMaxLLHs, argmax


""" Calculate array of LLH values for different energy bins 
    Exclude zenith binning for anisotropy study
"""
def getLLH_anisotropy(grid, llhTables, Sbin, Cbin, theta, phi, tankx, tanky, tankz, Dbins):

    t0 = time.time()

    # 4-D hists:  ( E x S x D x C )
    nG = int(len(grid))
    nE = int(len(llhTables))
    ntanks = int(len(tankx))
    RecoMaxLLHs = np.array([-np.inf for i in range(nG)], dtype=np.double)
    argmax = np.zeros(nG, dtype=np.int32)
    Dbin = np.zeros(ntanks, dtype=np.int32)

    loopcode = """

        // Starting parameters
        double Z = 1947;
        double ex = -sin(theta) * cos(phi);
        double ey = -sin(theta) * sin(phi);
        double ez = -cos(theta);

        // Loop over every grid position
        for (int i=0; i<nG; i++) {

            // Calculate and bin closest approach distances
            double X = grid(i,0);
            double Y = grid(i,1);

            for (int j=0; j<ntanks; j++) {
                double hx = X - tankx(j);
                double hy = Y - tanky(j);
                double hz = Z - tankz(j);
                double s = ex*hx + ey*hy + ez*hz;
                double x1 = tankx(j) + s*ex;
                double y1 = tanky(j) + s*ey;
                double z1 = tankz(j) + s*ez;
                double dist = sqrt(pow(x1-X,2) + pow(y1-Y,2) + pow(z1-Z,2));
                int bin = 0;
                while (dist > Dbins(bin+1)) {
                    bin += 1;
                }
                Dbin(j) = bin;
            }

            // Calculate llh for every energy bin, return most likely
            for (int e=0; e<nE; e++) {
                double llh = 0.0;
                for (int p=0; p<ntanks; p++) {
                    llh += llhTables(e, Sbin(p), Dbin(p), Cbin(p));
                }
                if (llh > RecoMaxLLHs(i)) {
                    RecoMaxLLHs(i) = llh;
                    argmax(i) = e;
                }
            }
        }
        return_val = 0;
        """

    weave.inline(loopcode, ['nG','grid','theta','phi','ntanks','tankx','tanky','tankz','Dbins','nE','llhTables','RecoMaxLLHs','Sbin','Dbin','Cbin','argmax'], type_converters=weave.converters.blitz)

    #print 'C: ', time.time()-t0

    return RecoMaxLLHs, argmax


""" Theoretically for use with MCMC """
def getLLH2(hists, varList, binList):
    phist, fhist = hists[0], hists[1]
    # 5-D hists: ( C x S x Z x E x D )
    Ebin  = int(digitize([varList[0]], binList[0])[0] - 1)
    Zbin  = int(digitize([varList[1]], binList[1])[0] - 1)
    Sbin  = digitize(varList[2], binList[2]).astype(int32) - 1
    Dbin  = digitize(varList[3], binList[3]).astype(int32) - 1
    Cbin  = digitize(varList[4], binList[4]).astype(int32) - 1
    tanks = len(Sbin)

    pllh, fllh = 0.0, 0.0
    loopcode = """
        for (int p=0; p<tanks; p++) {
            pllh += phist(Cbin(p), Sbin(p), Zbin, Ebin, Dbin(p));
            fllh += fhist(Cbin(p), Sbin(p), Zbin, Ebin, Dbin(p));
        }
        return_val = 0;
        """
    weave.inline(loopcode, ['phist','fhist','pLLHs','fLLHs','Zbin','Sbin','Dbin','Cbin','tanks'], type_converters=weave.converters.blitz)

    return pLLHs, fLLHs

