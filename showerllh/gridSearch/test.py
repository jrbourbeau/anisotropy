#!/usr/bin/env python

import scipy.weave as weave
import numpy as np

def testFunc(grid, hist, binNames, binVals, theta, phi, tankxyz, Dbins):

    gridx, gridy = np.transpose(grid)
    tankx, tanky, tankz = np.transpose(tankxyz)
    nDims = len(binNames)
    temp_list = []

    RecoMaxLLHs = np.array([-np.inf for i in range(len(gridx))],dtype=np.double)
    argmax = np.zeros(len(gridx), dtype=np.int32)

    argnames  = ['gridx','gridy','hist','binNames','binVals','theta','phi']
    argnames += ['nDims','RecoMaxLLHs','argmax']
    argnames += ['tankx','tanky','tankz','Dbins','temp_list']

    loopcode = """

        // Starting parameters
        double Z = 1947;
        double ex = -sin(theta) * cos(phi);
        double ey = -sin(theta) * sin(phi);
        double ez = -cos(theta);

        // Loop over every grid position
        for (int i=0; i<Ngridx[0]; i++) {

            // Calculate and bin closest approach distances
            double X = gridx[i];
            double Y = gridy[i];

            for (int j=0; j<Ntankx[0]; j++) {
                double hx = X - tankx[j];
                double hy = Y - tanky[j];
                double hz = Z - tankz[j];
                double s = ex*hx + ey*hy + ez*hz;
                double x1 = tankx[j] + s*ex;
                double y1 = tanky[j] + s*ey;
                double z1 = tankz[j] + s*ez;
                double dist = sqrt(pow(x1-X,2) + pow(y1-Y,2) + pow(z1-Z,2));
                int bin = 0;
                while (dist > Dbins[bin+1])
                    bin += 1;
                binVals[j*nDims + 3] = bin;
            }

            // Calculate llh for every energy bin, return most likely
            int idx, c, binVal, dim;
            for (int e=0; e<Nhist[0]; e++) {
                double llh = 0.0;
                for (int p=0; p<Ntankx[0]; p++) {
                    binVals[p*nDims] = e;
                    idx = 0;
                    c = 1;
                    for (int k=nDims; k-- > 0; ) {
                        binVal = binVals[p*nDims + k];
                        idx += c * binVal;
                        dim = Nhist[k];
                        c *= dim;
                    }
                    llh += hist[idx];
                }
                if (llh > RecoMaxLLHs[i]) {
                    RecoMaxLLHs[i] = llh;
                    argmax[i] = e;
                }
            }
        }
        return_val = 0;
        """

    weave.inline(loopcode, argnames)
    return RecoMaxLLHs, argmax


if __name__ == "__main__":

    # Real parameters
    gridFile = '/data/user/fmcnally/ShowerLLH/resources/IT73_grid.npy'
    grids = np.load(gridFile)
    grids = grids.item()
    grid = grids[0]
    llhFile = '/data/user/fmcnally/ShowerLLH/resources/LLHTables.npy'
    llhInfo = np.load(llhFile)
    llhInfo = llhInfo.item()
    bins = llhInfo['bins']
    llhTable = llhInfo['llhtables']
    Dbins = bins['D']

    # Fake parameters
    ntanks = 5
    Zbin = 0
    Sbin = np.random.randint(0, 4, ntanks)
    Cbin = np.random.randint(0, 20, ntanks)
    Sbin = [i for i in Sbin]
    Cbin = [i for i in Cbin]
    tankxyz = np.random.uniform(-500, 500, (ntanks,3))
    theta = 0.1
    phi = 0.2
    binNames = ['E','Z','S','D','C']
    binVals = np.array([[0,Zbin,Sbin[p],0,Cbin[p]] for p in range(len(Sbin))])
    binVals = binVals.flatten()
    binVals = [i for i in binVals]

    t = testFunc(grid, llhTable['proton'], binNames, binVals, theta, phi, tankxyz, Dbins)
    print t


