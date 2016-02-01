#!/usr/bin/env python 

################################################################################
# Runs the grid search on the given files and writes to hdf5 output
################################################################################

from icecube import dataio, toprec, ShowerLLH
from I3Tray import *
from icecube.tableio import I3TableWriter
from icecube.hdfwriter import I3HDFTableService

import numpy as np
import sys, time, glob

sys.path.append('/home/fmcnally/ShowerLLH/segments')
from modules import FindLoudestStation, LoudestStationOnEdge, LargestTankCharges
from modules import GetStations, PruneIceTop, moveMCPrimary

def maker(outFile, fileList):

    # Starting parameters
    recoPulses = 'CleanedHLCTankPulses'
    config = 'IT73'
    resourcedir = '/net/user/fmcnally/ShowerLLH/resources/'
    compDict = {'P':'proton', 'He':'helium', 'O':'oxygen', 'Fe':'iron'}

    # Load grids for grid scan, llh tables, and bins
    grids = np.load(resourcedir + config+'_BinSpots.npy')
    grids = grids.item()
    LLHTables = np.load(resourcedir + 'LLHTables_ztest.npy')
    LLHTables = LLHTables.item()
    binDict = np.load(resourcedir + 'ShowerLLH_bins_ztest.npy')
    binDict = binDict.item()

    # Keys to write to frame
    keys  = []
    keys += ['ShowerPlane', 'ShowerPlaneParams']
    keys += ['LoudestStation', 'LoudestOnEdge']
    keys += ['SaturationList', 'SaturatedOnEdge']
    keys += ['Q1', 'Q2', 'Q3', 'Q4']
    keys += [recoPulses, 'NStations']
    for comp in LLHTables.keys():
        keys += ['ShowerLLH_'+compDict[comp]]
    keys += ['ShowerLLHParams']
    keys += ['MCPrimary']


    t0 = time.time()

    tray = I3Tray()
    tray.AddModule('I3Reader', FileNameList=fileList)
    hdf = I3HDFTableService(outFile)

    #====================================================================
    # Clean up events

    tray.AddModule(PruneIceTop,
            it_stream = 'top_hlc_clusters')

    tray.AddModule(moveMCPrimary)


    #====================================================================
    # Cut information

    tray.AddModule(GetStations,
            InputITpulses = recoPulses,
            output = 'NStations')

    tray.AddModule(FindLoudestStation,
            InputITpulses = recoPulses,
            SaturationValue = 600,
            output = 'SaturationList')

    tray.AddModule(LoudestStationOnEdge,
            InputLoudestStation = 'LoudestStation',
            config = config,
            output = 'LoudestOnEdge')

    tray.AddModule(LoudestStationOnEdge,
            InputLoudestStation = 'SaturationList',
            config = config,
            output = 'SaturatedOnEdge')

    tray.AddModule(LargestTankCharges,
            ITpulses = recoPulses)


    #====================================================================
    # Run the reconstruction

    tray.AddModule(ShowerLLH.GridLLH_C_ztest, 'ShowerLLHGrid_C_ztest',
            LLHTables = LLHTables,
            binDict = binDict,
            grids = grids, 
            recoPulses = recoPulses)

    #====================================================================
    # Finish

    tray.AddModule(I3TableWriter, tableservice=hdf, keys=keys, 
            SubEventStreams=['top_hlc_clusters'])

    tray.Execute()
    tray.Finish()

    print "Time taken: ", time.time() - t0



if __name__ == "__main__":

    outFile  = sys.argv[1]
    fileList = sys.argv[2:]
    maker(outFile, fileList)







