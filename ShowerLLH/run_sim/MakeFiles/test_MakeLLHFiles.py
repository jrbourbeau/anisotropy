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

load('libsmallshower-filter')
sys.path.append('/home/fmcnally/ShowerLLH/segments')
from laputop_standard_traysegment import LaputopStandard
from laputop_smallshower_traysegment import LaputopSmallShower
from modules import FindLoudestStation, LoudestStationOnEdge, LargestTankCharges
from modules import GetStations, moveSmall, PruneIceTop, moveMCPrimary

def maker(config, outFile, fileList):

    # Starting parameters
    recoPulses = 'CleanedHLCTankPulses'
    resourcedir = '/net/user/fmcnally/ShowerLLH/resources/'
    compDict = {'P':'proton', 'He':'helium', 'O':'oxygen', 'Fe':'iron'}
    compDict['G'] = 'gamma'
    geoDict = {'IT59':'IC59','IT73':'IC79','IT81':'IC86','IT81-II':'IC86'}
    icgeo = geoDict[config]

    # Load grids for grid scan, llh tables, and bins
    grids = np.load(resourcedir + config+'_BinSpots.npy')
    grids = grids.item()
    LLHTables = np.load(resourcedir + 'LLHTables.npy')
    LLHTables = LLHTables.item()
    binDict = np.load(resourcedir + 'ShowerLLH_bins.npy')
    binDict = binDict.item()

    # Keys to write to frame
    keys  = []
    keys += ['ShowerPlane', 'ShowerPlaneParams']
    keys += ['Laputop', 'LaputopParams']
    keys += ['LoudestStation', 'LoudestOnEdge']
    keys += ['SaturationList', 'SaturatedOnEdge']
    keys += ['Q1', 'Q2', 'Q3', 'Q4']
    keys += [recoPulses, 'NStations']
    for e in LLHTables.keys():
        keys += ['ShowerLLH_'+compDict[e], 'ShowerLLHParams_'+compDict[e]]
    #keys += ['ShowerLLHParams']
    keys += ['MCPrimary']


    t0 = time.time()

    tray = I3Tray()
    tray.AddModule('I3Reader', FileNameList=fileList)
    hdf = I3HDFTableService(outFile)

    #====================================================================
    # Clean up events

    tray.AddModule(PruneIceTop, it_stream = 'ice_top')

    tray.AddModule(moveMCPrimary)

    tray.AddModule('I3IcePickModule<I3SmallShowerFilter>',
            FilterGeometry = icgeo,
            TopPulseKey = recoPulses,
            NStationResultName = 'SmallShowerNStation',
            DecisionName = 'SmallShowerFilterPassed')


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
    # Run the reconstructions

    tray.AddSegment(LaputopStandard, 'Laputop',
            pulses = recoPulses)

    tray.AddSegment(LaputopSmallShower, 'LaputopSmall',
            pulses = recoPulses) 

    #tray.AddModule(moveSmall)

    tray.AddModule(ShowerLLH.GridLLH_C, 'ShowerLLHGrid_C',
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

    config   = sys.argv[1]
    outFile  = sys.argv[2]
    fileList = sys.argv[3:]
    maker(config, outFile, fileList)







