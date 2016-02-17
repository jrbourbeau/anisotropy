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
from laputop.laputop_standard_traysegment import LaputopStandard
from laputop.laputop_smallshower_traysegment import LaputopSmallShower
from i3modules import GetStations, moveSmall, PruneIceTop
from i3modules import FindLoudestStation, LoudestStationOnEdge
from i3modules import LargestTankCharges

def maker(config, outFile, fileList):

    # Starting parameters
    recoPulses = 'CleanedHLCTankPulses'
    track = 'ShowerPlane'
    resourcedir = '/net/user/fmcnally/ShowerLLH/resources/'
    compDict = {'P':'proton', 'He':'helium', 'O':'oxygen', 'Fe':'iron'}

    if config == 'IT73':
        it_stream = 'top_hlc_clusters'
        ITgeo = 'IT73'
        ICgeo = 'IC79'
        filter_mask = 'FilterMask'
    if config == 'IT81':
        it_stream = 'ice_top'
        ITgeo = 'IT81'
        ICgeo = 'IC86'
        filter_mask = 'FilterMask'
    if config == 'IT81-II':
        it_stream = 'ice_top'
        ITgeo = 'IT81'
        ICgeo = 'IC86'
        filter_mask = 'QFilterMask'

    # Load grids for grid scan, llh tables, and bins
    grids = np.load(resourcedir + ITgeo+'_BinSpots.npy')
    grids = grids.item()
    LLHTables = np.load(resourcedir + 'LLHTables.npy')
    LLHTables = LLHTables.item()
    binDict = np.load(resourcedir + 'ShowerLLH_bins.npy')
    binDict = binDict.item()

    # Keys to write to frame
    keys  = []
    keys += ['I3EventHeader', filter_mask]
    keys += ['ShowerPlane', 'ShowerPlaneParams']
    keys += ['LaputopStandard', 'LaputopStandardParams']
    keys += ['NStations']
    for comp in LLHTables.keys():
        keys += ['ShowerLLH_'+compDict[comp]]
    keys += ['ShowerLLHParams']

    t0 = time.time()

    tray = I3Tray()
    tray.AddModule('I3Reader', FileNameList=fileList)
    hdf = I3HDFTableService(outFile)

    #====================================================================
    # Clean up events

    tray.AddModule(PruneIceTop, 
            it_stream = it_stream)

    tray.AddModule('I3IcePickModule<I3SmallShowerFilter>',
            FilterGeometry = ICgeo,
            TopPulseKey = recoPulses,
            NStationResultName = 'SmallShowerNStation',
            DecisionName = 'SmallShowerFilterPassed')


    #====================================================================
    # ShowerLLH cuts

    tray.AddModule(GetStations,
            InputITpulses = recoPulses,
            output = 'NStations')

    tray.AddModule(FindLoudestStation,
            InputITpulses = recoPulses,
            SaturationValue = 600,
            output = 'SaturationList')

    tray.AddModule(LoudestStationOnEdge,
            InputLoudestStation = 'LoudestStation',
            config = ITgeo,
            output = 'LoudestOnEdge')

    tray.AddModule(LargestTankCharges,
            ITpulses = recoPulses)

    def IceTop_cuts(frame):
        if frame.Has(track):
            c1 = (np.cos(frame[track].dir.zenith) >= 0.8)   # zenith cut
            c2 = (not frame['LoudestOnEdge'].value)         # edge cut
            c3 = (frame['Q1'].value >= 6)                   # charge cut
            return c1*c2*c3

    tray.AddModule(IceTop_cuts)


    #====================================================================
    # Run the reconstructions

    # Laputop already run in IT81-II level2
    if config != 'IT81-II':
        tray.AddSegment(LaputopStandard, 'LaputopStandard',
                pulses = recoPulses)
        tray.AddSegment(LaputopSmallShower, 'LaputopSmallShower',
                pulses = recoPulses)

    tray.AddModule(moveSmall)

    tray.AddModule(ShowerLLH.GridLLH_C, 'ShowerLLHGrid',
            LLHTables = LLHTables,
            binDict = binDict,
            grids = grids, 
            recoPulses = recoPulses)

    #====================================================================
    # Finish

    tray.AddModule(I3TableWriter, tableservice=hdf, keys=keys, 
            SubEventStreams=[it_stream])

    tray.Execute()
    tray.Finish()

    print "Time taken: ", time.time() - t0


if __name__ == "__main__":

    config   = sys.argv[1]
    outFile  = sys.argv[2]
    fileList = sys.argv[3:]
    maker(config, outFile, fileList)



