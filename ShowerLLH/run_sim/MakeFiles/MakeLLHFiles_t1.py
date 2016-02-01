#!/usr/bin/env python 

################################################################################
# Runs the grid search on the given files and writes to hdf5 output
################################################################################

from icecube import dataio, toprec, ShowerLLH
from icecube import topeventcleaning as tec
from I3Tray import *
from icecube.tableio import I3TableWriter
from icecube.hdfwriter import I3HDFTableService
import sys, time, pickle
#load('libsmallshower-filter')
LLHTables = ShowerLLH.SupportFunctions.LLHTables

sys.path.append('/home/fmcnally/ShowerLLH/useful')
from LLHFunctions import simList

sys.path.append('/home/fmcnally/ShowerLLH/segments')
#from laputop_standard_traysegment import LaputopStandard
#from laputop_smallshower_traysegment import LaputopSmallShower
#from modules_L3cuts import FindLoudestStation, LoudestStationOnEdge, LargestTankCharges, GetStations

if __name__ == "__main__":

    # Input parameters
    #config = sys.argv[1]    # IT59 or IT73
    #sim = sys.argv[2]       # See 'maker' for list of simulations used
    #group = sys.argv[3]

    config = 'IT73'
    sim = '7273'
    startnum = 1
    testnum = 10

    recoPulses = 'CleanedHLCTankPulses'
    files, name = simList(config, sim, startnum=startnum, testnum=testnum)
    name = '/home/fmcnally/ShowerLLH/run_sim/test.hdf5'

    # Keys to write to frame
    keys =  ['ShowerPlane', 'ShowerPlaneParams']
    #keys += ['ShowerCombined', 'ShowerCombinedParams']
    #keys += ['SmallShowerCombined', 'SmallShowerCombinedParams']
    keys += ['Laputop', 'LaputopParams']
    keys += ['LoudestAverageStation', 'LoudestStation', 'LoudestOnEdge']
    keys += ['SaturationList', 'SaturatedOnEdge']
    keys += ['Q1', 'Q2', 'Q3', 'Q4']
    keys += ['NStations']
    keys += [recoPulses]
    keys += ['ShowerLLH_proton', 'ShowerLLH_iron', 'ShowerLLHParams']
    keys += ['ShowerLLH_helium', 'ShowerLLH_oxygen']
    keys += ['MCPrimary']

    # Load grids for grid scan, llh tables, and edges for cuts
    gridList = [[],[],[]]
    resourcedir = '/net/user/fmcnally/ShowerLLH/resources/'
    fl = open(resourcedir+'BinSpots.pkl', 'rb')
    grids = pickle.load(fl)
    fl.close()

    gridList[0] = grids[config]['coarse']
    gridList[1] = grids[config]['middle']
    gridList[2] = grids[config]['fine']
    LLHlist = LLHTables(resourcedir + config+'_LLHTables.hdf5')

    def moveSmall(frame):
        if frame['SmallShowerFilterPassed'].value:
            frame.Delete('Laputop')
            frame.Delete('LaputopParams')
            frame['Laputop'] = frame['LaputopSmall']
            frame['LaputopParams'] = frame['LaputopSmallParams']

    def PruneIceTop(frame):
        if frame['I3EventHeader'].sub_event_stream != 'top_hlc_clusters':
            frame.Delete('ShowerCOG')
            frame.Delete('ShowerCombined')
            frame.Delete('ShowerCombinedParams')
            frame.Delete('ShowerPlane')
            frame.Delete('ShowerPlaneParams')

    t0 = time.time()

    tray = I3Tray()
    tray.AddModule('I3Reader', 'reader', FileNameList=files)
    hdf = I3HDFTableService(name)

    #====================================================================
    # Clean up events

    tray.AddModule(PruneIceTop, 'pruner')

    tray.AddModule(tec.modules.KeepOnlyLargestEvent, 'keep_only_largest_event',
            Pulses = recoPulses)

    tray.AddModule(tec.modules.MergeIIIT, 'merged_ii_it',
            IceTopReco='ShowerCombined',
            InIceReco = 'LineFit')

    #tray.AddModule('I3IcePickModule<I3SmallShowerFilter>', 'smallfilter',
    #        FilterGeometry = 'IC79',
    #        TopPulseKey = recoPulses,
    #        NStationResultName = 'SmallShowerNStation',
    #        DecisionName = 'SmallShowerFilterPassed')

    #====================================================================
    # Cut information

    #tray.AddModule(GetStations, 'NStations',
    #        InputITpulses = recoPulses,
    #        OutputName = 'NStations')
    
    #tray.AddModule(FindLoudestStation, "FindLoud",
    #        InputITpulses = recoPulses,
    #        SaturationValue = 600,
    #        OutputSaturatedStationsName = 'SaturationList')

    #tray.AddModule(LoudestStationOnEdge, "SaturatedStationOnEdge",
    #        InputLoudestStation = 'SaturationList',
    #        OutputName = 'SaturatedOnEdge',
    #        WhichDetector = config)

    #tray.AddModule(LargestTankCharges, "LargestTankCharge",
    #        ITpulses = recoPulses)

    #tray.AddModule(LoudestStationOnEdge, "LoudStationOnEdge",
    #        InputLoudestStation = 'LoudestStation',
    #        WhichDetector = config,
    #        OutputName = 'LoudestOnEdge')


    #====================================================================
    # Run the reconstructions

    #tray.AddSegment(LaputopStandard, 'Laputop', 
    #        pulses = recoPulses)

    #tray.AddSegment(LaputopSmallShower, 'LaputopSmall', 
    #        pulses = recoPulses) 

    #tray.AddModule(moveSmall, 'moveSmall')

    tray.AddModule(ShowerLLH.GridLLH_C, 'ShowerLLHGrid_C', 
            LLHTables = LLHlist,
            gridList = gridList, 
            recoPulses = recoPulses)

    tray.AddModule(ShowerLLH.GridLLH_fast, 'ShowerLLHGrid_fast', 
            LLHTables = LLHlist,
            gridList = gridList, 
            recoPulses = recoPulses)
 
    #====================================================================
    # Finish

    tray.AddModule(I3TableWriter, 'writer', tableservice=hdf, keys=keys, 
            SubEventStreams=['merged_ii_it'])

    tray.AddModule('TrashCan', 'yeswecan')
    tray.Execute()
    tray.Finish()

    print "Time taken: ", time.time() - t0

