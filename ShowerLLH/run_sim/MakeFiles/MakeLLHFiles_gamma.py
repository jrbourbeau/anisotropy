#!/usr/bin/env python 

################################################################################
# Runs the grid search on the given files and writes to hdf5 output
################################################################################
from icecube import dataio, dataclasses, toprec, ShowerLLH, InIceLLH
from icecube import phys_services as i3ps
from I3Tray import *
from icecube.tableio import I3TableWriter
from icecube.hdfwriter import I3HDFTableService

import numpy as np
import sys, time, glob

#load('libsmallshower-filter')
sys.path.append('/home/zgriffith/ShowerLLH/segments')
sys.path.append('/home/zgriffith/ShowerLLH/useful')
import geo
from laputop_standard_traysegment import LaputopStandard
from laputop_smallshower_traysegment import LaputopSmallShower
from modules import FindLoudestStation, LoudestStationOnEdge, LargestTankCharges
from modules import GetStations, moveSmall, PruneIceTop, moveMCPrimary
from module import MergeIIIT

def maker(outFile, outFile_i3, fileList):

    # Starting parameters
    pulses_ii = 'OfflinePulses'
    recoPulses = 'CleanedHLCTankPulses'
    track_ii = 'PoleMuonLlhFit'
    track_it = 'ShowerPlane'
    config = 'IT81'
    #resourcedir = '/net/user/zgriffith/ShowerLLH/resources/'
    resourcedir  = '/net/user/fmcnally/ShowerLLH/resources/'
    compDict = {'P':'proton', 'He':'helium', 'O':'oxygen', 'Fe':'iron', 'G':'gamma'}

    
    i3name = outFile_i3

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
    keys += ['I3EventHeader']
    keys += ['Q1', 'Q2', 'Q3', 'Q4']
    keys += [recoPulses, 'NStations']
    for comp in LLHTables.keys():
        keys += ['ShowerLLH_'+compDict[comp]]
    keys += ['ShowerLLHParams']
    keys += ['MCPrimary', 'pulse_time', 'pulse_charges']
    keys += ['hlc_count', 'd_str', 'pulse_z', 'tres', 'hits']
    keys += ['MaskedOfflinePulses', 'FilterMask', 'prescale_event']
    t0 = time.time()

    tray = I3Tray()
    tray.AddModule('I3Reader', FileNameList=fileList)
    hdf = I3HDFTableService(outFile)
    #====================================================================
    # Clean up events

    tray.AddModule(PruneIceTop,
            it_stream = 'ice_top')

    '''
    tray.AddModule(tec.modules.MergeIIIT, 'merged_ii_it',
            IceTopReco='ShowerCombined',
            InIceReco = 'LineFit')
    '''

    tray.AddModule(MergeIIIT, 'merged_ii_it',
        IceTopReco = track_it,          # Previously ShowerCombined
        InIceReco  = track_ii,          # Previously LineFit
        IceTopPulses = recoPulses,
        InIcePulses  = pulses_ii)

    tray.AddModule(lambda frame:frame['I3EventHeader'].sub_event_stream == 'merged_ii_it' or frame['I3EventHeader'].sub_event_stream == 'icetop_only')
    #tray.AddModule(lambda frame:frame.Has('Laputop'))

    #tray.AddModule(moveMCPrimary)

    """
    tray.AddModule('I3IcePickModule<I3SmallShowerFilter>',
            FilterGeometry = 'IC86',
            TopPulseKey = recoPulses,
            NStationResultName = 'SmallShowerNStation',
            DecisionName = 'SmallShowerFilterPassed')
    """

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
    '''
    tray.AddSegment(LaputopStandard, 'Laputop',
            pulses = recoPulses)

    tray.AddSegment(LaputopSmallShower, 'LaputopSmall',
            pulses = recoPulses) 

    tray.AddModule(moveSmall)
    '''
    tray.AddModule(ShowerLLH.GridLLH_C, 'ShowerLLHGrid_C',
            LLHTables = LLHTables,
            binDict = binDict,
            grids = grids, 
            recoPulses = recoPulses)
    
   #=====================================================================
    # Add Calculated Parameters

    def stats(frame):
        if frame['I3EventHeader'].sub_event_stream == 'merged_ii_it':

            allpulses       = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, 'MaskedOfflinePulses')
            geo             = frame['I3Geometry']
            laputop         = frame['Laputop']
            tres         = []
            time         = []
            charges      = []
            d            = []
            d2           = []
            dpos         = []
            hlcs         = []
            pulse_z      = []
            edge_list = [1,2,3,4,5,6,7,13,14,21,22,30,31,40,41,50,51,59,60,67,68,72,73,74,75,76,77,78]

            for omkey, pulses in allpulses:
                omgeo    = geo.omgeo[omkey]
                if omkey[0] not in edge_list:
                    d0       = i3ps.I3Calculator.closest_approach_distance(laputop, omgeo.position)
                    d2      += [d0*d0]
                    d       += [d0] # closest approach distance to DOM
                tres    += [i3ps.I3Calculator.time_residual(laputop, omgeo.position, p.time) for p in pulses] # time res wrt Cherenkov cone
                charges += [p.charge for p in pulses]
                pulse_z += [omgeo.position.z for p in pulses]
                dpos    += [i3ps.I3Calculator.closest_approach_position(laputop, omgeo.position)]
                hlcs    += [(p.flags & dataclasses.I3RecoPulse.PulseFlags.LC) for p in pulses] # HLC hits
                time    += [p.time for p in pulses]

            frame['hits']           = dataclasses.I3VectorInt(hlcs)
            frame['tres']           = dataclasses.I3VectorDouble(tres)
            frame['pulse_z']        = dataclasses.I3VectorDouble(pulse_z)
            frame['hlc_count']      = dataclasses.I3Double(hlcs.count(1))
            frame['d_str']          = dataclasses.I3Double(np.min(d))
            frame['pulse_time']     = dataclasses.I3VectorDouble(time)
            frame['pulse_charges']  = dataclasses.I3VectorDouble(charges)
            frame['prescale_event'] = dataclasses.I3VectorBool([frame['FilterMask']['FilterMinBias_11'].prescale_passed])
    tray.AddModule(stats, 'stats')






    #====================================================================
    # Finish
    tray.AddModule(I3TableWriter, tableservice=hdf, keys=keys, SubEventStreams = ['merged_ii_it','icetop_only'])
    #tray.AddModule(I3TableWriter, tableservice=hdf_it_only, keys=keys, SubEventStreams = ['icetop_only'])

    #tray.AddModule(lambda frame:frame['I3EventHeader']=='merged_ii_it')
    #tray.AddModule(lambda frame:frame['I3EventHeader']=='icetop_only')
    #tray.AddModule("I3Writer", 'i3writer',
    #  DropOrphanStreams=[icetray.I3Frame.DAQ],
    #  filename=i3name)
  

    tray.Execute()
    tray.Finish()

    print "Time taken: ", time.time() - t0



if __name__ == "__main__":

    outFile  = sys.argv[1]
    outFile_i3 = sys.argv[2]
    fileList = sys.argv[3:]
    maker(outFile, outFile_i3,fileList)







