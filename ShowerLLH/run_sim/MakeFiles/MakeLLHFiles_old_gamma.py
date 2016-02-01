#!/usr/bin/env python 

################################################################################
# Runs the grid search on the given files and writes to hdf5 output
################################################################################
from icecube import icetray, dataio, dataclasses,ShowerLLH, InIceLLH#,toprec
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
#from laputop_fixed_traysegment import LaputopFixed
#from laputop_smallfixed_traysegment import LaputopSmallFixed

from modules import FindLoudestStation, LoudestStationOnEdge, LargestTankCharges
from modules import GetStations, PruneIceTop, moveMCPrimary #,moveSmall
from module import MergeIIIT

def maker(outFile, outFile_i3, fileList):

    # Starting parameters
    pulses_ii = 'OfflinePulses'
    recoPulses = 'IceTopHLCTankPulses'
    track_ii = 'PoleMuonLinefit'
    track_it = 'ShowerPlane'
    config = 'IT73'
    #resourcedir = '/data/user/zgriffith/ShowerLLH/resources/'
    resourcedir  = '/data/user/fmcnally/ShowerLLH/resources/'
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
    keys += ['LaputopStandard', 'LaputopStandardParams']
    keys += ['LoudestStation', 'LoudestOnEdge']
    keys += ['SaturationList', 'SaturatedOnEdge']
    keys += ['I3EventHeader']
    keys += ['Q1', 'Q2', 'Q3', 'Q4']
    keys += [recoPulses, 'NStations']
    for comp in LLHTables.keys():
        keys += ['ShowerLLH_'+compDict[comp]]
        #keys += ['ShowerLLH_lc_'+compDict[comp]]
    keys += ['ShowerLLHParams', 'slc_count', 'mjd_time']
    keys += ['MCPrimary', 'pulse_x', 'pulse_y', 'pulse_time', 'pulse_charges']
    keys += ['hlc_count', 'd_pulse', 'd_str', 'pulse_z', 'tres', 'hits', 'qtot', 'om_row']
    keys += ['MaskedOfflinePulses', 'FilterMask', 'mjd_time']
    t0 = time.time()

    tray = I3Tray()
    tray.AddModule('I3Reader', FileNameList=fileList)
    hdf = I3HDFTableService(outFile)

    def moveSmall(frame):
        isSmall = False
        if frame.Has('IsSmallShower'):
            if frame['IsSmallShower'].value == True:
                isSmall = True
        elif frame.Has('SmallShowerFilterPassed'):
            isSmall = True
        if isSmall:
            frame.Delete('LaputopStandard')
            frame.Delete('LaputopStandardParams')
            frame['LaputopStandard'] = frame['LaputopSmallShower']
            frame['LaputopStandardParams'] = frame['LaputopSmallParams']

    #====================================================================
    # Clean up events

    tray.AddModule(PruneIceTop,
            it_stream = 'top_hlc_clusters')

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

    tray.AddSegment(LaputopStandard, 'LaputopStandard',
            pulses = recoPulses)

    tray.AddSegment(LaputopSmallShower, 'LaputopSmallShower',
            pulses = recoPulses) 

    tray.AddModule(moveSmall)

    tray.AddModule(ShowerLLH.GridLLH_C, 'ShowerLLHGrid_C',
            LLHTables = LLHTables,
            binDict = binDict,
            grids = grids,
            recoPulses = recoPulses)
    #tray.AddModule(ShowerLLH.GridLLH_C_laputop_core, 'ShowerLLHGrid_C_laputop_core',
    #        LLHTables = LLHTables,
    #        binDict = binDict,
    #        grids = grids,
    #        recoPulses = recoPulses) 
    
   #=====================================================================
    # Add Calculated Parameters

    def stats(frame):

        if frame['I3EventHeader'].sub_event_stream == 'merged_ii_it':

            allpulses       = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, 'OfflinePulses')
            geo             = frame['I3Geometry']
            laputop         = frame['LaputopStandard']
            tres         = []
            time         = []
            charges      = []
            d            = []
            d2           = []
            dpos         = []
            hlcs         = []
            pulse_x      = []
            pulse_y      = []
            pulse_z      = []
            d_pulse      = []
            om_row       = []
            slcs         = 0
            edge_list = [1,2,3,4,5,6,7,13,14,21,22,30,31,40,41,50,51,59,60,67,68,72,73,74,75,76,77,78]
            top_list  = [1,2,3,4,5,6]

            #for i in range(86):
            #    omgeo =  geo.omgeo[icetray.OMKey(i+1,1)]
            #    if i+1 not in edge_list:
            #        d0       = i3ps.I3Calculator.closest_approach_distance(laputop, omgeo.position)
            #        d       += [np.abs(d0)] # closest approach distance to DOM

            for omkey, pulses in allpulses:
                omgeo    = geo.omgeo[omkey]
                tres    += [i3ps.I3Calculator.time_residual(laputop, omgeo.position, p.time) for p in pulses] # time res wrt Cherenkov cone
                charges += [p.charge for p in pulses]
                pulse_x += [omgeo.position.x for p in pulses]
                pulse_y += [omgeo.position.y for p in pulses]
                pulse_z += [omgeo.position.z for p in pulses]
                d_pulse += [i3ps.I3Calculator.closest_approach_distance(laputop, omgeo.position) for p in pulses]
                hlcs    += [(p.flags & dataclasses.I3RecoPulse.PulseFlags.LC) for p in pulses] # HLC hits
                time    += [p.time for p in pulses]
                d_close  = i3ps.I3Calculator.closest_approach_distance(laputop, omgeo.position)
                om_row  += [omkey[1] for p in pulses]

                if omkey[1] in top_list:
                    for p in pulses:
                        if (p.flags & dataclasses.I3RecoPulse.PulseFlags.LC) == 0:
                            if np.abs(d_close) < 200 and p.time > 14800 and p.time < 17500:
                                slcs += 1
            '''
            frame['hits']           = dataclasses.I3VectorInt(hlcs)
            frame['tres']           = dataclasses.I3VectorDouble(tres)
            frame['d_pulse']        = dataclasses.I3VectorDouble(d_pulse)
            frame['pulse_x']        = dataclasses.I3VectorDouble(pulse_x)
            frame['pulse_y']        = dataclasses.I3VectorDouble(pulse_y)
            frame['pulse_z']        = dataclasses.I3VectorDouble(pulse_z)
            frame['pulse_time']     = dataclasses.I3VectorDouble(time)
            frame['pulse_charges']  = dataclasses.I3VectorDouble(charges)
            frame['om_row']         = dataclasses.I3VectorInt(om_row)
            '''

            frame['mjd_time']       = dataclasses.I3Double(frame['I3EventHeader'].start_time.mod_julian_day_double)
            frame['hlc_count']      = dataclasses.I3Double(hlcs.count(1))
            #frame['d_str']          = dataclasses.I3Double(np.min(d))
            frame['slc_count']      = dataclasses.I3Double(slcs)
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







