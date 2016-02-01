#!/usr/bin/env python 

################################################################################
# Runs the grid search on the given files and writes to hdf5 output
################################################################################

from icecube import dataio, dataclasses, toprec, ShowerLLH
from icecube import phys_services as i3ps
from icecube import topeventcleaning as tec
from I3Tray import *
from icecube.tableio import I3TableWriter
from icecube.hdfwriter import I3HDFTableService
import sys, time, pickle, glob
import numpy as np
load('libsmallshower-filter')
LLHTables = ShowerLLH.SupportFunctions.LLHTables

sys.path.append('/home/zgriffith/ShowerLLH/useful')
from LLHFunctions import simList

sys.path.append('/home/santander/ShowerLLH/segments')
from laputop_standard_traysegment import LaputopStandard
from laputop_smallshower_traysegment import LaputopSmallShower

#from laputop_fixed_traysegment import LaputopFixed as LaputopStandard
#from laputop_smallfixed_traysegment import LaputopSmallFixed as LaputopSmallShower
from modules_L3cuts import FindLoudestStation, LoudestStationOnEdge, LargestTankCharges, GetStations

import geo

if __name__ == "__main__":

    # Input parameters
    config = sys.argv[1]    # IT59 or IT73
    sim = sys.argv[2]       # See 'maker' for list of simulations used
    group = sys.argv[3]

    recoPulses = 'CleanedHLCTankPulses'
    #files, name, i3name = simList(config, sim, group=group, startnum=0, testnum=10)
    #files, name, i3name = simList(config, sim, group=group, startnum=0, testnum=1)
    files, name, i3name = simList(config, sim, group=group)

    # Keys to write to frame
    keys =  ['ShowerPlane', 'ShowerPlaneParams']
    #keys += ['ShowerCombined', 'ShowerCombinedParams']
    #keys += ['SmallShowerCombined', 'SmallShowerCombinedParams']
    keys += ['LaputopStandard', 'LaputopStandardParams']
    keys += ['LoudestAverageStation', 'LoudestStation', 'LoudestOnEdge']
    keys += ['SaturationList', 'SaturatedOnEdge']
    keys += ['Q1', 'Q2', 'Q3', 'Q4']
    keys += [recoPulses, 'NStations']
    keys += ['ShowerLLH_proton', 'ShowerLLH_iron', 'ShowerLLHParams']
    keys += ['ShowerLLH_helium', 'ShowerLLH_oxygen']
    keys += ['MCPrimary', 'pulse_x', 'pulse_y', 'pulse_time', 'pulse_charges']
    keys += ['hlc_count', 'd_str', 'pulse_z', 'tres', 'hits', 'qtot', 'd_pulse', 'om_row']

    # Load grids for grid scan, llh tables, and edges for cuts
    gridList = [[],[],[]]
    resourcedir = '/net/user/santander/ShowerLLH/resources/'
    fl = open(resourcedir+'BinSpots.pkl', 'rb')
    grids = pickle.load(fl)
    fl.close()

    gridList[0] = grids[config]['coarse']
    gridList[1] = grids[config]['middle']
    gridList[2] = grids[config]['fine']
    LLHlist = LLHTables(resourcedir + config+'_LLHTables.hdf5')

    def getRecoShowers(frame):
      return (frame['I3EventHeader'].sub_event_stream == 'merged_ii_it')

    def moveSmall(frame):
        if frame['SmallShowerFilterPassed'].value:
            frame.Delete('LaputopStandard')
            frame.Delete('LaputopParams')
            frame['LaputopStandard'] = frame['LaputopSmallShower']
            frame['LaputopParams'] = frame['LaputopSmallShowerParams']

    def PruneIceTop(frame):
        if frame['I3EventHeader'].sub_event_stream != 'top_hlc_clusters':
            frame.Delete('ShowerCOG')
            frame.Delete('ShowerCombined')
            frame.Delete('ShowerCombinedParams')
            frame.Delete('ShowerPlane')
            frame.Delete('ShowerPlaneParams')

    icgeo = geo.loadGeometryTxt("/net/user/santander/gamma/geometry/IC79_Outline.dat")
    itgeo = geo.loadGeometryTxt("/net/user/santander/gamma/geometry/IT73_Outline.dat")

    def isCoincidentEvent(frame):
      if frame.Has('Laputop'):
        primary = frame['Laputop']
        partDir = primary.dir
        partPos = primary.pos

        return geo.coincidentCutDirPos(primary, primary, itgeo, icgeo)

      else:
        return False

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

    tray.AddModule('I3IcePickModule<I3SmallShowerFilter>', 'smallfilter',
            FilterGeometry = 'IC79',
            TopPulseKey = recoPulses,
            NStationResultName = 'SmallShowerNStation',
            DecisionName = 'SmallShowerFilterPassed')
    
    #====================================================================
    # Cut information

    tray.AddModule(GetStations, 'NStations',
            InputITpulses = recoPulses,
            OutputName = 'NStations')

    tray.AddModule(FindLoudestStation, "FindLoud",
            InputITpulses = recoPulses,
            SaturationValue = 600,
            OutputSaturatedStationsName = 'SaturationList')

    tray.AddModule(LoudestStationOnEdge, "SaturatedStationOnEdge",
            InputLoudestStation = 'SaturationList',
            OutputName = 'SaturatedOnEdge',
            WhichDetector = config)

    tray.AddModule(LargestTankCharges, "LargestTankCharge",
            ITpulses = recoPulses)

    tray.AddModule(LoudestStationOnEdge, "LoudStationOnEdge",
            InputLoudestStation = 'LoudestStation',
            WhichDetector = config,
            OutputName = 'LoudestOnEdge')


    #====================================================================
    # Run the reconstructions
    

    #tray.AddSegment(LaputopStandard, 'Laputop', 
    #        pulses = recoPulses)

    #tray.AddSegment(LaputopSmallShower, 'LaputopSmall', 
    #        pulses = recoPulses) 
    tray.AddModule(ShowerLLH.GridLLH_C, 'ShowerLLHGrid', 
            LLHTables = LLHlist,
            gridList = gridList, 
            recoPulses = recoPulses)

    tray.AddModule(getRecoShowers, "reccut")       
    
    tray.AddSegment(LaputopStandard, 'LaputopStandard', 
            pulses = recoPulses)
    
    tray.AddSegment(LaputopSmallShower, 'LaputopSmallShower', 
            pulses = recoPulses) 

    tray.AddModule(moveSmall, 'moveSmall')

    
    #tray.AddModule(isCoincidentEvent, "coinc")

    def stats(frame):

        if frame['I3EventHeader'].sub_event_stream == 'merged_ii_it':

            allpulses       = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, 'MaskedOfflinePulses')
            geo             = frame['I3Geometry']
            laputop         = frame['LaputopStandard']
            tres         = []
            time         = []
            charges      = []
            d            = []
            d2           = []
            dpos         = 0
            hlcs         = []
            pulse_x      = []
            pulse_y      = []
            pulse_z      = []
            d_pulse      = []
            om_row       = []
            slcs         = 0
            edge_list = [2,3,4,5,6,8,13,15,21,23,30,32,40,41,50,51,59,60,67,68,72,73,74,75,76,77,78] #IC79 outer strings
            non_list  = [1,7,14,22,31,79,80] #string numbers not in IC79
            top_list  = [1,2,3,4,5,6]

            for i in range(86):
                if i+1 not in non_list:
                    omgeo =  geo.omgeo[icetray.OMKey(i+1,1)]
                    if i+1 not in edge_list:
                        d0       = i3ps.I3Calculator.closest_approach_distance(laputop, omgeo.position)
                        d       += [np.abs(d0)] # closest approach distance to DOM

            for omkey, pulses in allpulses:
                omgeo    = geo.omgeo[omkey]
                #tres    += [i3ps.I3Calculator.time_residual(laputop, omgeo.position, p.time) for p in pulses] # time res wrt Cherenkov cone
                charges += [p.charge for p in pulses]
                pulse_x += [omgeo.position.x for p in pulses]
                pulse_y += [omgeo.position.y for p in pulses]
                pulse_z += [omgeo.position.z for p in pulses]
                d_pulse += [i3ps.I3Calculator.closest_approach_distance(laputop, omgeo.position) for p in pulses]
                hlcs    += [(p.flags & dataclasses.I3RecoPulse.PulseFlags.LC) for p in pulses] # HLC hits
                time    += [p.time for p in pulses]
                om_row  += [omkey[1] for p in pulses]

                if omkey[1] in top_list:
                    for p in pulses:
                        if (p.flags & dataclasses.I3RecoPulse.PulseFlags.LC) == 0:
                            if np.abs(i3ps.I3Calculator.closest_approach_distance(laputop, omgeo.position)) < 200:
                                if p.time > 14800 and p.time < 17500:
                                    slcs += 1

            frame['hits']           = dataclasses.I3VectorInt(hlcs)
            #frame['tres']           = dataclasses.I3VectorDouble(tres)
            frame['pulse_x']        = dataclasses.I3VectorDouble(pulse_x)
            frame['pulse_y']        = dataclasses.I3VectorDouble(pulse_y)
            frame['pulse_z']        = dataclasses.I3VectorDouble(pulse_z)
            frame['hlc_count']      = dataclasses.I3Double(hlcs.count(1))
            frame['d_str']          = dataclasses.I3Double(np.min(d))
            frame['slc_count']      = dataclasses.I3Double(slcs)
            frame['pulse_time']     = dataclasses.I3VectorDouble(time)
            frame['pulse_charges']  = dataclasses.I3VectorDouble(charges)
            frame['d_pulse']        = dataclasses.I3VectorDouble(d_pulse)
            frame['om_row']         = dataclasses.I3VectorInt(om_row)

    tray.AddModule(stats, 'stats')

    def recocut(frame):
        if frame.Has('ShowerLLH_proton'):
            if frame['ShowerLLH_proton'].energy != 0:
                print frame['ShowerLLH_proton'].energy
                return True
        print 'reject'
        return False
    tray.AddModule(recocut, 'recocut')
    #tray.AddModule(lambda frame: frame['ShowerLLHParams'].RecoRanCut == 1, 'recocut')
    #====================================================================
    # Finish

    tray.AddModule(I3TableWriter, 'writer', tableservice=hdf, keys=keys, 
            SubEventStreams=['merged_ii_it'])

    
    #tray.AddModule("Delete", "del",
    #  keys = ["MMCTrackList", "NFE*", "TWNFEMergedPulsesHLC", "I3MCTree"]
    #)
    
    #tray.AddModule("I3Writer", 'i3writer',
    #  filename=i3name,
    #  DropOrphanStreams=[icetray.I3Frame.DAQ]
    #)          

    tray.AddModule('TrashCan', 'yeswecan')
    tray.Execute()
    tray.Finish()

    print "Time taken: ", time.time() - t0

