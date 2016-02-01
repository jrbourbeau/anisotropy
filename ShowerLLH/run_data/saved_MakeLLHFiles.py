#!/usr/bin/env python 

################################################################################
# Runs the grid search on the given files and writes to hdf5 output
################################################################################

from icecube import dataclasses,dataio, toprec, ShowerLLH
from icecube import phys_services as i3ps
from I3Tray import *
from icecube.tableio import I3TableWriter
from icecube.hdfwriter import I3HDFTableService

import numpy as np
import sys, time, glob

load('libsmallshower-filter')
sys.path.append('/home/zgriffith/ShowerLLH/segments')
from laputop_snowcorr_traysegment import LaputopStandard
from laputop_smallshower_snowcorr_traysegment import LaputopSmallShower
from modules import FindLoudestStation, LoudestStationOnEdge, LargestTankCharges
from modules import GetStations, PruneIceTop #, moveSmall

sys.path.append('/home/zgriffith/ShowerLLH/useful')
import geo as geo_functions

def maker(config, outFile, fileList):

    # Starting parameters
    recoPulses = 'CleanedHLCTankPulses'
    track = 'ShowerPlane'
    resourcedir = '/data/user/fmcnally/ShowerLLH/resources/'
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
    keys += ['Q1', 'Q2', 'Q3', 'Q4']
    keys += [recoPulses, 'NStations']
    for comp in LLHTables.keys():
        keys += ['ShowerLLH_'+compDict[comp]]
    keys += ['ShowerLLHParams']
    keys += ['hlc_count', 'slc_count', 'hlc_noise', 'hlc_signal', 'd_str', 'mjd_time']

    #keys recorded only if you want individual hit information
    keys += ['pulse_x', 'pulse_y', 'pulse_z']
    keys += ['hits', 'tres', 'pulse_time', 'pulse_charges', 'd_pulse']

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
            frame['LaputopStandardParams'] = frame['LaputopSmallShowerParams']

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
    def IceTop_cuts(frame):
        if frame.Has(track):
            c1 = (np.cos(frame[track].dir.zenith) >= 0.8)   # zenith cut
            #c2 = (not frame['LoudestOnEdge'].value)         # edge cut
            #c3 = (frame['Q1'].value >= 6)                   # charge cut
            return c1

    #tray.AddModule(IceTop_cuts)
    tray.AddModule(GetStations,
        InputITpulses = recoPulses,
        output = 'NStations')

    #====================================================================
    # Run the reconstructions
    tray.AddModule(ShowerLLH.GridLLH_C, 'ShowerLLHGrid',
            LLHTables = LLHTables,
            binDict = binDict,
            grids = grids, 
            recoPulses = recoPulses)

    # Laputop already run in IT81-II level2
    
    if config != 'IT81-II':
        tray.AddSegment(LaputopStandard, 'LaputopStandard',
                pulses = recoPulses)
        tray.AddSegment(LaputopSmallShower, 'LaputopSmallShower',
                pulses = recoPulses)

    tray.AddModule(moveSmall)

    #====================================================================
    #In Ice Pulse Information
    def pulse_info(frame):

        allpulses       = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, 'InIcePulses')
        geo             = frame['I3Geometry']
        laputop         = frame['LaputopStandard']
        tres         = []
        time         = []
        charges      = []
        d            = []
        d2           = []
        charges      = []
        dpos         = []
        hlcs         = []
        pulse_x      = []
        pulse_y      = []
        pulse_z      = []
        d_pulse      = []
        slcs         = 0
        hlc_noise    = 0
        hlc_signal   = 0
        edge_list = [1,2,3,4,5,6,7,13,14,21,22,30,31,40,41,50,51,59,60,67,68,72,73,74,75,76,77,78]
        top_list  = [1,2,3,4,5,6]


        for i in range(78):
            omgeo =  geo.omgeo[icetray.OMKey(i+1,1)]
            if i+1 not in edge_list:
                d0  = geo_functions.distanceToInnerString(laputop.pos.x, laputop.pos.y, laputop.dir.zenith, laputop.dir.azimuth, omgeo.position[0], omgeo.position[1])
                d  += [np.abs(d0)] # closest approach distance to DOM
 
        for omkey, pulses in allpulses:
            omgeo    = geo.omgeo[omkey]
            d_pulse += [i3ps.I3Calculator.closest_approach_distance(laputop, omgeo.position) for p in pulses]
            tres    += [i3ps.I3Calculator.time_residual(laputop, omgeo.position, p.time) for p in pulses] # time res wrt Cherenkov cone
            charges += [p.charge for p in pulses]
            pulse_x += [omgeo.position.x for p in pulses]
            pulse_y += [omgeo.position.y for p in pulses]
            pulse_z += [omgeo.position.z for p in pulses]
            dpos    += [i3ps.I3Calculator.closest_approach_position(laputop, omgeo.position)]
            hlcs    += [(p.flags & dataclasses.I3RecoPulse.PulseFlags.LC) for p in pulses] # HLC hits
            time    += [p.time for p in pulses]
            d_close  = i3ps.I3Calculator.closest_approach_distance(laputop, omgeo.position)
            for p in pulses:
                if (p.flags & dataclasses.I3RecoPulse.PulseFlags.LC) != 0:
                    if np.abs(d0) > 200 or p.time < 14800 or p.time > 17500:
                        hlc_noise  += 1
                    if np.abs(d0) < 200 and p.time > 14800 and p.time < 17500:
                        hlc_signal += 1

            if omkey[1] in top_list:
                for p in pulses:
                    if (p.flags & dataclasses.I3RecoPulse.PulseFlags.LC) == 0:
                        if np.abs(d0) < 200 and p.time > 14800 and p.time < 17500:
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
        '''
        frame['hlc_count']      = dataclasses.I3Double(hlcs.count(1))
        frame['slc_count']      = dataclasses.I3Double(slcs)
        frame['hlc_noise']      = dataclasses.I3Double(hlc_noise)
        frame['hlc_signal']      = dataclasses.I3Double(hlc_signal)
        frame['mjd_time']       = dataclasses.I3Double(frame['I3EventHeader'].start_time.mod_julian_day_double)
        frame['d_str']          = dataclasses.I3Double(np.min(d))
    tray.AddModule(pulse_info, 'PulseInfo')
    #====================================================================
    # Finish


    tray.AddModule(I3TableWriter, tableservice=hdf, keys=keys, 
            SubEventStreams=[it_stream, "nullsplit"])

    tray.Execute()
    tray.Finish()

    print "Time taken: ", time.time() - t0


if __name__ == "__main__":

    config   = sys.argv[1]
    outFile  = sys.argv[2]
    fileList = sys.argv[3:]
    maker(config, outFile, fileList)



