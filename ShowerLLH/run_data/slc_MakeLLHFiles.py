#!/usr/bin/env python 

################################################################################
# Runs the grid search on the given files and writes to hdf5 output
################################################################################
from icecube import icetray, dataio, dataclasses,ShowerLLH, InIceLLH, toprec
from icecube import phys_services as i3ps
from icecube import topeventcleaning as tec
from I3Tray import *
from icecube.tableio import I3TableWriter
from icecube.hdfwriter import I3HDFTableService
'''
import traceback

class TracePrints(object):
  def __init__(self):    
    self.stdout = sys.stdout
  def write(self, s):
    self.stdout.write("Writing %r\n" % s)
    traceback.print_stack(file=self.stdout)

sys.stdout = TracePrints()
'''

import numpy as np
import sys, time, glob, simFunctions
import argparse

load('libsmallshower-filter')
sys.path.append('/home/zgriffith/ShowerLLH/segments')
sys.path.append('/home/zgriffith/ShowerLLH/useful')
import geo as geo_functions

from laputop_traysegment import Laputop
from laputop_smallshower_traysegment import LaputopSmallShower

from modules import FindLoudestStation, LoudestStationOnEdge, LargestTankCharges
from modules import GetStations, PruneIceTop, moveMCPrimary, moveSmall

if __name__ == "__main__":

    p = argparse.ArgumentParser(
            description='Runs ShowerLLH over a given fileList')
    p.add_argument('-f', '--files', dest='files', nargs='*',
            help='Files to run over')
    p.add_argument('-c', '--config', dest='config',
            help='Detector configuration')
    p.add_argument('--gridFile', dest='gridFile',
            help='File containing locations for iterative grid search')
    p.add_argument('--llhFile', dest='llhFile',
            help='File with llh tables for reconstruction')
    p.add_argument('-o', '--outFile', dest='outFile',
            help='Output file')
    args = p.parse_args()

    # Starting parameters
    recoPulses  = 'CleanedHLCTankPulses'
    track_it    = 'ShowerPlane'
    resourcedir = '/data/user/zgriffith/ShowerLLH/resources/'

    if args.config == 'IT73':
        it_stream = 'top_hlc_clusters'
        ICgeo = 'IC79'
    if args.config in ['IT81','IT81-II']:
        it_stream = 'ice_top'
        ITgeo = 'IT81'
        ICgeo = 'IC86'
        filter_mask = 'QFilterMask'

    # Load grids for grid scan, llh tables, and bins
    grids = np.load(args.gridFile)
    grids = grids.item()
    LLHTables = np.load(args.llhFile)
    LLHTables = LLHTables.item()
    binDict = LLHTables['bins']
    LLHTables = LLHTables['llhtables']

    # Keys to write to frame
    keys  = ['I3EventHeader']
    for comp in LLHTables.keys():
        keys += ['ShowerLLH_'+comp, 'maxLLH_'+comp]
    keys += ['ShowerLLHParams']
    keys += ['LaputopStandard', 'LaputopStandardParams']
    keys += ['hlc_count', 'slc_count','mjd_time']
    #keys recorded only if you want individual hit information
    keys += ['pulse_x', 'pulse_y'] 
    keys += ['pulse_z', 'pulse_distance', 'pulse_time']
    keys += ['hits', 'pulse_charges']
    keys += ['off_slc_count']
    #keys += [recoPulses, pulses_ii]

    t0 = time.time()

    tray = I3Tray()
    tray.AddModule('I3Reader', FileNameList=args.files)
    hdf = I3HDFTableService(args.outFile)

    #====================================================================
    # Clean up events

    tray.AddModule(PruneIceTop,
            it_stream=it_stream)

    '''
    tray.AddModule('I3IcePickModule<I3SmallShowerFilter>',
            FilterGeometry = ICgeo,
            TopPulseKey = recoPulses,
            NStationResultName = 'SmallShowerNStation',
            DecisionName = 'SmallShowerFilterPassed')
    '''
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
            config = 'IT81',#args.config,
            output = 'LoudestOnEdge')

    tray.AddModule(LoudestStationOnEdge,
            InputLoudestStation = 'SaturationList',
            config = 'IT81', #args.config,
            output = 'SaturatedOnEdge')

    tray.AddModule(LargestTankCharges,
            ITpulses = recoPulses)

    def IceTop_cuts(frame):
        if frame.Has(track):
            c1 = (np.cos(frame[track].dir.zenith) >= 0.8)   # zenith cut
    #tray.AddModule(IceTop_cuts)

    #====================================================================
    # Run the reconstructions
    tray.AddModule(ShowerLLH.GridLLH, 'ShowerLLHGrid',
            LLHTables = LLHTables,
            binDict = binDict,
            grids = grids,
            recoPulses = recoPulses)
    tray.AddModule(moveSmall)
    if args.config != 'IT81-II':
        def deleteLaputop(frame):
            if frame.Has('LaputopStandard'):
                frame.Delete('LaputopStandard')
                frame.Delete('LaputopStandardParams')
            if frame.Has('LaputopSmallShower'):
                frame.Delete('LaputopSmallShower')
                frame.Delete('LaputopSmallShowerParams')
        tray.AddModule(deleteLaputop, 'deleter')
        tray.AddSegment(Laputop, 'Laputop',
                pulses = recoPulses)
        tray.AddSegment(LaputopSmallShower, 'LaputopSmallShower',
                pulses = recoPulses) 
        tray.AddModule(moveSmall)

    #====================================================================
    def isCoincidentEvent(frame, config):
        icgeo,itgeo = simFunctions.getGeos(config)
        keep       = False
        if frame.Has('LaputopStandard') and frame.Has('ShowerLLH_proton') and frame.Has('Q1') and frame.Has('LoudestStation'):
            reco_dir   = frame['LaputopStandard']
            reco_core  = frame['ShowerLLH_proton']
            lap_params = frame['LaputopStandardParams']
            n_stations = frame['NStations']
            loudest    = frame['LoudestStation']
            Q1         = frame['Q1']
            if reco_dir.fit_status == 0 and reco_core.pos.x != 0 and reco_core.energy > 10**6:
                if lap_params.beta > 1.4 and lap_params.beta < 9.5:
                    if n_stations > 4 and loudest != 1 and Q1 > 6.0: 
                        keep = geo_functions.coincidentCutDirPos(reco_dir, reco_core, itgeo, icgeo)
        print(keep)
        return keep

    tray.AddModule(isCoincidentEvent,'Coincidence', config = 'IT81')

    #====================================================================
    #In Ice Pulse Information
    def pulse_info(frame):
        print('hey')
        geo             = frame['I3Geometry']
        #laputop         = frame['Laputop']
        laputop         = frame['LaputopStandard']
        non_list     = []
        d            = []
        edge_list = [1,2,3,4,5,6,7,13,14,21,22,30,31,40,41,50,51,59,60,67,68,72,73,74,75,76,77,78]
        if args.config == 'IT73':
            non_list += [1,7,14,22,31,79,80] #string numbers not in IC79

        for i in range(78):
            if i+1 not in non_list:
                omgeo =  geo.omgeo[icetray.OMKey(i+1,1)]
                if i+1 not in edge_list:
                    d0  = geo_functions.distanceToInnerString(laputop.pos.x, laputop.pos.y, laputop.dir.zenith, laputop.dir.azimuth, omgeo.position[0], omgeo.position[1], omgeo.position[2])
                    d  += [np.abs(d0)] # closest approach distance to DOM

        allpulses       = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, 'InIcePulses')
        tres         = []
        time         = []
        charges      = []
        dpos         = []
        hlcs         = []
        pulse_x      = []
        pulse_y      = []
        pulse_z      = []
        d_pulse      = []
        slcs         = 0
        off_slcs     = 0
        hlc_noise    = 0
        hlc_signal   = 0

        for omkey, pulses in allpulses:
            omgeo    = geo.omgeo[omkey]
            d_pulse  = i3ps.I3Calculator.closest_approach_distance(laputop, omgeo.position)
            dpos    += [d_pulse for p in pulses]
            tres    += [i3ps.I3Calculator.time_residual(laputop, omgeo.position, p.time) for p in pulses] # time res wrt Cherenkov cone
            time    += [p.time for p in pulses]
            charges += [p.charge for p in pulses]
            pulse_x += [omgeo.position.x for p in pulses]
            pulse_y += [omgeo.position.y for p in pulses]
            pulse_z += [omgeo.position.z for p in pulses]
            hlcs    += [(p.flags & dataclasses.I3RecoPulse.PulseFlags.LC) for p in pulses] # HLC hits
            top_list  = [1,2,3,4,5,6]
            if omkey[1] in top_list:
                for p in pulses:
                    if (p.flags & dataclasses.I3RecoPulse.PulseFlags.LC) == 0:
                        if np.abs(d_pulse) < 200 and p.time > 14800 and p.time < 16500:
                           slcs += 1
                        if np.abs(d_pulse) < 200 and p.time > 0 and p.time < 1700:
                           off_slcs += 1

        frame['hlc_count']      = dataclasses.I3Double(hlcs.count(1))
        frame['slc_count']      = dataclasses.I3Double(slcs)
        frame['off_slc_count']  = dataclasses.I3Double(off_slcs)
        '''
        frame['hits']           = dataclasses.I3VectorInt(hlcs)
        frame['tres']           = dataclasses.I3VectorDouble(tres)
        frame['pulse_x']        = dataclasses.I3VectorDouble(pulse_x)
        frame['pulse_y']        = dataclasses.I3VectorDouble(pulse_y)
        frame['pulse_z']        = dataclasses.I3VectorDouble(pulse_z)
        frame['pulse_time']     = dataclasses.I3VectorDouble(time)
        frame['pulse_charges']  = dataclasses.I3VectorDouble(charges)
        frame['pulse_distance'] = dataclasses.I3VectorDouble(dpos)
        frame['mjd_time']       = dataclasses.I3Double(frame['I3EventHeader'].start_time.mod_julian_day_double)
        frame['d_str']          = dataclasses.I3Double(np.min(d))
        '''

    tray.AddModule(pulse_info, 'PulseInfo')


    #====================================================================
    # Finish
    tray.AddModule(I3TableWriter, tableservice=hdf, keys=keys, SubEventStreams = [it_stream])
  

    tray.Execute()
    tray.Finish()

    print "Time taken: ", time.time() - t0

