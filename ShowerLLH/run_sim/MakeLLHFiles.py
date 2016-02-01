#!/usr/bin/env python 

################################################################################
# Runs the grid search on the given files and writes to hdf5 output
################################################################################
from icecube import icetray, dataio, dataclasses, toprec, MuonGun #,ShowerLLH
from icecube import phys_services as i3ps
from icecube import topeventcleaning as tec
from I3Tray import *
from icecube.tableio import I3TableWriter
from icecube.hdfwriter import I3HDFTableService

import numpy as np
import sys, time, glob
import argparse, simFunctions

load('libsmallshower-filter')
sys.path.append('/home/zgriffith/ShowerLLH/segments')
sys.path.append('/home/zgriffith/ShowerLLH/useful')
import geo as geo_functions

from laputop_traysegment import Laputop
from laputop_smallshower_traysegment import LaputopSmallShower

from modules import FindLoudestStation, LoudestStationOnEdge, LargestTankCharges
from modules import GetStations, PruneIceTop, moveMCPrimary, moveSmall, rmPoleMuon
from module import MergeIIIT_simple

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
    pulses_ii   = 'MaskedOfflinePulses'
    recoPulses  = 'CleanedHLCTankPulses'
    #track_ii    =  'PoleMuonLlhFit' # 'LineFit'
    track_it    = 'CleanedHLCTankPulses'#'ShowerPlane'
    resourcedir = '/data/user/zgriffith/ShowerLLH/resources/'

    if args.config == 'IT73':
        it_stream = 'top_hlc_clusters'
        icgeo = 'IC79'
    if args.config in ['IT81','IT81-II']:
        it_stream = 'ice_top'
        icgeo = 'IC86'

    # Load grids for grid scan, llh tables, and bins
    grids = np.load(args.gridFile)
    grids = grids.item()
    LLHTables = np.load(args.llhFile)
    LLHTables = LLHTables.item()
    binDict = LLHTables['bins']
    LLHTables = LLHTables['llhtables']

    # Keys to write to frame
    keys  = ['I3EventHeader']
    keys += ['ShowerPlane', 'ShowerPlaneParams']
    keys += ['LoudestStation', 'LoudestOnEdge']
    keys += ['SaturationList', 'SaturatedOnEdge']
    keys += ['Q1', 'Q2', 'Q3', 'Q4']
    keys += [recoPulses, pulses_ii, 'NStations']
    for comp in LLHTables.keys():
        keys += ['ShowerLLH_'+comp, 'maxLLH_'+comp]
    keys += ['Laputop','LaputopParams']
    keys += ['MCPrimary', 'ShowerLLHParams']
    keys += ['hlc_count', 'slc_count', 'd_str', 'mjd_time']

    #keys recorded only if you want individual hit information
    #keys += ['pulse_x', 'pulse_y', 'pulse_z', 'pulse_distance', 'all_slc']
    #keys += ['hits', 'tres', 'pulse_time', 'pulse_charges', 'd_pulse']

    #keys += ['I3MCTree','muon_energy','I3TriggerHierarchy']

    t0 = time.time()

    tray = I3Tray()
    tray.AddModule('I3Reader', FileNameList=args.files)
    hdf = I3HDFTableService(args.outFile)

    #====================================================================
    # Clean up events

    tray.AddModule(PruneIceTop,
            it_stream=it_stream)

    tray.AddModule(moveMCPrimary)
    #tray.AddModule(rmPoleMuon)


    #====================================================================
    # Merge events
    tray.AddModule(MergeIIIT_simple, 'merged_ii_it',
        IceTopReco   = track_it,          # Previously ShowerCombined
        #InIceReco    = track_ii,          # Previously LineFit
        IceTopPulses = recoPulses,
        InIcePulses  = pulses_ii)

    tray.AddModule(lambda frame:frame['I3EventHeader'].sub_event_stream == 'merged_ii_it' or frame['I3EventHeader'].sub_event_stream == 'icetop_only')
    #tray.AddModule(lambda frame:frame['I3EventHeader'].sub_event_stream == 'ice_top')
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
            config = args.config,
            output = 'LoudestOnEdge')

    tray.AddModule(LoudestStationOnEdge,
            InputLoudestStation = 'SaturationList',
            config = args.config,
            output = 'SaturatedOnEdge')

    tray.AddModule(LargestTankCharges,
            ITpulses = recoPulses)
    #====================================================================
    # Run the reconstructions
    '''
    tray.AddModule(ShowerLLH.GridLLH, 'ShowerLLHGrid',
            LLHTables = LLHTables,
            binDict = binDict,
            grids = grids,
            recoPulses = recoPulses)
    '''

    if True:#args.config != 'IT81-II':
        def deleteLaputop(frame):
            if frame.Has('LaputopStandard'):
                frame.Delete('LaputopStandard')
                frame.Delete('LaputopStandardParams')
            if frame.Has('LaputopSmallShower'):
                frame.Delete('LaputopSmallShower')
                frame.Delete('LaputopSmallShowerParams')
            if frame.Has('Laputop'):
                frame.Delete('Laputop')
                frame.Delete('LaputopParams')
            if frame.Has('LaputopSmall'):
                frame.Delete('LaputopSmall')
                frame.Delete('LaputopSmallParams')
        tray.AddModule(deleteLaputop, 'deleter')
        tray.AddSegment(Laputop, 'Laputop',
                pulses = recoPulses,
                snowfactor = 2.25)
        tray.AddSegment(LaputopSmallShower, 'LaputopSmall',
                pulses = recoPulses,
                snowfactor = 2.25)
        tray.AddModule(moveSmall)

    tray.AddModule(moveSmall)


    tray.AddModule("I3ParaboloidFitter","paraboloid")(
        ("SeedService","LaputopToprecSeed"),
        ("LogLikelihood","LaputopToprecLike2"),
        ("VertexStepSize",5.0*I3Units.m),
        ("MaxMissingGridPoints",1),
        ("GridpointVertexCorrection","seedprep"),
        ("Minimizer","Laputopminuit"),
    )  

    #====================================================================
    #In Ice Pulse Information
    def pulse_info(frame):
        geo             = frame['I3Geometry']
        laputop         = frame['Laputop']
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

        if frame['I3EventHeader'].sub_event_stream == 'merged_ii_it':
            allpulses       = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, 'MaskedOfflinePulses')
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
            hlc_noise    = 0
            hlc_signal   = 0

            for omkey, pulses in allpulses:
                omgeo    = geo.omgeo[omkey]
                d_pulse  = i3ps.I3Calculator.closest_approach_distance(laputop, omgeo.position)
                dpos    += [d_pulse for p in pulses]
                tres    += [i3ps.I3Calculator.time_residual(laputop, omgeo.position, p.time) for p in pulses] # time res wrt Cherenkov cone
                #time    += [p.time for p in pulses]
                #charges += [p.charge for p in pulses]
                #pulse_x += [omgeo.position.x for p in pulses]
                #pulse_y += [omgeo.position.y for p in pulses]
                #pulse_z += [omgeo.position.z for p in pulses]
                hlcs    += [(p.flags & dataclasses.I3RecoPulse.PulseFlags.LC) for p in pulses] # HLC hits

                if (p.flags & dataclasses.I3RecoPulse.PulseFlags.LC) != 0:
                    if np.abs(d_pulse) > 200 or p.time < 14800 or p.time > 17500:
                        hlc_noise  += 1
                    if np.abs(d_pulse) < 200 and p.time > 14800 and p.time < 17500:
                        hlc_signal += 1
    
                top_list  = [1,2,3,4,5,6]
                if omkey[1] in top_list:
                    for p in pulses:
                        if (p.flags & dataclasses.I3RecoPulse.PulseFlags.LC) == 0:
                            if np.abs(d_pulse) < 200 and p.time > 14800 and p.time < 16500:
                               slcs += 1
            '''            
            frame['hits']           = dataclasses.I3VectorInt(hlcs)
            frame['tres']           = dataclasses.I3VectorDouble(tres)
            frame['pulse_x']        = dataclasses.I3VectorDouble(pulse_x)
            frame['pulse_y']        = dataclasses.I3VectorDouble(pulse_y)
            frame['pulse_z']        = dataclasses.I3VectorDouble(pulse_z)
            frame['pulse_time']     = dataclasses.I3VectorDouble(time)
            frame['pulse_charges']  = dataclasses.I3VectorDouble(charges)
            frame['pulse_distance'] = dataclasses.I3VectorDouble(dpos)
            '''

            frame['hlc_count']      = dataclasses.I3Double(hlcs.count(1))
            frame['all_slc']      = dataclasses.I3Double(hlcs.count(0))
            frame['slc_count']      = dataclasses.I3Double(slcs)
            frame['hlc_noise']      = dataclasses.I3Double(hlc_noise)
            frame['hlc_signal']     = dataclasses.I3Double(hlc_signal)
            frame['mjd_time']       = dataclasses.I3Double(frame['I3EventHeader'].start_time.mod_julian_day_double)
            frame['d_str']          = dataclasses.I3Double(np.min(d))
        else:
            frame['mjd_time']       = dataclasses.I3Double(frame['I3EventHeader'].start_time.mod_julian_day_double)
            frame['hlc_count']      = dataclasses.I3Double(np.nan)
            frame['slc_count']      = dataclasses.I3Double(np.nan)
            frame['hlc_noise']      = dataclasses.I3Double(np.nan)
            frame['hlc_signal']     = dataclasses.I3Double(np.nan)
            frame['d_str']          = dataclasses.I3Double(np.min(d))
            frame['all_slc']        = dataclasses.I3Double(np.nan)
            

    #tray.AddModule(pulse_info, 'PulseInfo')

    def icecube_muon_energies(frame):
        energies = []
        surface = MuonGun.Cylinder(1000,500)
        for p in MuonGun.muons_at_surface(frame, surface):
            energies.append(p.energy)

        frame['muon_energy'] =  dataclasses.I3VectorDouble(energies)

    #tray.AddModule(icecube_muon_energies, 'muonEnergies')

    def isCoincidentEvent(frame, config):
        icgeo,itgeo = simFunctions.getGeos(config)
        if frame.Has('Laputop') and frame.Has('ShowerLLH_proton') and frame.Has('d_str'):
            reco_dir = frame['Laputop']
            reco_core = frame['ShowerLLH_proton']
            d_str = frame['d_str']
            if not np.isnan(reco_dir.pos.x) and reco_core.pos.x != 0 and d_str.value < 60:
                return geo_functions.coincidentCutDirPos(reco_dir, reco_core, itgeo, icgeo)
            else:
                return False
        else:
            return False

    #tray.AddModule(isCoincidentEvent,'Coincidence', config = 'IT81')

    def candidates(frame):
        I3EventHeader = frame['I3EventHeader']
        print(I3EventHeader.run)
    #tray.AddModule(candidates,'candidates')
    #====================================================================
    # Finish
    tray.AddModule(I3TableWriter, tableservice=hdf, keys=keys, SubEventStreams = ['merged_ii_it', 'icetop_only'])
    '''
    tray.AddModule("I3Writer", 'i3writer',
        filename=args.outFile,
        DropOrphanStreams=[icetray.I3Frame.DAQ]
    ) 
    '''

    tray.Execute()
    tray.Finish()

    print "Time taken: ", time.time() - t0

