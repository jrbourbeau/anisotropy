#!/usr/bin/env python 

################################################################################
# Runs the grid search on the given files and writes to hdf5 output
################################################################################

from icecube import dataio, toprec, ShowerLLH
from I3Tray import *
from icecube.tableio import I3TableWriter
from icecube.hdfwriter import I3HDFTableService

import numpy as np
import time, glob, argparse

import dataFunctions as df
from i3modules import GetStations, PruneIceTop
from i3modules import FindLoudestStation, LoudestStationOnEdge
from i3modules import LargestTankCharges, moveSmall


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
    p.add_argument('--old', dest='old',
            default=False, action='store_true',
            help='Option to run old (slower) python version')
    p.add_argument('-o', '--outFile', dest='outFile',
            help='Output file')
    args = p.parse_args()

    # Starting parameters
    recoPulses  = df.recoPulses(args.config)
    it_stream   = df.it_stream(args.config)
    filter_mask = df.filter_mask(args.config)
    it_geo = df.it_geo(args.config)
    recoTrack = df.recoTrack(args.config)

    # Load grids for grid scan, llh tables, and bins
    grids = np.load(args.gridFile)
    grids = grids.item()
    LLHTables = np.load(args.llhFile)
    LLHTables = LLHTables.item()
    binDict = LLHTables['bins']
    LLHTables = LLHTables['llhtables']

    # Keys to write to frame
    keys  = []
    keys += ['I3EventHeader', filter_mask]
    keys += ['ShowerPlane', 'ShowerPlaneParams']
    keys += ['LoudestStation', 'LoudestOnEdge']
    keys += ['SaturationList', 'SaturatedOnEdge']
    keys += ['Q1', 'Q2', 'Q3', 'Q4']
    keys += ['NStations']
    for comp in LLHTables.keys():
        keys += ['ShowerLLH_'+comp, 'ShowerLLHParams_'+comp]

    t0 = time.time()

    tray = I3Tray()
    tray.AddModule('I3Reader', FileNameList=args.files)
    hdf = I3HDFTableService(args.outFile)

    #====================================================================
    # Clean up events

    tray.AddModule(PruneIceTop, 
            it_stream = it_stream)

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
            config = it_geo,
            output = 'LoudestOnEdge')

    tray.AddModule(LoudestStationOnEdge,
            InputLoudestStation = 'SaturationList',
            config = it_geo,
            output = 'SaturatedOnEdge')

    tray.AddModule(LargestTankCharges,
            ITpulses = recoPulses)

    if args.config in ['IT81-III','IT81-IV']:
        tray.AddModule(moveSmall)

    def IceTop_cuts(frame):
        passCuts = False
        if frame.Has(recoTrack):
            c1 = (np.cos(frame[recoTrack].dir.zenith) >= 0.8)   # zenith cut
            c2 = (not frame['LoudestOnEdge'].value)             # edge cut
            c3 = (frame['Q1'].value >= 6)                       # charge cut
            passCuts = c1*c2*c3
        return passCuts

    #tray.AddModule(IceTop_cuts)


    #====================================================================
    # Run the reconstruction

    # C++ implementation for speed
    if not args.old:
        gridx, gridy = np.transpose(grids[0])
        for comp in LLHTables.keys():
            llh = list(LLHTables[comp].flatten())
            tray.AddModule('ShowerLLH', 'ShowerLLH_'+comp,
                    recoPulses=recoPulses,
                    comp=comp,
                    ebins=list(binDict[0][1]),
                    sbins=list(binDict[1][1]),
                    dbins=list(binDict[2][1]),
                    cbins=list(binDict[3][1]),
                    gridx=gridx, gridy=gridy,
                    directionFit=recoTrack,
                    llhTable=llh)

    # Python implementation more flexible with binning
    else:
        tray.AddModule(ShowerLLH.GridLLH, 'ShowerLLHGrid',
                LLHTables = LLHTables,
                binDict = binDict,
                grids = grids,
                recoPulses = recoPulses,
                directionFit = recoTrack)

    #====================================================================
    # Finish

    tray.AddModule(I3TableWriter, tableservice=hdf, keys=keys, 
            SubEventStreams=[it_stream])

    tray.Execute()
    tray.Finish()

    print "Time taken: ", time.time() - t0




