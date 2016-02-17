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

import simFunctions_IT as simFunctions
from i3modules import PruneIceTop, moveMCPrimary


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
    p.add_argument('--test', dest='test',
            default=False, action='store_true',
            help='Option for running test sample off cluster')
    p.add_argument('--cpp', dest='cpp',
            default=False, action='store_true',
            help='Run in (static) c++ mode')
    args = p.parse_args()

    # Starting parameters
    recoPulses = simFunctions.recoPulses(args.config)
    it_stream  = simFunctions.it_stream(args.config)

    # Load grids for grid scan, llh tables, and bins
    grids = np.load(args.gridFile)
    grids = grids.item()
    LLHTables = np.load(args.llhFile)
    LLHTables = LLHTables.item()
    binDict = LLHTables['bins']
    LLHTables = LLHTables['llhtables']

    # Keys to write to frame
    keys  = []
    keys += ['I3EventHeader']
    for comp in LLHTables.keys():
        keys += ['ShowerLLH_'+comp, 'ShowerLLHParams_'+comp]
    keys += ['MCPrimary']

    t0 = time.time()

    tray = I3Tray()
    tray.AddModule('I3Reader', FileNameList=args.files)
    hdf = I3HDFTableService(args.outFile)

    #====================================================================
    # Clean up events

    tray.AddModule(PruneIceTop,
            it_stream=it_stream)

    tray.AddModule(moveMCPrimary)

    #====================================================================
    # Run the reconstruction

    # C++ mode faster
    if args.cpp:
        gridx, gridy = np.transpose(grids[0])
        for comp in LLHTables.keys():
            llh = list(LLHTables[comp].flatten())
            tray.AddModule('ShowerLLH', 'ShowerLLH_'+comp,
                    recoPulses=recoPulses,
                    comp=comp,
                    ebins=list(binDict[0][1]),
                    #zbins=list(binDict[1][1]),
                    sbins=list(binDict[1][1]),
                    dbins=list(binDict[2][1]),
                    cbins=list(binDict[3][1]),
                    gridx=gridx, gridy=gridy,
                    llhTable=llh)

    # Python implementation fast enough to run on simulation
    # (and more flexible for binning tests)
    else:
        tray.AddModule(ShowerLLH.GridLLH, 'ShowerLLHGrid',
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





