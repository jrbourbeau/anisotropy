#!/usr/bin/env python 

################################################################################
# Runs the grid search on the given files and writes to hdf5 output
################################################################################

import argparse, time

from icecube import icetray, dataio, toprec
from I3Tray import I3Tray
from icecube.tableio import I3TableWriter
from icecube.hdfwriter import I3HDFTableService

import simFunctions_IT as simFunctions


if __name__ == "__main__":

    p = argparse.ArgumentParser(
            description='Extracts MC primary information from fileList')
    p.add_argument('-c', '--config', dest='config',
            help='Detector configuration')
    p.add_argument('-f', '--files', dest='files', nargs='*',
            help='Files to run over')
    p.add_argument('-o', '--outFile', dest='outFile',
            help='Output file')
    args = p.parse_args()

    # Keys to write to frame
    keys = ['MCPrimary']
    subeventstream = simFunctions.null_stream(args.config)

    t0 = time.time()

    tray = I3Tray()
    tray.AddModule('I3Reader', 'reader', FileNameList=args.files)
    hdf = I3HDFTableService(args.outFile)
    tray.AddModule(I3TableWriter, 'writer', tableservice=hdf, keys=keys, 
            SubEventStreams=[subeventstream])
    tray.Execute()
    tray.Finish()

    print "Time taken: ", time.time() - t0

