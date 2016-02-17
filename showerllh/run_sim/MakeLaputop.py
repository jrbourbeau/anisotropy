#!/usr/bin/env python 

################################################################################
# Runs the grid search on the given files and writes to hdf5 output
################################################################################

from icecube import dataio, toprec
from I3Tray import *
from icecube.tableio import I3TableWriter
from icecube.hdfwriter import I3HDFTableService

import numpy as np
import time, glob, argparse

from i3modules import GetStations, PruneIceTop, moveMCPrimary
from laputop.laputop_standard_traysegment import LaputopStandard
from laputop.laputop_smallshower_traysegment import LaputopSmallShower


if __name__ == "__main__":

    p = argparse.ArgumentParser(
            description='Runs ShowerLLH over a given fileList')
    p.add_argument('-f', '--files', dest='files', nargs='*',
            help='Files to run over')
    p.add_argument('-c', '--config', dest='config',
            help='Detector configuration')
    p.add_argument('-o', '--outFile', dest='outFile',
            help='Output file')
    args = p.parse_args()

    # Starting parameters
    ## NOTE: NOT PREPARED FOR IT59
    recoPulses = 'CleanedHLCTankPulses'
    if args.config == 'IT73':
        it_stream = 'top_hlc_clusters'
    if args.config in ['IT81','IT81-II']:
        it_stream = 'ice_top'

    # Keys to write to frame
    keys  = ['I3EventHeader']
    keys += ['LaputopStandard', 'LaputopStandardParams']
    keys += ['LaputopSmallShower', 'LaputopSmallShowerParams']

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
    # Run the reconstructions

    def lapCheck(frame):
        if frame.Has('LaputopStandard') or frame.Has('LaputopSmallShower'):
            return False

    tray.AddModule(lapCheck)

    tray.AddSegment(LaputopStandard, 'LaputopStandard',
            pulses = recoPulses)

    tray.AddSegment(LaputopSmallShower, 'LaputopSmallShower',
            pulses = recoPulses)

    #tray.AddModule(moveSmall)

    #====================================================================
    # Finish

    tray.AddModule(I3TableWriter, tableservice=hdf, keys=keys, 
            SubEventStreams=[it_stream])

    tray.Execute()
    tray.Finish()

    print "Time taken: ", time.time() - t0





