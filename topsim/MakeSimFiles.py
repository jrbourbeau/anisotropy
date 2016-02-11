#!/usr/bin/env python 

from icecube import icetray, dataio, toprec
from icecube import dataclasses as dc
from icecube.hdfwriter import I3HDFTableService
from icecube.tableio import I3TableWriter
from I3Tray import *

import sys, time, glob
import numpy as np

from i3modules import GetStations, moveMCPrimary


def maker(config, outFile, fileList):

    # Starting parameters
    recoPulses = 'CleanedHLCTankPulses'
    if config == 'IT73':
        it_stream = 'top_hlc_clusters'
    if config == 'IT81':
        it_stream = 'ice_top'

    # Keys to write to frame
    keys = []
    keys += ['I3EventHeader']
    keys += ['ShowerPlane', 'ShowerPlaneParams']
    keys += ['NStations']
    keys += ['MCPrimary']

    t0 = time.time()

    tray = I3Tray()
    tray.AddModule('I3Reader', FileNameList=fileList)
    hdf = I3HDFTableService(outFile)

    tray.AddModule(moveMCPrimary)

    tray.AddModule(GetStations,
            InputITpulses = recoPulses,
            output = 'NStations')

    tray.AddModule(I3TableWriter, tableservice=hdf, keys=keys,
            SubEventStreams=[it_stream])

    tray.Execute()
    tray.Finish()

    print "Time taken: ", time.time() - t0



if __name__ == "__main__":

    config   = sys.argv[1]
    outFile  = sys.argv[2]
    fileList = sys.argv[3:]
    maker(config, outFile, fileList)







