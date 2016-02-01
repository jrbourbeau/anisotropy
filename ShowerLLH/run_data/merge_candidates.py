#!/usr/bin/env python 

################################################################################
# Runs the grid search on the given files and writes to hdf5 output
################################################################################
from icecube import icetray, dataio
from I3Tray import *
from icecube.tableio import I3TableWriter
from icecube.hdfwriter import I3HDFTableService

import numpy as np
import sys, time, glob, simFunctions

if __name__ == "__main__":

    # Starting parameters
    recoPulses  = 'CleanedHLCTankPulses'
    track_it    = 'ShowerPlane'
    resourcedir = '/data/user/zgriffith/ShowerLLH/resources/'

    prefix = '/data/user/zgriffith/ShowerLLH/IT81-II_data/files/burn_sample/'
    files = glob.glob(prefix+'*.i3')
    t0 = time.time()

    tray = I3Tray()
    tray.AddModule('I3Reader', FileNameList=files)

    #====================================================================
    # Finish
    tray.AddModule("I3Writer", 'i3writer',
        filename=prefix+'Candidates.i3',
        DropOrphanStreams=[icetray.I3Frame.DAQ]
    )  

    tray.Execute()
    tray.Finish()

    print "Time taken: ", time.time() - t0

