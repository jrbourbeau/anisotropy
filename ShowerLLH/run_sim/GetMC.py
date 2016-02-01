#!/usr/bin/env python 

################################################################################
# Runs the grid search on the given files and writes to hdf5 output
################################################################################

from icecube import icetray, tableio, dataclasses, dataio, toprec
from I3Tray import I3Tray
from icecube.tableio import I3TableWriter
from icecube.hdfwriter import I3HDFTableService
import sys, time
sys.path.append('/home/santander/ShowerLLH/useful')
from LLHFunctions import simList

if __name__ == "__main__":

    # Input parameters
    config = sys.argv[1]         # IT59 or IT73
    sim = sys.argv[2]
    group = sys.argv[3]

    files, name = simList(config, sim, group=group)

    # Keys to write to frame
    keys = ['MCPrimary']
    name = name.replace('.', '_MC.')

    t0 = time.time()

    tray = I3Tray()
    tray.AddModule('I3Reader', 'reader', FileNameList=files)
    hdf = I3HDFTableService(name)
    tray.AddModule(I3TableWriter, 'writer', tableservice=hdf, keys=keys, 
            SubEventStreams=['nullsplit'])
    tray.AddModule('TrashCan', 'yeswecan')
    tray.Execute()
    tray.Finish()

    print "Time taken: ", time.time() - t0

