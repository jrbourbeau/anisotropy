#!/usr/bin/env python 

################################################################################
# Runs the grid search on the given files and writes to hdf5 output
################################################################################

from icecube import icetray, dataio, toprec
from I3Tray import I3Tray
from icecube.tableio import I3TableWriter
from icecube.hdfwriter import I3HDFTableService
import sys, time

def maker(outFile, fileList):

    # Keys to write to frame
    keys = ['MCPrimary']
    outFile = outFile.replace('.', '_MC.')

    t0 = time.time()

    tray = I3Tray()
    tray.AddModule('I3Reader', 'reader', FileNameList=fileList)
    hdf = I3HDFTableService(outFile)
    tray.AddModule(I3TableWriter, 'writer', tableservice=hdf, keys=keys, 
            SubEventStreams=['nullsplit'])
    tray.AddModule('TrashCan', 'yeswecan')
    tray.Execute()
    tray.Finish()

    print "Time taken: ", time.time() - t0


if __name__ == "__main__":

    outFile  = sys.argv[1]
    fileList = sys.argv[2:]
    maker(outFile, fileList)
