#!/usr/bin/env python 

#############################################################################
# Makes the probability tables using AddToHist and records as .npy files
#############################################################################

from icecube import icetray, dataio
from I3Tray import I3Tray
import myGlobals as my
import numpy as np
import argparse, time, AddToHist


if __name__ == "__main__":

    # Global variables setup for path names
    my.setupShowerLLH(verbose=False)

    p = argparse.ArgumentParser(
            description='Builds binned histograms for use with ShowerLLH')
    p.add_argument('-f', '--files', dest='files', nargs='*',
            help='Input filelist to run over')
    p.add_argument('-b', '--bintype', dest='bintype',
            default='standard',
            choices=['standard','nozenith','logdist'],
            help='Option for a variety of preset bin values')
    p.add_argument('-o', '--outFile', dest='outFile',
            help='Output filename')
    args = p.parse_args()

    # Starting parameters
    recoPulses = 'CleanedHLCTankPulses'
    binFile = '%s/ShowerLLH_bins.npy' % my.llh_resource

    # Import bins
    binDict = np.load(binFile)
    binDict = binDict.item()
    binDict = binDict[args.bintype]

    # Execute
    t0 = time.time()
    tray = I3Tray()
    tray.AddModule('I3Reader', 'reader', FileNameList=args.files)
    tray.AddModule(AddToHist.fillHist, 'fillHist', binDict=binDict,
            recoPulses=recoPulses, outFile=args.outFile)
    tray.Execute()
    tray.Finish()

    print "Time taken: ", time.time() - t0





