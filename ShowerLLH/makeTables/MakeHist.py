#!/usr/bin/env python 

#############################################################################
# Makes the probability tables using AddToHist and records as .npy files
#############################################################################

from icecube import icetray, dataio
from I3Tray import I3Tray
import time, sys, AddToHist
import numpy as np

def maker(outFile, fileList):

    # Starting parameters
    resourcedir = '/net/user/zgriffith/ShowerLLH/resources/'
    recoPulses = 'CleanedHLCTankPulses'

    # Import bins
    binDict = np.load(resourcedir + 'ShowerLLH_bins.npy')
    binDict = binDict.item()

    t0 = time.time()
    tray = I3Tray()
    tray.AddModule('I3Reader', 'reader', FileNameList=fileList)
    tray.AddModule(AddToHist.fillHist, 'fillHist', binDict=binDict,
            recoPulses=recoPulses, outFile=outFile)
    tray.Execute()
    tray.Finish()

    print "Time taken: ", time.time() - t0



if __name__ == "__main__":

    outFile  = sys.argv[1]
    fileList = sys.argv[2:]
    maker(outFile, fileList)















