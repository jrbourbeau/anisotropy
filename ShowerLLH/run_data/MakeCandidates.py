#!/usr/bin/env python 

################################################################################
# Runs the grid search on the given files and writes to hdf5 output
################################################################################
from icecube import icetray, dataio, dataclasses
from I3Tray import *
from icecube.tableio import I3TableWriter
from icecube.hdfwriter import I3HDFTableService

import numpy as np
import sys, time, glob, simFunctions
import argparse

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
    args = p.parse_args()

    # Starting parameters
    recoPulses  = 'CleanedHLCTankPulses'
    track_it    = 'ShowerPlane'
    resourcedir = '/data/user/zgriffith/ShowerLLH/resources/'

    if args.config == 'IT73':
        it_stream = 'top_hlc_clusters'
        ICgeo = 'IC79'
    if args.config in ['IT81','IT81-II']:
        it_stream = 'ice_top'
        ITgeo = 'IT81'
        ICgeo = 'IC86'
        filter_mask = 'QFilterMask'

    t0 = time.time()

    tray = I3Tray()
    tray.AddModule('I3Reader', FileNameList=args.files)

    #====================================================================
    def isCandidate(frame):
        time = np.load('/data/user/zgriffith/datasets/candidate_times.npy')
        if frame['I3EventHeader'].start_time.mod_julian_day_double in time:
            return True
        else:
            return False  

    tray.AddModule(isCandidate,'Candidate')

    #====================================================================
    # Finish
    tray.AddModule("I3Writer", 'i3writer',
        filename=args.outFile,
        DropOrphanStreams=[icetray.I3Frame.DAQ]
    )  

    tray.Execute()
    tray.Finish()

    print "Time taken: ", time.time() - t0

