#!/usr/bin/env python

from icecube import icetray, dataio, dataclasses
from I3Tray import *

import numpy as np
import time

load("libDOMcalibrator")
load("libDomTools") 
load("libFeatureExtractor") 
load("libicepick")
load("libdaq-decode")
load("libBadDomList")

class getDST(icetray.I3ConditionalModule):

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter('outFile', 'Output file name', '')
        self.AddParameter('recoPulses', 'Pulses used in reconstruction', '')
        self.AddParameter('dstReco', 'Directional dst reconstruction', '')
        self.AddParameter('nFiles', 'Number of files run over', '')

    def Configure(self):
        self.outFile = self.GetParameter('outFile')
        self.recoPulses = self.GetParameter('recoPulses')
        self.dstReco = self.GetParameter('dstReco')
        self.nFiles = self.GetParameter('nFiles')
        self.FilterMinBias = 'FilterMinBias_09'

        # Create dictionary for storage
        self.d = {}
        keys = ['nchannel','dstZenith','dstAzimuth','zen','azi','energy','type']
        for key in keys:
            self.d[key] = []

    def Physics(self, frame):

        # Make sure we have the information needed
        checks = ['I3MCTree', 'FilterMask', self.recoPulses, self.dstReco]
        if not all([frame.Has(check) for check in checks]):
            self.PushFrame(frame)
            return

        # Make sure an actual primary triggered the reconstruction
        c0 = (len(frame['I3MCTree'].primaries) != 0)
        c1 = (frame['FilterMask'][self.FilterMinBias].condition_passed == 1)
        if not c0 or not c1:
            self.PushFrame(frame)
            return

        # Extract MCPrimary information - take most energetic if multiple
        Emax = 0
        for primary in frame['I3MCTree'].primaries:
            if primary.energy > Emax:
                Emax = primary.energy
                MCPrimary = primary

        # Calculate nchannel
        dstPulses = frame[self.recoPulses]
        dstPulses = dstPulses.apply(frame)
        nchannel = len(dstPulses.keys())

        # Record desired information
        self.d['nchannel'].append(nchannel)
        self.d['dstZenith'].append(frame[self.dstReco].dir.zenith)
        self.d['dstAzimuth'].append(frame[self.dstReco].dir.azimuth)
        self.d['zen'].append(MCPrimary.dir.zenith)
        self.d['azi'].append(MCPrimary.dir.azimuth)
        self.d['energy'].append(MCPrimary.energy)
        #self.d['type'].append(MCPrimary.pdgEncoding)
        self.d['type'].append(MCPrimary.type)

        self.PushFrame(frame)
        return

    def Finish(self):
        for key in self.d.keys():
            self.d[key] = np.asarray(self.d[key])
        self.d['nFiles'] = self.nFiles
        np.save(self.outFile, self.d)
        return



def maker(outFile, fileList):

    nFiles = len([f for f in fileList if 'GCD' not in f])

    t0 = time.time()
    tray = I3Tray()
    tray.AddModule('I3Reader', FileNameList=fileList)

    ##========================================================================
    ## Get the same pulse series used in DST

    tray.AddModule("QConverter", "convert")

    tray.AddModule( "I3DOMLaunchCleaning", "LaunchCleaning_09") (
        ( "InIceInput", "InIceRawData" ), # Default
        ( "IceTopInput", "IceTopRawData" ), # Default
        ( "InIceOutput", "CleanInIceRawData" ), # Default
        ( "IceTopOutput", "CleanIceTopRawData" ), # Default
        ( "FirstLaunchCleaning", True ), # Default
        ( "CleanedKeysList", "" ), # Default
        ( "CleanedKeys", [OMKey( 38, 59 ),  # Blackberry
                          OMKey( 6, 11 )   # Discworld - Meteor DOM
                          ])
        )

    tray.AddModule( "I3LCCleaning", "InIceLCCleaning_09" ) (
        ( "InIceInput", "CleanInIceRawData" ),
        ( "InIceOutput", "HLCInIceRawData")
        )

    tray.AddModule( "I3DOMcalibrator", "CalibrateInIce_09" ) (
        ( "InputRawDataName", "HLCInIceRawData" ), # ! Use cleaned raw data
        ( "OutputATWDDataName", "CalibratedATWD" ), # Default
        ( "OutputFADCDataName", "CalibratedFADC" ), # Default
        ( "OutputToFile", False ), # Default
        ( "OutputFileName", "DOMcalibrator.root" ), # Default
        ( "CalibrationMode", 0 ), # Default
        ( "SubtractBaseline", False ), # Default
        ( "SubtractTransitTime", True ), # Default (new default)
        ( "CalibrateDataWithSLC", False ), # Default
        ( "KeepCStampRedundantInfo", False ), # Default
        ( "CorrectPedestalDroop", False ), # Default
        ( "CorrectPedestalDroopDualTau", True ), # Correct for electronics droop
        ( "ATWDSaturationLevel", 1022 ), # Default
        ( "FADCSaturationLevel", 1022 ), # Default
        ( "FADCTimeOffset", 0 ), # keep timing wrong to simulate pole correctly
        )

    # Extract pulses for *muon filter*
    tray.AddModule( "I3FeatureExtractor", "Features_09") (
        ( "RawReadoutName", "HLCInIceRawData" ), # ! Use cleand raw data
        ( "CalibratedFADCWaveforms", "CalibratedFADC" ), # Default
        ( "CalibratedATWDWaveforms", "CalibratedATWD" ), # Default
        ( "InitialPulseSeriesReco", "MuonPulseSeriesReco" ), 
        ( "InitialHitSeriesReco", "MuonHitSeriesReco" ), # Default
        ( "DisableHitSeries", True ), # Default
        ( "MaxNumHits", 0 ), # Default
        ( "FastPeakUnfolding", -1 ), # Don't extract multiple pulses
        ( "FastFirstPeak", 7 ), # ! Baseline from first 3 bins
        ( "MinSpeWidth", 4 ), # Default
        ( "ExclusionSize", 5 ), # Ensure backwards compatibility to online
        ( "MaxSpeWidth", 20 ), # Default
        ( "ADCThreshold", 1.1 ), # Threshold at 1.1 times the hardware setting
        ( "PMTTransit", 2 ),  # Ensure backwards compatibility to online
        )

    tray.AddModule( "I3TimeWindowCleaning<I3RecoPulse>", "TimeWindow_09" ) (
        ( "InputResponse", "MuonPulseSeriesReco" ), # Extracted pulse series
        ( "OutputResponse", "TWCMuonPulseSeriesReco" ), # Cleaned pulse series
        ( "TimeWindow", 6000 * I3Units.ns ), # ! 6 usec time window
        )

    ##========================================================================
    ## Extract DST information

    tray.AddModule(getDST, 'getDST',
        outFile=outFile,
        recoPulses='TWCMuonPulseSeriesReco',
        dstReco='PoleMuonLlhFit',
        nFiles=nFiles)

    tray.Execute()
    tray.Finish()


if __name__ == "__main__":

    config   = sys.argv[1]
    outFile  = sys.argv[2]
    fileList = sys.argv[3:]
    maker(outFile, fileList)

