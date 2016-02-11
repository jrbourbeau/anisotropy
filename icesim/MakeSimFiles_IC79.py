#!/usr/bin/env python

from icecube import icetray, dataio, dataclasses
from I3Tray import *

import numpy as np
import time

load("libDOMcalibrator")
load("libDomTools") 
load("libFeatureExtractor") 


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
        self.FilterMinBias = 'FilterMinBias_10'

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

    tray.AddModule("I3DOMLaunchCleaning","4_BadDomCleaning")( 
        ("InIceInput","InIceRawData"), 
        ("IceTopInput","IceTopRawData"), 
        ("InIceOutput","CleanInIceRawDataL1"), 
        ("IceTopOutput","CleanIceTopRawDataL1"), 
        ("FirstLaunchCleaning",False),
        ("CleanedKeysList",""),
        ("IcePickServiceKey",""),
        ("CleanedKeys",[OMKey(38,59),# Blackberry 
                        OMKey(6,11), # Discworld - Meteor DOM
                        OMKey(68,42)  # Krabba
                        ]) 
        )

    tray.AddModule("I3LCCleaning","6_InIceLCClean") (
        ("IcePickServiceKey", ""),
        ("InIceInput", "CleanInIceRawDataL1"), 
        ("InIceOutput", "HLCInIceRawDataL1") 
        ) 

    tray.AddModule("I3DOMcalibrator","8_StdDomcal")( 
        ("InputRawDataName","HLCInIceRawDataL1"), 
        ("OutputFADCDataName","CalibratedFADC"),
        ("OutputATWDDataName","CalibratedATWD"),
        ("ATWDSaturationLevel",900),
        ("FADCSaturationLevel",900),
        ("CalibrateDataWithSLC",False),
        ("CalibrationMode",0),
        ("CorrectPedestalDroopDualTau",True), 
        ("CorrectPedestalDroop",False),
        ("IcePickServiceKey",""),
        ("KeepCstampRedundantInfo",False),
        ("OutputToFile",False),
        ("SubtractBaseline",False),
        ("FADCTimeOffset",0),
        ("SubtractTransitTime",True),
        ) 

    tray.AddModule( "I3FeatureExtractor", "9_Features" ) ( 
        ( "RawReadoutName", "HLCInIceRawDataL1" ), # ! Use clean data 
        ( "CalibratedATWDWaveforms", "CalibratedATWD" ), # Default 
        ( "CalibratedFADCWaveforms", "CalibratedFADC" ), # Default 
        ( "InitialPulseSeriesReco", "MuonPulseSeriesReco" ), # ! Different Name 
        ( "InitialHitSeriesReco", "InitialHitSeriesReco"), # Just for i3moni.
        ( "DisableHitSeries", False ), # ! On, just for I3moni
        ( "MaxNumHits", 0 ), # Default 
        ( "FastPeakUnfolding", -1 ), # ! Multi pulse extraction 
        ( "FastFirstPeak", 7 ), # ! FADC and baseline from first 3 bin
        ( "ADCThreshold", 1.8 ), # ! We pulse above 0.25 pe
        ( "ExclusionSize", 1 ),
        ( "MaxSPEWidth", 20 ), # Default 
        ( "MinSPEWidth", 4 ), # Default 
        ( "PMTTransit", -1 ), # Default 
        ( "UseNewDiscThreshold", False), # ! Use the old thereshold calculation
        ( "TinyThreshold", 0.00 ), # to be confirmed by data/MC comparison
        ( "IcePickServiceKey","")
        ) 

    tray.AddModule("I3TimeWindowCleaning<I3RecoPulse>", "10_TimeWindowClean") ( 
        ( "InputResponse", "MuonPulseSeriesReco" ), 
        ( "OutputResponse", "TWCMuonPulseSeriesReco" ), 
        ( "TimeWindow", 6000*I3Units.ns ),
        ( "IcePickServiceKey","")
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

