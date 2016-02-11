#!/usr/bin/env python 

from icecube import icetray, dataio, toprec, dst
from I3Tray import *

import sys, time, glob
import numpy as np


class getDST(icetray.I3ConditionalModule):

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter('outFile', 'Output file name', '')
        self.AddParameter('config', 'Detector configuration', '')

    def Configure(self):
        self.outFile = self.GetParameter('outFile')
        self.config  = self.GetParameter('config')
        self.d = {}
        for key in ['nchannel','reco2','zen','azi','energy','type']:
            self.d[key] = []

        self.FilterMinBias = 'FilterMinBias_'
        if config == 'IC86':
            self.dst = 'I3DST11'
            self.FilterMask = 'FilterMask'
            self.FilterMinBias += '11'
        if config == 'IC86-II':
            self.dst = 'I3DST12Reco'
            self.FilterMask = 'QFilterMask'
            self.FilterMinBias += '12'
        if config == 'IC86-III':
            self.dst = 'I3DST13Reco'
            self.FilterMask = 'QFilterMask'
            self.FilterMinBias += '13'


    def Physics(self, frame):

        # Make sure we have the information needed
        checks = ['I3MCTree', self.dst, self.FilterMask]
        if not all([frame.Has(check) for check in checks]):
            self.PushFrame(frame)
            return

        # Make sure an actual primary triggered the reconstruction
        c0 = (len(frame['I3MCTree'].primaries) != 0)
        # Passed min bias filter
        c1 = (frame[self.FilterMask][self.FilterMinBias].condition_passed == 1)
        if not c0 or not c1:
            self.PushFrame(frame)
            return

        # Extract MCPrimary information - take most energetic if multiple
        Emax = 0
        for primary in frame['I3MCTree'].primaries:
            if primary.energy > Emax:
                Emax = primary.energy
                MCPrimary = primary

        # Record desired information
        self.d['nchannel'].append(frame[self.dst].ndom)
        self.d['reco2'].append(frame[self.dst].reco_2)
        self.d['zen'].append(MCPrimary.dir.zenith)
        self.d['azi'].append(MCPrimary.dir.azimuth)
        self.d['energy'].append(MCPrimary.energy)
        self.d['type'].append(MCPrimary.type)

        self.PushFrame(frame)
        return

    def Finish(self):
        for key in self.d.keys():
            self.d[key] = np.asarray(self.d[key])
        np.save(self.outFile, self.d)
        return
        

def maker(config, outFile, fileList):

    t0 = time.time()
    tray = I3Tray()
    tray.AddModule('I3Reader', FileNameList=fileList)
    tray.AddModule(getDST, outFile=outFile, config=config)
    tray.Execute()
    tray.Finish()

    print "Time taken: ", time.time() - t0



if __name__ == "__main__":

    config   = sys.argv[1]
    outFile  = sys.argv[2]
    fileList = sys.argv[3:]
    maker(config, outFile, fileList)







