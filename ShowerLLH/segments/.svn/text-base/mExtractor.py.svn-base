#!/usr/bin/env python

#############################################################################
# Records the energy, zenith, distances, snow depths, and charges for the
# creation of the LLH histogram
#############################################################################

from icecube import icetray, tableio, dataio, toprec, simclasses
from icecube import dataclasses as dc
from numpy import *
import time

class InIceLLH_millipede(icetray.I3Module):

    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        self.AddOutBox('OutBox')
        self.AddParameter('millipede', 'Name of millipede particle vect', '')

    def Configure(self):
        self.millipede = self.GetParameter('millipede')
        pass

    def Geometry(self, frame):
        self.geometry = frame['I3Geometry']
        self.m = {}
        self.time = 0.0
        self.PushFrame(frame)

    def Physics(self, frame):

        t0 = time.time()
        if not frame.Has(self.millipede):
            self.PushFrame(frame)
            return

        # Get millipede energy losses and length along track
        mill = frame[self.millipede]
        Eloss, length  = [],[]
        for item in mill:
            Eloss  += [item.energy]
            length += [item.length]
        Eloss = array(Eloss)
        length = cumsum(length)

        # strip zeros at beginning and end
        if Eloss.max()!=0:
            firstloss = (Eloss==Eloss[Eloss>0][0])
            start = firstloss.argmax()
            lastloss = (Eloss==Eloss[Eloss>0][-1])
            # make sure loss is unique, so we are finding the true end
            assert sum(lastloss)==1
            stop = lastloss.argmax()+1
            Eloss = Eloss[start:stop]
            length = length[start:stop] - length[0] + 10

        # Write to frame
        frame['millipede_Eloss']  = dc.I3VectorDouble(Eloss)
        frame['millipede_length'] = dc.I3VectorDouble(length)
        self.time += time.time() - t0
        self.PushFrame(frame)
        return

    def Finish(self):
        print 'Extractor time:', self.time
        return

