#!/usr/bin/env python

#############################################################################
# Records the energy, zenith, distances, snow depths, and charges for the
# creation of the LLH histogram
#############################################################################

from icecube import icetray, phys_services
from icecube import dataclasses as dc
from numpy import *

class fillHist(icetray.I3Module):

    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        self.AddOutBox('OutBox')
        self.AddParameter('binDict','Input LLH bins', 0)
        self.AddParameter('recoPulses','Input RecoPulseSeries to use for reconstruction', "CleanedHLCTankPulses")
        self.AddParameter('outFile', 'Input name for output file', 'test')

    def Configure(self):
        self.binDict = self.GetParameter('binDict')
        self.recoPulses = self.GetParameter('recoPulses')
        self.outFile = self.GetParameter('outFile')
        pass

    def Geometry(self, frame):
        self.geometry = frame['I3Geometry']
        self.om2index = {}
        self.tankpositions = []
        self.snowheight = []
        n = {}
        for key in self.binDict.keys():
            n[key[0]] = len(self.binDict[key]) - 1
        self.binLengths = (n['E'],n['Z'],n['S'],n['D'],n['C'])
        self.hist = zeros(self.binLengths)
        self.PushFrame(frame)

    def DetectorStatus(self, frame):
        self.status = frame['I3DetectorStatus']
        index = 0
        for number, station in self.geometry.stationgeo:
            for tank in station:
                om1, om2 = tank.omkey_list
                dom = self.status.dom_status
                # don't add tank unless at least one dom is in dom_status
                if om1 not in dom.keys() or om2 not in dom.keys():
                    continue
                # make sure doms are on
                elif dom[om1].pmt_hv > 5e-7 and dom[om2].pmt_hv > 5e-7:
                    self.om2index[om1] = index
                    self.om2index[om2] = index
                    self.tankpositions.append(tank.position)
                    #get snow height, if it changes there will be a gcd file
                    self.snowheight.append(tank.snowheight)
                    index += 1

        if len(self.tankpositions) < 5: # or some arbitrary low number
            raise ValueError("Not enough tanks found!")

        self.snowheight = array(self.snowheight)
        self.PushFrame(frame)

    def Physics(self, frame):

        # Load LLH bins   (E x Z x S x D x C)
        bins = {}
        for key in self.binDict.keys():
            bins[key[0]] = self.binDict[key]
        SDCedges = (bins['S'], bins['D'], bins['C'])

        # Get shower information
        try:  # the first frame doesn't have a particle?
            MC = frame['I3MCTree'].primaries[0]
            VEMpulses = frame[self.recoPulses]
            if VEMpulses.__class__ == dc.I3RecoPulseSeriesMapMask:
                VEMpulses = VEMpulses.apply(frame)
        except KeyError:
            self.PushFrame(frame)
            return

        ClosestDist = phys_services.I3Calculator.closest_approach_distance
        dists = array([ClosestDist(MC, xyz) for xyz in self.tankpositions])
        # Calculate charges in tanks
        VEMcharges = zeros(self.snowheight.shape)
        for om, pulses in VEMpulses:
            try:
                VEMcharges[self.om2index[om]] = pulses[0].charge
            except KeyError:    # skip OMKey(39,62)
                continue

        # Pick out the tanks with valid charges
        goodcharges = logical_not(isnan(VEMcharges))
        snows = self.snowheight[goodcharges]
        dists = dists[goodcharges]
        charges = VEMcharges[goodcharges]

        # Bin shower
        Ebin = digitize([log10(MC.energy)], bins['E'])[0] - 1
        Zbin = digitize([MC.dir.zenith], bins['Z'])[0] - 1
        SDChist = histogramdd((snows, dists, charges), SDCedges)[0]
        self.hist[Ebin,Zbin] += SDChist
        self.PushFrame(frame)
        return

    def Finish(self):
        if self.outFile:
            save(self.outFile, self.hist)
        return







