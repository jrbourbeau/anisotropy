#!/usr/bin/env python

#############################################################################
# Records the energy, zenith, distances, snow depths, and charges for the
# creation of the LLH histogram - mimics C++ code in showerllh project
#############################################################################

from icecube import icetray, phys_services
from icecube import dataclasses as dc
import numpy as np

class fillHist(icetray.I3Module):

    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        self.AddOutBox('OutBox')
        self.AddParameter('binDict','Input LLH bins', 0)
        self.AddParameter('recoPulses','Input RecoPulseSeries to use for reconstruction', "CleanedHLCTankPulses")
        self.AddParameter('outFile', 'Input name for output file', 'test')

    def Configure(self):
        self.bins = self.GetParameter('binDict')
        self.recoPulses = self.GetParameter('recoPulses')
        self.outFile = self.GetParameter('outFile')
        pass

    def Geometry(self, frame):
        self.geometry = frame['I3Geometry']
        self.om2index = {}
        self.tankpositions = []
        self.snowheight = []
        self.PushFrame(frame)

    def Calibration(self, frame):
        self.calibration = frame['I3Calibration']
        self.PushFrame(frame)

    def DetectorStatus(self, frame):
        self.status = frame['I3DetectorStatus']
        self.dom = self.status.dom_status
        index = 0
        for number, station in self.geometry.stationgeo:
            for tank in station:
                omList = tank.omkey_list
                if not any([om in self.dom.keys() for om in omList]):
                    continue
                for om in omList:
                    self.om2index[om] = index
                self.tankpositions.append(tank.position)
                self.snowheight.append(tank.snowheight)
                index += 1

        if len(self.tankpositions) < 5: # or some arbitrary low number
            raise ValueError("Not enough tanks found!")

        self.hist = np.zeros([len(self.bins[i][1])-1 for i in self.bins])
        self.snowheight = np.array(self.snowheight)
        self.PushFrame(frame)

    def Physics(self, frame):

        # Get shower information
        try:  # the first frame doesn't have a particle?
            MC = frame['I3MCTree'].primaries[0]
            VEMpulses = frame[self.recoPulses]
            if VEMpulses.__class__ == dc.I3RecoPulseSeriesMapMask:
                VEMpulses = VEMpulses.apply(frame)
        except KeyError:
            self.PushFrame(frame)
            return

        # Closest approach distances
        ClosestDist = phys_services.I3Calculator.closest_approach_distance
        dists = np.array([ClosestDist(MC, xyz) for xyz in self.tankpositions])

        # Saturation values
        sat_lg_pe = 90000.
        vemcal_map = self.calibration.vem_cal

        # Calculate charges in tanks
        VEMcharges = np.zeros(self.snowheight.shape)
        saturated = np.zeros(self.snowheight.shape, dtype='bool')
        for om, pulses in VEMpulses:

            # Skip OMKey(39,62)
            if om not in self.om2index.keys():
                continue

            # Only record first pulse
            charge = pulses[0].charge
            # CleanedHLCTankPulses has one charge per tank
            idx = self.om2index[om]
            VEMcharges[idx] = charge

            # Check saturation
            vemCalib = vemcal_map[om]
            pe_per_vem = vemCalib.pe_per_vem / vemCalib.corr_factor
            lg_sat = sat_lg_pe / pe_per_vem
            gain = self.dom[om].dom_gain_type
            if (charge > lg_sat) and (gain == dc.I3DOMStatus.Low):
                saturated[idx] = True

        # Pick out the tanks with valid charges
        goodcharges = np.logical_not(np.isnan(VEMcharges))
        # Set saturated DOMs to fixed value
        binNames = [self.bins[i][0] for i in sorted(self.bins)]
        c_idx = binNames.index('C')
        cbins = self.bins[c_idx][1]
        # Make sure they'll fall in last bin
        VEMcharges[goodcharges*saturated] = 1.01 * cbins[-2]

        unbinned_vals = {}
        unbinned_vals['E'] = [np.log10(MC.energy)] * goodcharges.sum()
        unbinned_vals['Z'] = [MC.dir.zenith] * goodcharges.sum()
        unbinned_vals['S'] = self.snowheight[goodcharges]
        unbinned_vals['D'] = dists[goodcharges]
        unbinned_vals['C'] = VEMcharges[goodcharges]

        # Bin shower
        all_unbinned = [unbinned_vals[k] for k in binNames]
        all_edges = [self.bins[i][1] for i in sorted(self.bins)]
        event_hist = np.histogramdd(all_unbinned, all_edges)[0]
        self.hist += event_hist
        self.PushFrame(frame)
        return

    def Finish(self):
        # Record count tables and bin values in outFile
        if self.outFile:
            d = {}
            d['bins']   = self.bins
            d['counts'] = self.hist
            np.save(self.outFile, d)
        return







