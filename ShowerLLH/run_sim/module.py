#!/usr/bin/env python

from icecube import icetray
#from icecube.InIceLLH import InIceLLHFitParams
from icecube import dataclasses as dc
from icecube.phys_services import I3Calculator as calc
from numpy import *
from copy import deepcopy
#from SupportFunctions import b2e, getLLH
import time

"""
class GetLLH(icetray.I3ConditionalModule):

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddOutBox('OutBox')
        self.AddParameter('binDict','Input LLH bins', 0)
        self.AddParameter('LLHTables', 'Input LLH tables', 0)
        self.AddParameter('millipede', 'Extracted data from millipede OR \
                Name of millipede particle vector', '')
        self.AddParameter('pulses', 'Name of pulses to use for llh', '')
        self.AddParameter('track', 'Name of track to use for llh', '')
        self.AddParameter('e_seed', 'Name of track to pull energy from \
                (defaults to ShowerLLH)', True)
        self.AddParameter('TRange', 'array of [tmin,tmax] time residual bounds\
                to use for Tres llh','')
        self.AddParameter('DRange', 'array of [dmin,dmax] track-DOM distance\
                bounds to use with Tres LLH','')

    def Configure(self):
        self.binDict = self.GetParameter('binDict')
        self.LLHTables = self.GetParameter('LLHTables')
        self.millipede = self.GetParameter('millipede')
        self.pulses = self.GetParameter('pulses')
        self.track = self.GetParameter('track')
        self.e_seed = self.GetParameter('e_seed')
        self.tmin = self.GetParameter('TRange')[0]
        self.tmax = self.GetParameter('TRange')[1]
        self.dmin = self.GetParameter('DRange')[0]
        self.dmax = self.GetParameter('DRange')[1]
        self.t0 = time.time()
        pass

    def Geometry(self, frame):
        self.geometry = frame['I3Geometry']
        self.PushFrame(frame)

    def Physics(self, frame):

        # Make sure frame has MPEFit and is downgoing
        if not frame.Has(self.track):
            self.PushFrame(frame)
            return
        if frame[self.track].dir.zenith > pi/2:
            self.PushFrame(frame)
            return

        ### Not sure if this is necessary anymore ###
        # Get MCPrimary information
        #if frame.Has('I3MCTree'):
        #    MCPrimary = frame['I3MCTree'].primaries
        #    if len(MCPrimary) > 1:
        #        print 'Warning, there are multiple seed particles'
        #    MCPrimary = MCPrimary[0]
        #elif frame.Has('MCPrimary'):
        #    MCPrimary = frame['MCPrimary']
        #else:
        #    self.PushFrame(frame)
        #    return

        # Limit energy range to IceTop simulation
        #if (log10(MCPrimary.energy)<4) or (log10(MCPrimary.energy)>9.5):
        #    self.PushFrame(frame)
        #    return

        # Get energy seed information (if applicable)
        seed_check = False
        if self.e_seed:
            if self.e_seed == True:
                p_seed = 'ShowerLLH_proton'
                f_seed = 'ShowerLLH_iron'
            else:
                p_seed, f_seed = self.e_seed, self.e_seed

            if frame.Has(p_seed):
                seed_check = True
                P_seed = log10(frame[p_seed].energy)
                F_seed = log10(frame[f_seed].energy)
                pEbin = int(digitize([P_seed], self.binDict['Ebins'])[0] - 1)
                fEbin = int(digitize([F_seed], self.binDict['Ebins'])[0] - 1)

        # Get millipede energy losses and length along track
        if isinstance(self.millipede, str):
            if frame.Has('millipede_Eloss'):
                Eloss  = asarray(frame['millipede_Eloss'])
                length = asarray(frame['millipede_length'])
            else:
                self.PushFrame(frame)
                return
        else:
            RunID    = frame['I3EventHeader'].run_id
            SubRunID = frame['I3EventHeader'].sub_run_id
            EventID  = frame['I3EventHeader'].event_id
            event_name = '%i %i %i' % (RunID, SubRunID, EventID)
            if event_name not in self.millipede.keys():
                self.PushFrame(frame)
                return
            Eloss  = self.millipede[event_name]['Eloss']
            length = self.millipede[event_name]['length']

        # Calculate energy proxy - dEdX
        EDep = Eloss.sum()
        if EDep<=0.1:
            EDep = 0.11
        index = (Eloss>0)
        if len(Eloss[index])>0:
            length_tot = length[index][-1] - length[index][0] + 10
        else:
            length_tot = max(length)
        dEdX = EDep / length_tot
        if dEdX<=0.1:
            dEdX = 0.11

        # Calculate time residual
        pulses = dc.I3RecoPulseSeriesMap.from_frame(frame, self.pulses)
        #pulses = frame[self.pulses]
        track = frame[self.track]
        geo = self.geometry.omgeo
        tres, dist = [],[]
        for dom in pulses.keys():
            if pulses[dom]:
                pos = geo[dom].position
                t = pulses[dom][0].time
                tres += [calc.time_residual(track, pos, t)]
                dist += [calc.closest_approach_distance(track, pos)]
        tres = array(tres)
        dist = array(dist)

        # Restrict to DOMs in a certain range
        tInd = (tres>self.tmin) & (tres<self.tmax)
        dInd = (dist>self.dmin) & (dist<self.dmax)
        tres = tres[(tInd)&(dInd)]
        dist = dist[(tInd)&(dInd)]

        # Bin shower
        Zbin = int(digitize([track.dir.zenith], self.binDict['Zbins'])[0] - 1)
        dEdXbin = int(digitize([dEdX], self.binDict['dEdXbins'])[0] - 1)
        ebin = digitize(Eloss, self.binDict['ebins']).astype(int32) - 1
        lbin = digitize(length, self.binDict['lbins']).astype(int32) - 1
        if len(tres)!=0:
            tbin = digitize(tres, self.binDict['tbins']).astype(int32) - 1
            dbin = digitize(dist, self.binDict['dbins']).astype(int32) - 1

        # Load LLH tables     (comp x E x Z x dEdX x e/t x l/d)
        hists = self.LLHTables
        hList = hists.keys()
        # Only include TD if we have pulses within our range
        if len(tres)==0:
            hList.remove('TD')

        # Organize binned information into list
        binList = {}
        for h in hList:
            binList[h] = [Zbin, dEdXbin]
            if h == 'EL':
                binList[h] += [ebin, lbin]
            if h == 'TD':
                binList[h] += [tbin, dbin]

        if seed_check:
            seedList = {'p':{},'f':{}}
            for h in hList:
                seedList['p'][h] = binList[h][:]
                seedList['p'][h].insert(0, pEbin)
                seedList['f'][h] = binList[h][:]
                seedList['f'][h].insert(0, fEbin)

        # Find most likely energy bin and store as track
        track, maxLLH, maxLLH_pseed, maxLLH_fseed = {},{},{},{}
        for h in hList:
            compList = hists[h].keys()
            compList.sort()
            track[h], maxLLH[h] = {},{}
            maxLLH_pseed[h], maxLLH_fseed[h] = {},{}

            for e in compList:
                # Setup I3Particle with position and direction
                track[h][e] = dc.I3Particle()
                track[h][e].fit_status = dc.I3Particle.OK
                maxLLH[h][e], Ebin = getLLH(hists[h][e], binList[h])
                if seed_check:
                    maxLLH_pseed[h][e] = getLLH(hists[h][e], seedList['p'][h])
                    maxLLH_fseed[h][e] = getLLH(hists[h][e], seedList['f'][h])
                track[h][e].energy = b2e(self.binDict['Ebins'], Ebin)

        ## Setup InIceLLH parameters for output ##
        params = {}
        for h in hList:
            params[h] = InIceLLHFitParams()
            params[h].RecoRanCut = bool(1)
            params[h].maxLLH_proton   = float(maxLLH[h]['P'])
            #params[h].maxLLH_helium   = float(maxLLH[h]['He'])
            #params[h].maxLLH_nitrogen = float(maxLLH[h]['N'])
            #params[h].maxLLH_oxygen = float(maxLLH[h]['O'])
            #params[h].maxLLH_aluminum = float(maxLLH[h]['Al'])
            params[h].maxLLH_iron     = float(maxLLH[h]['Fe'])

            if seed_check:
                params[h].maxLLH_proton_pseed = float(maxLLH_pseed[h]['P'])
                params[h].maxLLH_proton_fseed = float(maxLLH_fseed[h]['P'])
                params[h].maxLLH_iron_pseed = float(maxLLH_pseed[h]['Fe'])
                params[h].maxLLH_iron_fseed = float(maxLLH_fseed[h]['Fe'])


        ## WRITE TO FRAME ##
        for h in hList:
            frame['InIceLLH_'+h+'_proton'] = track[h]['P']
            #frame['InIceLLH_'+h+'_helium'] = track[h]['He']
            #frame['InIceLLH_'+h+'_nitrogen'] = track[h]['N']
            #frame['InIceLLH_'+h+'_oxygen'] = track[h]['O']
            #frame['InIceLLH_'+h+'_aluminum'] = track[h]['Al']
            frame['InIceLLH_'+h+'_iron'] = track[h]['Fe']
            frame['InIceLLHParams_'+h] = params[h]


        self.PushFrame(frame)
        return
"""


class MergeIIIT(icetray.I3PacketModule):

    def __init__(self, ctx):
        icetray.I3PacketModule.__init__(self, ctx, icetray.I3Frame.DAQ)
        self.AddParameter("IceTopReco", "IceTop reconstruction to look for", "")
        self.AddParameter("InIceReco", "InIce reconstruction to look for", "")
        self.AddParameter("IceTopPulses", "IceTop pulses to look for", "")
        self.AddParameter("InIcePulses", "InIce pulses to look for", "")
        self.AddOutBox("OutBox")

    def Configure(self):
        self._icetop_reco = self.GetParameter("IceTopReco")
        #self._inice_reco = self.GetParameter("InIceReco")
        self._icetop_pulses = self.GetParameter("IceTopPulses")
        self._inice_pulses = self.GetParameter("InIcePulses")
        if len(self._icetop_reco) == 0:# or len(self._inice_reco) == 0:
            raise RuntimeError("Requires both IceTop and InIce reconstructions")
        #if len(self._inice_reco) == 0:
        #    raise RuntimeError("Requires InIce reconstruction")
        if len(self._icetop_pulses) == 0 or len(self._inice_pulses) == 0:
            raise RuntimeError("Requires both IceTop and InIce pulses")

    def _get_event_size(self, frame, pulses_name):
        # Account for I3RecoPulseSeriesMap and I3RecoPulseSeriesMapMask objects
        try:
            pulses = frame[pulses_name]
            pulse_values = pulses.values()
        except AttributeError:
            pulses = frame[pulses_name].apply(frame)
            pulse_values = pulses.values()
        npulses = 0
        for pulse_series in pulse_values:
            npulses += len(pulse_series)
        return npulses

    def FramePacket(self, frames):

        if len(frames) <= 1:            # only Q frame
            for frame in frames:
                self.PushFrame(frame)
            return

        self.PushFrame(frames[0])       # Push Q frame

        icetop_frame, inice_frame, merged_frame = False, False, False
        icetop_largest_size, inice_largest_size = 0,0
        event_id = frames[1]['I3EventHeader'].event_id

        for frame in frames[1:]:
            
            # Get largest icetop frame
            if self._icetop_pulses in frame:
                event_size = self._get_event_size(frame, self._icetop_pulses)
                if event_size > icetop_largest_size:
                    icetop_frame = frame
                    icetop_largest_size = event_size
            # Get largest inice frame
            if self._inice_pulses in frame:#self._inice_reco in frame:
                event_size = self._get_event_size(frame, self._inice_pulses)
                if event_size > inice_largest_size:
                    inice_frame = frame
                    inice_largest_size = event_size
            # Make sure merged frame doesn't already exist
            #if self._icetop_reco in frame and self._inice_reco in frame:
            #    merged_frame = True

            self.PushFrame(frame)
        # Create merged frame
        if icetop_frame and inice_frame: # and not merged_frame:
            new_frame = icetray.I3Frame(icetray.I3Frame.Physics)
            new_frame.merge(inice_frame)
            eh = dc.I3EventHeader(new_frame['I3EventHeader'])
            new_frame.Delete('I3EventHeader')
            eh.sub_event_id = 0
            eh.sub_event_stream = self.name
            new_frame.merge(icetop_frame)
            new_frame.Delete('I3EventHeader')
            new_frame['I3EventHeader'] = eh
            self.PushFrame(new_frame)

        elif icetop_frame and not inice_frame:
            new_frame = icetray.I3Frame(icetray.I3Frame.Physics)
            new_frame.merge(icetop_frame)
            eh = dc.I3EventHeader(new_frame['I3EventHeader'])
            eh.sub_event_id = 0
            eh.sub_event_stream = 'icetop_only' 
            new_frame.Delete('I3EventHeader')
            new_frame['I3EventHeader'] = eh
            self.PushFrame(new_frame)

        #Currently not wanted
        '''
        if inice_frame and not icetop_frame:
            new_frame = icetray.I3Frame(icetray.I3Frame.Physics)
            new_frame.merge(inice_frame)
            eh = dc.I3EventHeader(new_frame['I3EventHeader'])
            new_frame.Delete('I3EventHeader')
            eh.sub_event_id = 0
            eh.sub_event_stream = 'icecube_only' 
            new_frame['I3EventHeader'] = eh
            self.PushFrame(new_frame)
        '''

class MergeIIIT_simple(icetray.I3PacketModule):

    def __init__(self, ctx):
        icetray.I3PacketModule.__init__(self, ctx, icetray.I3Frame.DAQ)
        self.AddParameter("IceTopReco", "IceTop reconstruction to look for", "")
        self.AddParameter("IceTopPulses", "IceTop pulses to look for", "")
        self.AddParameter("InIcePulses", "InIce pulses to look for", "")
        self.AddOutBox("OutBox")

    def Configure(self):
        self._icetop_reco = self.GetParameter("IceTopReco")
        self._icetop_pulses = self.GetParameter("IceTopPulses")
        self._inice_pulses = self.GetParameter("InIcePulses")
        if len(self._icetop_reco) == 0: 
            raise RuntimeError("Requires IceTop reconstruction")
        if len(self._icetop_pulses) == 0 or len(self._inice_pulses) == 0:
            raise RuntimeError("Requires both IceTop and InIce pulses")

    def _get_event_size(self, frame, pulses_name):
        # Account for I3RecoPulseSeriesMap and I3RecoPulseSeriesMapMask objects
        try:
            pulses = frame[pulses_name]
            pulse_values = pulses.values()
        except AttributeError:
            pulses = frame[pulses_name].apply(frame)
            pulse_values = pulses.values()
        npulses = 0
        for pulse_series in pulse_values:
            npulses += len(pulse_series)
        return npulses

    def FramePacket(self, frames):

        if len(frames) <= 1:            # only Q frame
            for frame in frames:
                self.PushFrame(frame)
            return

        self.PushFrame(frames[0])       # Push Q frame

        icetop_frame, inice_frame, merged_frame = False, False, False
        icetop_largest_size, inice_largest_size = 0,0
        event_id = frames[1]['I3EventHeader'].event_id

        for frame in frames[1:]:
            
            # Get largest icetop frame
            if self._icetop_pulses in frame:
                event_size = self._get_event_size(frame, self._icetop_pulses)
                if event_size > icetop_largest_size:
                    icetop_frame = frame
                    icetop_largest_size = event_size
            # Get largest inice frame
            if self._inice_pulses in frame:
                event_size = self._get_event_size(frame, self._inice_pulses)
                if event_size > inice_largest_size:
                    inice_frame = frame
                    inice_largest_size = event_size
            # Make sure merged frame doesn't already exist
            #if self._icetop_reco in frame and self._inice_reco in frame:
            #    merged_frame = True

            self.PushFrame(frame)
        # Create merged frame
        if icetop_frame and inice_frame: # and not merged_frame:
            new_frame = icetray.I3Frame(icetray.I3Frame.Physics)
            new_frame.merge(inice_frame)
            eh = dc.I3EventHeader(new_frame['I3EventHeader'])
            new_frame.Delete('I3EventHeader')
            eh.sub_event_id = 0
            eh.sub_event_stream = self.name
            new_frame.merge(icetop_frame)
            new_frame.Delete('I3EventHeader')
            new_frame['I3EventHeader'] = eh
            self.PushFrame(new_frame)

        elif icetop_frame and not inice_frame:
            new_frame = icetray.I3Frame(icetray.I3Frame.Physics)
            new_frame.merge(icetop_frame)
            eh = dc.I3EventHeader(new_frame['I3EventHeader'])
            eh.sub_event_id = 0
            eh.sub_event_stream = 'icetop_only' 
            new_frame.Delete('I3EventHeader')
            new_frame['I3EventHeader'] = eh
            self.PushFrame(new_frame)

        #Currently not wanted
        '''
        if inice_frame and not icetop_frame:
            new_frame = icetray.I3Frame(icetray.I3Frame.Physics)
            new_frame.merge(inice_frame)
            eh = dc.I3EventHeader(new_frame['I3EventHeader'])
            new_frame.Delete('I3EventHeader')
            eh.sub_event_id = 0
            eh.sub_event_stream = 'icecube_only' 
            new_frame['I3EventHeader'] = eh
            self.PushFrame(new_frame)
        '''
