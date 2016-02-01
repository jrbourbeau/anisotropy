## TODO at some point : all cleaners of launches or pulses should output Masks..

from icecube import icetray,dataio,dataclasses,toprec, phys_services, simclasses
import math

#########################
# Python snowservice :  #
#########################
import csv
class ChangeSnowHeight(icetray.I3Module):
    def __init__(self, ctx):
        icetray.I3Module.__init__(self,ctx)
        self.AddParameter('filename','Filename w/ snow values for each tank', 0)
        self.AddOutBox("OutBox")

    def Configure(self):
        self.snowfile = self.GetParameter('filename')
        self.newheights = self._readFile(self.snowfile)

    def Geometry(self, frame):
        if 'I3Geometry' in frame:
            geom = frame['I3Geometry']
            stageo = geom.stationgeo
            for e,st in stageo:
                updated_heights = dataclasses.I3StationGeo()
                if not self.newheights.has_key(e):
                    print 'Did not find station ' , e, ' in new snowheight dict'
                    continue
                #ok we have I3TankGeo here... look it up
                ## assume first is A and second is B
                st[0].snowheight = max((self.newheights[e][0],0.))
                st[1].snowheight = max((self.newheights[e][1],0.))
                updated_heights.append(st[0])
                updated_heights.append(st[1])
                stageo[e] = updated_heights
            del frame['I3Geometry']
            frame['I3Geometry'] = geom
        else:
            print 'No geometry found'
        self.PushFrame(frame,"OutBox")


    def _readFile(self, filename):
        file = open(filename)
        ## error catching should go here...

        newheights = {}
        for row in csv.reader(file,delimiter=' '):
            # Deal with Takao's files:
            if len(row) > 1 :        ## skips first line
                height = float(row[len(row)-1])
                tank = row[0]
                station = int(tank[:2])
                id = tank[2]
                if station not in newheights:
                    newheights[station] = [0,0]
                if id == 'A':
                    newheights[station][0] = height
                elif id == 'B':
                    newheights[station][1] = height
                else:
                    print 'unknown tankid'
        return newheights

###############################################################################
# 0) BUG in Thinned MC : unsaturated tanks close to the core (within a
# station!!)              #
#                        --> NOT seen in non-thinned MC (overlap Ebin)
#                        #
#                        --> NOT seen in data
#                        #
# Showed up in log(Q) vs R(to true core) plot as a separate distribution :
# #
# lower half of crocodile mouth
# #
# ==> Based on that plot : E0 dependent cut which removes this unnatural mouth
# and 'Fixes' MC  #
# ==> SHOULD be fixed in dethinning for topsimulator though...
# #
###############################################################################
class cleanBadThinning(icetray.I3Module):
    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        self.AddParameter('InputITpulseName',
                'Which IceTop pulses to clean out',0)
        self.AddParameter('OutputPulseName',
                'Name of the resulting thinFix cleaned IT pulses',0)
        self.AddParameter('ExcludedPulseName',
                'Name of the resulting removed buggy IT pulses',0)
        self.AddOutBox("OutBox")

    def Configure(self):
        self.itpulses = self.GetParameter('InputITpulseName')
        self.outpulseName = self.GetParameter('OutputPulseName')
        self.exclpulseName = self.GetParameter('ExcludedPulseName')

    def Geometry(self, frame):
        if 'I3Geometry' in frame:
            geo = frame['I3Geometry']
            self.omg = geo.omgeo
        else:
            print 'No geometry found'
        self.PushFrame(frame,"OutBox")

    def Physics(self, frame):
        if self.itpulses in frame:
            itpulses = frame[self.itpulses]
            if itpulses.__class__ == dataclasses.I3RecoPulseSeriesMapMask:
                itpulses = itpulses.apply(frame)
            geo = frame["I3Geometry"]
            mcprim = frame["MCPrimary"]
            outputPulses = dataclasses.I3RecoPulseSeriesMap()
            excludedPulses = dataclasses.I3RecoPulseSeriesMap()
            for om,pulses in itpulses:
                  if om in geo.omgeo:
                      ompos = self.omg[om]
                      dist = math.sqrt( 
            (mcprim.pos.x - ompos.position.x)*(mcprim.pos.x - ompos.position.x)
           +(mcprim.pos.y - ompos.position.y)*(mcprim.pos.y - ompos.position.y))
                      for pulse in pulses:
                          exclude = False
                          if (dist < 200):
                              if math.log10(mcprim.energy) < 8.0 and math.log10(mcprim.energy) >= 7.0 :
                                  if(math.log10(pulse.charge) < (-3.35/95. * dist +3.35)):
                                      exclude = True
                                  else:
                                      exclude = False
                              elif math.log10(mcprim.energy) >= 8.0 and math.log10(mcprim.energy) < 8.5:
                                  if(math.log10(pulse.charge) < (-2.85/150. * dist +3.35)):
                                      exclude = True
                                  else:
                                      exclude = False
                              elif math.log10(mcprim.energy) >= 8.5 and math.log10(mcprim.energy) < 9.0:
                                  if(math.log10(pulse.charge) < (-2.4/185. * dist +3.5)):
                                      exclude = True
                              elif math.log10(mcprim.energy) >= 9.0:
                                  if(math.log10(pulse.charge) < (-1.7/200. * dist +3.55)):
                                      exclude = True
                              else:
                                  exclude = False
                          else:
                              exclude = False
                          ## put them in their maps
                          if exclude:
                              if om not in excludedPulses:
                                  excludedPulses[om] = dataclasses.I3RecoPulseSeries()
                              excludedPulses[om].append(pulse)
                          else:
                              if om not in outputPulses:
                                  outputPulses[om] = dataclasses.I3RecoPulseSeries()
                              outputPulses[om].append(pulse)
                  else:
                      print "OM not found in omgeo, wtf? : ", om

            # outputting!
            frame[self.outpulseName] = outputPulses
            frame[self.exclpulseName] = excludedPulses
        else:
            print "Missing inputITpulses : ", self.itpulses
        self.PushFrame(frame,"OutBox")


##########################################
# 1) Move MCPrimary from P to Q frame   #
#########################################
class moveMCPrimary(icetray.I3PacketModule):
    def __init__(self, ctx):
        icetray.I3PacketModule.__init__(self, ctx, icetray.I3Frame.DAQ)
        self.AddOutBox("OutBox")
   
    def Configure(self):
        pass

    def FramePacket(self, frames):
        if len(frames) <= 1:     # only q frame
            for frame in frames:
                self.PushFrame(frame)
            return

        qframe = frames[0]
        pframes = frames[1:]
        prim = 0
        prim_info = 0

        for frame in pframes:
            if 'MCPrimary' in frame:
                if prim != 0:
                    raise RuntimeError("Sorry, but this module is pretty dumb. It cannot treat cases where the MCPrimary appears in more than one P frame.")
                prim = frame['MCPrimary']
                if not 'MCPrimaryInfo' in frame:
                    print prim.energy, frame
                prim_info = frame['MCPrimaryInfo']
                del frame['MCPrimary']
                del frame['MCPrimaryInfo']

        qframe['MCPrimary'] = prim
        qframe['MCPrimaryInfo'] = prim_info

        self.PushFrame(qframe,"OutBox")
        for frame in pframes:
            self.PushFrame(frame,"OutBox")


################################################################################
# 2. A FilterMaskFilter that not only drops P frames, but the whole packet (to
# prevent modules   #
#    that work on Q frame objects to do so for Orphan Qframes)
#    #
#    because FilterMask in the Q frame for IT73/IC79 data (do not use for
#    future IT81/IC86 data) #
#    where FilterMask will be in the P frame from the beginning
#    #
################################################################################
class FilterMaskFilter(icetray.I3PacketModule):
    def __init__(self, ctx):
        icetray.I3PacketModule.__init__(self, ctx, icetray.I3Frame.DAQ)
        self.AddParameter("FilterNameList","(DAQ) with which filters do you want to keep",[])
        self.AddParameter("FilterResultName","Name of FilterMaskList in the Q frame","")
        self.AddOutBox("OutBox")

    def Configure(self):
        self.filterNames = self.GetParameter("FilterNameList")
        self.filterResultName = self.GetParameter("FilterResultName")
        self.nTotal = 0
        self.nKept = 0

    def FramePacket(self, frames):
        self.nTotal +=1
        filterResult=[]
        # frames : QPPPP
        qframe = frames[0]
        #print qframe
        if qframe.Has(self.filterResultName):
            filterResult = qframe[self.filterResultName]
            #print filterResult
            for name in self.filterNames:
                if name in filterResult:
                    #print filterResult[name]
                    if filterResult[name].condition_passed and filterResult[name].prescale_passed:
                        # push all frames
                        #print "Filter condition met, pushing all frames 
                        #in the packet"
                        self.PushFrame(qframe)
                        for pframe in frames[1:]:
                            self.PushFrame(pframe)
                        self.nKept += 1
                        return
    def Finish(self):
        # only for PacketModules (normal modules don't have queueing)
        self.FlushQueue()
        print 'Number of events processed : ', self.nTotal
        if self.nTotal > 0:
            print 'Number of events kept after Coincident Filter : %i (of total : %lf)' %(self.nKept,float(self.nKept)/float(self.nTotal))

####################################################################
## 3) My own little selector : NSTA > 5 OR NCH > 8 : output a bool #
####################################################################
def SelectStationsChannels(frame, InputITpulses, InputIIpulses, NStaThreshold, NChThreshold, outputName, RemoveEvents):
    nstation = 0
    nch = 0
    if InputITpulses in frame:
        vemPulses = frame[InputITpulses]
        if vemPulses.__class__ == dataclasses.I3RecoPulseSeriesMapMask:
            vemPulses = vemPulses.apply(frame)
        stationList = set([pulse.key().string for pulse in vemPulses])
        nstation = len(stationList)
    if InputIIpulses in frame:
        inIcePulses = frame[InputIIpulses]
        if inIcePulses.__class__ == dataclasses.I3RecoPulseSeriesMapMask:
            inIcePulses = inIcePulses.apply(frame)
        # only count HLC !!
        for om,pulses in inIcePulses:
            for pulse in pulses:
                #print pulse.flags, " LC? ", pulse.flags & dataclasses.I3RecoPulse.PulseFlags.LC
                ## use bit-wise operation in python to check if the last bit is 1 (LC)
                if (pulse.flags & dataclasses.I3RecoPulse.PulseFlags.LC) : 
                    nch += 1
                    break  # only count 1 pulse per DOM if it isn't an SLC pulse
    #outputName = "%s_NChNStaCut" % InputIIpulses
    toCut = (nstation > NStaThreshold) and (nch > NChThreshold) 
    if not RemoveEvents :
        frame[outputName] = icetray.I3Bool(toCut)
        return True
    else:
        return toCut



####################################################################
## 3a) Modification: just output NSTA                              #
####################################################################
def GetStations(frame, InputITpulses, OutputName):
    nstation = 0
    if InputITpulses in frame:
        vemPulses = frame[InputITpulses]
        if vemPulses.__class__ == dataclasses.I3RecoPulseSeriesMapMask:
            vemPulses = vemPulses.apply(frame)
        stationList = set([pulse.key().string for pulse in vemPulses])
        nstation = len(stationList)
    frame[OutputName] = icetray.I3Int(nstation)
    return



##############################################################
## 4. Find the loudest station and station with loudest tank #
##############################################################
class FindLoudestStation(icetray.I3Module):
    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        self.AddParameter('InputITpulses','Which IceTop Pulses to use',0)
        self.AddParameter('SaturationValue','Tanks with charges above this values (in VEM) are saturated',600)
        self.AddParameter('OutputSaturatedStationsName','Name of vector with saturated stations (replaces loudeststationoutput if specified','')
        self.AddOutBox("OutBox")

    def Configure(self):
        self.pulses = self.GetParameter('InputITpulses')
        self.satu = self.GetParameter('SaturationValue')
        self.listName = self.GetParameter('OutputSaturatedStationsName')

    def Physics(self, frame):
        if self.pulses in frame:
            vem = frame[self.pulses]
            if vem.__class__ == dataclasses.I3RecoPulseSeriesMapMask:
                vem = vem.apply(frame)
            loudPulse = 0
            staCharge = 0
            loudStaCharge = 0
            avStaCharge = 0
            loudStation1 = 0
            loudStation2 = 0
            loudStation3 = 0
            prevSta = 0
            sat_stations = []
            for key,series in vem:
                for tank in series: #will be one waveform anyway

                    # if NaN : rely on waveform in other tank, so skip
                    if tank.charge != tank.charge:
                        continue

                    if tank.charge > loudPulse:
                        loudPulse = tank.charge
                        loudStation1 = key.string
                    if key.string != prevSta:
                        staCharge = tank.charge

                    elif key.string == prevSta:
                        staCharge += tank.charge
                        if staCharge > loudStaCharge:
                            loudStaCharge = staCharge
                            loudStation2 = key.string
                        if staCharge/2 > avStaCharge:
                            avStaCharge = staCharge/2
                            loudStation3 = key.string
                    prevSta = key.string

                    # LG saturaton bookkeeping :
                    if tank.charge > self.satu:
                        sat_stations.append(key.string)
            if self.listName == '':
               #Now put stuff in Frame
               frame["StationWithLoudestPulse"] = dataclasses.I3Double(loudStation1)
               frame["LoudestStation"] = dataclasses.I3Double(loudStation2)
               frame["LoudestAverageStation"] = dataclasses.I3Double(loudStation3)
            else :
                sat_stations = set(sat_stations)
                # Below : merging loudestStation and saturated stations in one
                # list
                # --> Preferably not, because makes interpretation of number of
                # saturation stations
                # difficult, so two outputs here
                ### if len(sat_stations) < 1:         ###
                ### sat_station.append(loudStation1)  ###
                sta_list = dataclasses.I3VectorInt()
                for sta in sat_stations:
                    sta_list.append(sta)
                frame[self.listName] = sta_list
                frame["StationWithLoudestPulse"] = dataclasses.I3Double(loudStation1)
                frame["LoudestStation"] = dataclasses.I3Double(loudStation2)
                frame["LoudestAverageStation"] = dataclasses.I3Double(loudStation3)

        #else :
        #    print "WARN : no IT pulses %s found, so no LoudestStation" % self.pulses
        self.PushFrame(frame,"OutBox")


################################################################################
## 5. Is the loudest station OR any saturated station on the edge for some
## configuration?  #
################################################################################
class LoudestStationOnEdge(icetray.I3Module):
    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        self.AddParameter('InputLoudestStation','Which Loudest Station (or I3Vector of saturated stations)',0)
        self.AddParameter('WhichDetector','Which detectorConfig? IT40/59/73/81? ',0)
        self.AddParameter('OutputName','What is the name of the output Bool in the frame',0)
        self.AddOutBox("OutBox")

    def Configure(self):
        self.loudest = self.GetParameter('InputLoudestStation')
        self.config = self.GetParameter('WhichDetector')
        self.outputName = self.GetParameter('OutputName')

    def Physics(self, frame):
        if self.loudest in frame:
            loud = frame[self.loudest]  ## is an I3Double or I3VectorInt
            edge = False
            edgelist = []
            if self.config == "IT40":
                edgelist = [21,30,40,50,
                            59,67,74,73,72,78,
                            77,76,75,68,
                            60,52,53,
                            44,45,46,47,38,29]   #47 is questionable (makes cut harder)
            elif self.config == "IT59":
                edgelist = [2,3,4,5,6,
                            13,21,30,40,50,
                            59,67,74,73,72,78,
                            77,76,75,
                            68,60,52,53,44,45,
                            36,26,17,9]
            elif self.config == "IT73":
                edgelist = [2,3,4,5,6,
                            13,21,30,40,50,
                            59,67,74,73,72,78,
                            77,76,75,
                            68,60,51,41,32,23,15,8]
            elif self.config == "IT81":
                edgelist = [2,3,4,5,6,
                            13,21,30,40,50,
                            59,67,74,73,72,78,
                            77,76,75,
                            68,60,51,41,31,22,14,7,1]
            else:
                print "Unknown Config, Please choose IT40,IT59,IT73,IT81"
                sys.exit(-1)

            for station in edgelist:     # the borderlist
                if loud.__class__ == dataclasses.I3Double:
                    if loud.value == station:
                        #print "We're on the edge of the array"
                        edge = True
                elif loud.__class__ == dataclasses.I3VectorInt:
                    # This is the situation for a list of saturated stations
                    if station in loud:
                        #print "Saturated station spotted on the edge of the array"
                        edge = True
            output = icetray.I3Bool(edge)
            frame[self.outputName] = output
        #else :
        #    print "WARN : no inputName %s found" % self.loudest
        self.PushFrame(frame,"OutBox")


##############################
# 6. For InIce Hitcleaning : #
##############################
class CleanNoiseLC(icetray.I3Module):
    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        self.AddParameter('InputDOMLaunches','Input Launches Name','SeededRTInIceRawData')
        self.AddParameter('OutputDOMLaunches','NoiseLC cleaned Launches Name','CleanedHLCRawData')
        self.AddParameter('MinLCtimeDifference','Average time difference between a DOM and all its LC neighbours must be larger than this value',-300)
        self.AddOutBox("OutBox")

    def Configure(self):
        self.inputlaunchName = self.GetParameter('InputDOMLaunches')
        self.outputlaunchName = self.GetParameter('OutputDOMLaunches')
        self.minTimeDiff = self.GetParameter('MinLCtimeDifference')

    def DAQ(self, frame):
        if self.inputlaunchName in frame:
            inputLaunches = frame[self.inputlaunchName]
            outputLaunches = dataclasses.I3DOMLaunchSeriesMap()
            for om, launches in inputLaunches:
                countLaunch = 1
                for launch in launches:
                    if launch.lc_bit:            # HLC
                        # map is organized per key, key per string and per DOM
                        # For each DOM, calculate the average time difference
                        # with all its neighbours
                        # if early noise => average will be negative!
                        av_time = 0
                        n_neighbours = 0
                        for i in (-2,-1,1,2): # Next-to_Nearest and Nearest neighbours if they're there
                            if om.om+i < 0:
                                continue
                            neighb_om = icetray.OMKey(om.string,om.om+i) ## one downwards
                            if neighb_om in inputLaunches:
                                for neighb_launch in inputLaunches[neighb_om]:
                                    if neighb_launch.lc_bit: 
                                        if abs(neighb_launch.time-launch.time) < 1001.:  # otherwise 2nd Launch, no LC
                                            n_neighbours += 1
                                            av_time += (launch.time-neighb_launch.time)

                        if n_neighbours > 0:
                            av_time /= n_neighbours
                            if av_time > self.minTimeDiff:
                                # keep HLC DOM
                                if om not in outputLaunches :
                                    outputLaunches[om] = dataclasses.I3DOMLaunchSeries()
                                outputLaunches[om].append(launch)
                            #else :
                            #    print om, " not kept because averageTime is ",
                            #    av_time

            frame[self.outputlaunchName] = outputLaunches
        self.PushFrame(frame,"OutBox")

################################################################################
# 7. Calculate the Largest n Pulses and the one in the neighbour tank for the
# largest one (Q1b) #
################################################################################
class LargestTankCharges(icetray.I3Module):
    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        self.AddParameter('nPulses','Book the largest N pulses for TailCut (+neighbour of largest)',4)
        self.AddParameter('ITpulses','IT pulses Name',0)
        self.AddOutBox("OutBox")

    def Configure(self):
        self.nPulses = self.GetParameter('nPulses')
        self.VEM = self.GetParameter('ITpulses')
        self.counter=0

    def Physics(self, frame):
        if self.VEM in frame:
            tank_map = frame[self.VEM]
            if tank_map.__class__ == dataclasses.I3RecoPulseSeriesMapMask:
                tank_map = tank_map.apply(frame)

            ## Fast way, but does NOT deal with NaNs
            ##charge_map = [(pulses[0].charge, dom) for dom,pulses in tank_map]

            ## If NaN charge, use charge in other tank as best approximation
            charge_map = []
            for om,pulses in tank_map:
                for wave in pulses:
                    charge = -1
                    #if math.isnan(wave.charge): #only works beyond python 2.6
                    if wave.charge != wave.charge:
                        #print "Charge is NaN, so saturated HG without LG (or DOM 39-62)" 
                        # Other tank will be next or previous element in map,
                        # use charge of other tank as rough approx
                        omkey = icetray.OMKey()
                        omkey.string = om.string

                        # Way too complicated way to get the omkey of the
                        # neighbouring tank
                        if om.om == 61 or om.om == 62:
                            omkey.om = om.om+2
                            if omkey not in tank_map:
                                if om.om == 61:
                                    omkey.om = om.om+3
                                else:
                                    omkey.om = om.om+1
                        elif om.om == 63 or om.om == 64:
                            omkey.om = om.om-2
                            if omkey not in tank_map:
                                if om.om == 63:
                                    omkey.om = om.om-1
                                else:
                                    omkey.om = om.om-3
                        #print omkey
                        if omkey in tank_map:
                            pulses = tank_map[omkey]
                            charge = pulses[0].charge
                        else:
                            # in case some pulse cleaning cleaned one tank of
                            # the station
                            continue
                        #print "So we use the charge of the other as rough approx (to deal better with loudest sta) : ", charge
                    else:
                        charge = wave.charge
                    charge_map.append((charge,om))

            charge_map = sorted(charge_map, reverse=True)
            if len(charge_map) < 1:
                self.PushFrame(frame,"OutBox")
                return
            q1 = charge_map[0][0]
            q1_dom = charge_map[0][1]
            #print q1, q1_dom
            i = 0
            if q1_dom.om == 61: 
                if icetray.OMKey(q1_dom.string,q1_dom.om+3) in tank_map:
                    i = 3
                elif icetray.OMKey(q1_dom.string,q1_dom.om+2) in tank_map:
                    i = 2
            elif q1_dom.om == 62:
                if icetray.OMKey(q1_dom.string,q1_dom.om+1) in tank_map:
                    i = 1
                elif icetray.OMKey(q1_dom.string,q1_dom.om+2) in tank_map:
                    i = 2
            elif q1_dom.om == 63: 
                if icetray.OMKey(q1_dom.string,q1_dom.om-1) in tank_map:
                    i = -1
                elif icetray.OMKey(q1_dom.string,q1_dom.om-2) in tank_map:
                    i = -2
            elif q1_dom.om == 64:
                if icetray.OMKey(q1_dom.string,q1_dom.om-3) in tank_map:
                    i = -3
                elif icetray.OMKey(q1_dom.string,q1_dom.om-2) in tank_map:
                    i = -2
            if i != 0:
                q2_dom = icetray.OMKey(q1_dom.string,q1_dom.om+i)
                pulses_other = tank_map[q2_dom]
                pulse_other = pulses_other[0]
                q2 = pulse_other.charge
                frame['Q1b'] = dataclasses.I3Double(q2)

            bookedN = 0
            while (bookedN < self.nPulses) and (bookedN < charge_map.__len__()):
                name = "Q%i" % (bookedN+1)
                frame[name] = dataclasses.I3Double(charge_map[bookedN][0])
                bookedN += 1
        self.PushFrame(frame,"OutBox")


################################################################################
# 8. Ratio of hit stations to total nstations in circle around COG (with R to
# furthest tank) #
################################################################################
class StationDensity(icetray.I3Module):
    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        self.AddParameter('InputITpulses','Which IceTop Pulses to use',0)
        self.AddParameter('InputCOG','Which IceTop COG to use',0)
        self.AddOutBox("OutBox")

    def Configure(self):
        self.pulses = self.GetParameter('InputITpulses')
        self.cog = self.GetParameter('InputCOG')

    def Geometry(self, frame):
        if 'I3Geometry' in frame:
            geom = frame['I3Geometry']
            self.stageo = geom.stationgeo
        else:
            print 'No geometry found'
        self.PushFrame(frame,"OutBox")

    def Physics(self, frame):
        if self.cog in frame:
            cog = frame[self.cog]
            if 'I3EventHeader' in frame:
                evhead = frame['I3EventHeader']
                run = evhead.run_id
                event = evhead.event_id
                R_om = []
                if self.pulses in frame:
                    itvem = frame[self.pulses]
                    if itvem.__class__ == dataclasses.I3RecoPulseSeriesMapMask:
                        itvem = itvem.apply(frame)

                    largeDist = 0
                    nstation = 0
                    nTotalSta = 0
                    for sta, stageos in self.stageo:
                        #stageos = [tankA,tankB]
                        dist = (stageos[0].position.calc_distance(cog.pos) +
                                stageos[1].position.calc_distance(cog.pos))/2.
                        R_om.append((dist,sta))
                        # if station is hit => one of the DOMs in tank is hit
                        if icetray.OMKey(sta,61) in itvem or icetray.OMKey(sta,62) in itvem:
                            nstation +=1
                            if dist > largeDist:
                                largeDist = dist

                    #print sorted(R_om)
                    #print nstation, largeDist
                    withinRlist = []
                    for tup in sorted(R_om):
                        if tup[0] <= largeDist+0.1:   #+epsilon with e=0.1
                            withinRlist.append(tup[1])
                    #print withinRlist
                    #print set(withinRlist)
                    nTotalSta = len(set(withinRlist))
                    if nstation <= 0 or nTotalSta <= 0:
                        #print 'Oh no, nstation %i ntotalSta %i, dist %f' % (nstation,nTotalSta,largeDist)
                        staDens = -1.
                    else :
                        staDens = float(nstation)/float(nTotalSta)
                    if staDens > 1.:
                        print nTotalSta
                        print nstation, largeDist
                        print sorted(R_om)
                    #print 'StationDensity is ',staDens
                    frame["StationDensity"] = dataclasses.I3Double(staDens)
        self.PushFrame(frame,"OutBox")


###################################
# 9. The actual cuts and counters
###################################
class ContainmentCuts(icetray.I3Module):
    def __init__(self, context):
        icetray.I3Module.__init__(self,context)
        self.AddParameter('LoudestStationName','Which Loudest Station (OnEdge) To Check for whether its on the edge',0)
        self.AddParameter('CutTail','Keep events with largest charge above this value',0)
        self.AddParameter('KillerCut','Keep events with neighbouring tank of largest charge above this value',0)
        self.AddParameter('DensityCut','Keep events with a stationDensity higher than this value',0)
        self.AddParameter('RemoveEvents','Really remove the events or only count?',0)
        self.AddParameter('NChNStaCutName','Name of NCh + NSta cut in the frame',0)
        self.AddOutBox("OutBox")

    def Configure(self):
        self.loudstaname = self.GetParameter('LoudestStationName')
        self.tailcut = self.GetParameter('CutTail')
        self.killcut = self.GetParameter('KillerCut')
        self.denscut = self.GetParameter('DensityCut')
        self.remove = self.GetParameter('RemoveEvents')
        self.nch_nsta_name = self.GetParameter('NChNStaCutName')
        # a bunch of counters :
        self.nTotal = 0
        self.nNchNsta = 0
        self.nLoudSta = 0
        self.nLoudCha = 0
        self.nKillCut = 0
        self.nDens = 0



    def Physics(self, frame):
        self.nTotal += 1
        #  if 'CoincOfflinePulses_NChNStaCut' in frame:
        if self.nch_nsta_name in frame:
            ncha_nsta = frame[self.nch_nsta_name]
            if ncha_nsta.value:  # If event has MORE than Ncha AND Nsta:
                self.nNchNsta += 1
                if self.loudstaname in frame:
                    loudstaCut = frame[self.loudstaname]
                    if not loudstaCut.value : # if loudest station is NOT on the edge
                        self.nLoudSta += 1
                        if 'Q1' in frame:
                            q1 = frame['Q1']
                            if q1.value > self.tailcut:
                                self.nLoudCha += 1 
                                if 'Q1b' in frame:
                                    q1b = frame['Q1b']
                                    if q1b.value > self.killcut:
                                        self.nKillCut += 1
                                        if 'StationDensity' in frame:
                                            staDens = frame['StationDensity']
                                            if staDens.value > self.denscut:
                                                self.nDens += 1
                                                self.PushFrame(frame,"OutBox")
                                                return
        else :
            print 'Name %s not found in frame ' % (self.nch_nsta_name)
        if not self.remove:
            self.PushFrame(frame,"OutBox")

    def Finish(self):
        print 'Number of events after Coincident Filter : ', self.nTotal
        if self.nTotal > 0:
            print 'Number of events after NStaNCh cut : %i (of total : %lf, of previous : %lf)' % (self.nNchNsta,float(self.nNchNsta)/float(self.nTotal),float(self.nNchNsta)/float(self.nTotal))
            if self.nNchNsta > 0:
                print 'Number of events after LoudestStation cut : %i (of total : %lf, of previous : %lf)' % (self.nLoudSta,float(self.nLoudSta)/float(self.nTotal),float(self.nLoudSta)/float(self.nNchNsta))
                if self.nLoudSta > 0:
                    print 'Number of events after LoudestCharge cut : %i (of total : %lf, of previous : %lf)' % (self.nLoudCha,float(self.nLoudCha)/float(self.nTotal),float(self.nLoudCha)/float(self.nLoudSta))
                    if self.nLoudCha > 0:
                        print 'Number of events after Killercut cut : %i (of total : %lf, of previous : %lf)' % (self.nKillCut,float(self.nKillCut)/float(self.nTotal),float(self.nKillCut)/float(self.nLoudCha))
                        if self.nKillCut > 0:
                            print 'Number of events after StationDensity cut : %i (of total : %lf, of previous : %lf)' % (self.nDens,float(self.nDens)/float(self.nTotal),float(self.nDens)/float(self.nKillCut))


class QualityCuts(icetray.I3Module):
    def __init__(self, context):
        icetray.I3Module.__init__(self,context)
        self.AddParameter('RecoName','Name of the IceTop reconstruction',0)
        self.AddParameter('CutMinBeta','Keep events with beta higher than this value',0)
        self.AddParameter('CutMaxBeta','Keep events with beta lower than this value',0)
        self.AddParameter('CutCornerClippers','Keep events with InIceSize above this value',0)
        self.AddParameter('RemoveEvents','Really remove the events or only count?',0)
        self.AddOutBox("OutBox")

    def Configure(self):
        self.reconame = self.GetParameter('RecoName')
        self.maxbetacut = self.GetParameter('CutMaxBeta')
        self.minbetacut = self.GetParameter('CutMinBeta')
        self.iiSizeCut = self.GetParameter('CutCornerClippers')
        self.remove = self.GetParameter('RemoveEvents')
        # a bunch of counters :
        self.nTotal = 0
        self.nStatus = 0
        self.nBeta = 0
        self.nIIsize = 0
   
    def Geometry(self, frame):
        if 'I3Geometry' in frame:
            geom = frame['I3Geometry']
            self.omgeo = geom.omgeo
            # Calculate InIceSiz :
            # CONTAINMENT PART FROM FLAT-NTUPLE : doing this part in the
            # Geometry method
            # saves a lot of CPUtime
            strings = icetray.vector_double()
            #IC79
            strings.append(2);
            strings.append(6);
            strings.append(50);
            strings.append(74);
            strings.append(73);
            strings.append(78);
            strings.append(75);
            strings.append(41);
            self.x = icetray.vector_double()
            self.y = icetray.vector_double()
            self.avgtop = 0
            self.avgbottom = 0
            for s in strings:
                omk = icetray.OMKey(int(s),30)  #middle of string on the edge
                omk_top = icetray.OMKey(int(s),1)  #top of string on the edge
                omk_bottom = icetray.OMKey(int(s),60)  #bottom of string on the edge
                for om,omg in self.omgeo:
                    if om == omk:
                        self.x.append(omg.position.x)
                        self.y.append(omg.position.y)
                    if om == omk_top:
                        self.avgtop += omg.position.z
                    if om == omk_bottom:
                        self.avgbottom += omg.position.z

            self.avgtop /= len(strings)
            self.avgbottom /= len(strings)
            #Now vectors are filled

        else:
            print 'No geometry found'
        self.PushFrame(frame,"OutBox")


    def Physics(self, frame):
        self.nTotal += 1
        if self.reconame in frame:
            laputop = frame[self.reconame]
            if laputop.fit_status == dataclasses.I3Particle.OK:
                self.nStatus += 1
                params = self.reconame + 'Params'
                if params in frame:
                    lapu_params = frame[params]
                    if lapu_params.beta > self.minbetacut and lapu_params.beta < self.maxbetacut:
                        self.nBeta += 1
                        IIcont = phys_services.I3Cuts.containment_volume_size(laputop, self.x, self.y, self.avgtop, self.avgbottom)
                        if IIcont < self.iiSizeCut:
                            self.nIIsize += 1
                            self.PushFrame(frame,"OutBox")
                            return
            else :
                print 'Fit is not OK but : ', laputop.fit_status
        if not self.remove:
            self.PushFrame(frame,"OutBox")

    def Finish(self):
        print 'Number of events at beginning of QualityCuts : ', self.nTotal
        if self.nTotal > 0:
            print 'Number of events after Status cut : %i (of total : %lf, of previous : %lf)' % (self.nStatus,float(self.nStatus)/float(self.nTotal),float(self.nStatus)/float(self.nTotal))
            if self.nStatus > 0:
                print 'Number of events after Beta cut : %i (of total : %lf, of previous : %lf)' % (self.nBeta,float(self.nBeta)/float(self.nTotal),float(self.nBeta)/float(self.nStatus))
                if self.nBeta > 0:
                    print 'Number of events after IIsize cut : %i (of total : %lf, of previous : %lf)' % (self.nIIsize,float(self.nIIsize)/float(self.nTotal),float(self.nIIsize)/float(self.nBeta))


####################################
## Add Stations from these Pulses to an existing (or not yet) 
## ExcludedStationList
## Application : SRTexcludedPulses to ClusterCleanedExcludedStations...
####################################
class AddStations(icetray.I3Module):
    def __init__(self,context):
        icetray.I3Module.__init__(self, context)
        self.AddParameter('RecoPulsesName','Name of recoPulses in the frame',0)
        self.AddParameter('ExcludedStationsName','Name of Already ExcludedStations in the frame',0)
        self.AddParameter('OutputExcludedStationsName','Name of New ExcludedStations in the frame',0)
        self.AddOutBox('OutBox')

    def Configure(self):
        self.pulsename = self.GetParameter('RecoPulsesName')
        self.exclstationname = self.GetParameter('ExcludedStationsName')
        self.outputname = self.GetParameter('OutputExcludedStationsName')

    def Physics(self,frame):
        exclstations = 0
        if not self.exclstationname in frame:
            exclstations = dataclasses.I3VectorInt()
        else:
            exclstations = dataclasses.I3VectorInt(frame[self.exclstationname]) # Do a COPY !! 

        if self.pulsename in frame:
            vemPulses = frame[self.pulsename]
            if vemPulses.__class__ == dataclasses.I3RecoPulseSeriesMapMask:
                vemPulses = vemPulses.apply(frame)
            stationList = set([pulse.key().string for pulse in vemPulses])
            for sta in stationList:
                exclstations.append(sta)
        frame[self.outputname] = exclstations
        self.PushFrame(frame,"OutBox")


# Bottom 2 need some slight v4 adaptation, but are not used in L3, could be
# useful for other studies though
###############################################################################
# 10. LoudestNonFluctuating station : also books loudest N pulses, as doubles
# (less ROOT c++ work)##
################################################################################
class FindLoudestNonFluctStation(icetray.I3Module):
    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        self.AddParameter('nPulses','Book the largest N pulses for TailCut',4)
        self.AddParameter('ITpulses','IT pulses Name',0)
        self.AddOutBox("OutBox")

    def Configure(self):
        self.nPulses = self.GetParameter('nPulses')
        self.VEM = self.GetParameter('ITpulses')
        self.counter=0

    def Physics(self, frame):
        if frame.Has(self.VEM):
            VEM = frame.Get(self.VEM)   # This is a map<OMKey,RecoPulseSeries>, sorting by OMKey, so first strings, then OMs
            pulselist = []
            for om,pulseseries in VEM:
                for wave in pulseseries:
                    charge = -1
                    #if math.isnan(wave.charge): #only works beyond python 2.6
                    if wave.charge != wave.charge:
                        print "Charge is NaN, so saturated HG without LG (or DOM 39-62)" # Other tank will be next or previous element in map, use charge of other tank as rough approx
                        omkey = icetray.OMKey()
                        omkey.SetString(om.string)
                        print "String ", om.string, " OM ", om.om
                        # Way too complicated way to get the omkey of the
                        # neighbouring tank
                        if om.om == 61 or om.om == 62:
                            omkey.SetOM(om.om+2)
                            if not (VEM.has_key(omkey)):
                                if om.om == 61:
                                    omkey.SetOM(om.om+3)
                                else:
                                    omkey.SetOM(om.om+1)
                        elif om.om == 63 or om.om == 64:
                            omkey.SetOM(om.om-2)
                            if not (VEM.has_key(omkey)):
                                if om.om == 63:
                                    omkey.SetOM(om.om-1)
                                else:
                                    omkey.SetOM(om.om-3)
                        print omkey
                        if omkey in VEM:
                            pulses = VEM[omkey]
                            charge = pulses[0].charge
                        else:
                            continue
                        print "So we use the charge of the other as rough approx (to deal better with loudest sta) : ", charge
                    else:
                        charge = wave.charge
                    pulselist.append((charge,om.string,om.om))

            newList = sorted(pulselist, reverse=True)  # if two charges are the same, sort keeps the present order of the list, amazing.

            ################################################################
            ## New From TailStudy : Book 6 largest pulses (and cut on those)
            #################################################################
            # newList is a [ ( , , ), ( , , ), ... ] list of
            # (charge,string/station,DOM) tuples
            # and it's already ordened by charge!
            # BOOK first N
            bookedN = 0
            while (bookedN < self.nPulses) and (bookedN < newList.__len__()):
                name = "Q%i" % (bookedN+1)
                frame[name] = dataclasses.I3Double(newList[bookedN][0])
                bookedN += 1

            # DEBUGGING : print newList
            StationNotFound = True
            attempt = 0
            loudSta = -1
            ### NEW ###
            # Make distributions of distance of loudest pulse to neighbouring
            # tank in Qspace (charge), so keep this distance
            distanceLoudest = 0
            ###########
            while StationNotFound and (len(newList) > 1):
                attempt += 1
                loudSta = newList[0][1]  # newList[0], tuple with loudest Charge, [1] is the StationNumber
                print "Tank one is ", loudSta, ", Looking how far away (in the charge ordered list) tank two is"
                del(newList[0])       # remove loudest expected tank from list to make looping easier to look for neighbour

                ## Next : Loop over rest of newList, looking for the other Tank (thus same Sta), calculating the distance
                dist = 0
                for tuple in newList:
                    dist += 1
                    if tuple[1] == loudSta:  # must be in there somewhere (because cleaned in TEB, no lonely tanks!)
                        # we found the neighbouring tank
                        # dist == 1 : normal, standard
                        # dist == 2 : can happen sometimes, don't be too strict
                        # dist > 2 : neighbouring tank is tank is far away, so
                        # loudest tank was fluctuating upwards,
                        # cut away neighbouring tank and repeat procedure (max
                        # twice).
                        print "Bingo, other tank found, how far? ", dist
                        if dist >= 3:
                            del(newList[dist-1]) # or do it outside the loop??
                            print "This was a fluctuating station, repeat procedure"
                        else:
                            StationNotFound = False

                        if attempt == 1:
                            distanceLoudest = dist
                        #print "Dist is ", dist, " distance Loudest = ", distanceLoudest
                    #break
                print "This is attempt ", attempt 
                if attempt > 2:
                    StationNotFound = False
            #Put loudestStation in Frame as Double
            eventheader = frame["I3EventHeader"]
            #print "Event : ", eventheader.EventID
            frame["LoudestNonFluctSta"] = dataclasses.I3Double(loudSta)

            ### NEW ###
            # Put distance of loudestTank to its neighbour in the frame
            frame["QDistLoudest"] = dataclasses.I3Double(distanceLoudest)
            ###########
        else:
            print "WARN: IT pulses %s not found" % (self.VEM)

        self.PushFrame(frame,"OutBox")


###############################
## 11. Could be pretty useful
###############################
class LoudestStationToMC(icetray.I3Module):

    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        self.AddParameter('InputLoudestStation','Which Loudest Station',0)
        self.AddParameter('OutputName','What is the name of the output I3Double','DistanceToMC')
        self.AddOutBox("OutBox")

    def Configure(self):
        self.loudest = self.GetParameter('InputLoudestStation')
        self.outputName = self.GetParameter('OutputName')

    def Physics(self, frame):
        if frame.Has('MCPrimary') and frame.Has('I3Geometry') and frame.Has(self.loudest):
            mcprim = frame.Get('MCPrimary')
            loud = frame.Get(self.loudest)
            geo = frame.Get('I3Geometry')
            x = 0
            y = 0
            for omkey, g in geo.omgeo:
                pos =g.position
                if omkey.string == loud.value:
                    if omkey.om > 60:
                        #print pos.X ,",", pos.Y
                        x += pos.X
                        y += pos.Y
            position = dataclasses.I3Position(x/4,y/4,mcprim.GetZ())
            distance = position.CalcDistance(mcprim.GetPos())
            #        print distance
            #        for station in
            #        [21,30,49,50,59,67,74,78,77,76,75,68,60,52,53,44,45,46,38,39]:
            #        # the IC40 borderlist
            #            if loud == station:
            #                print "We're on the edge of the array"
            if x == 0 and y == 0:
                distance = -1
            output = dataclasses.I3Double(distance)
            print "You Chose ", self.outputName," as name for the output I3Double"
            # frame["LoudestStationToMC"] = output
            frame[self.outputName] = output
        else :
            print "WARN : Either no MCPrimary - no I3Geometry or no LoudestStation Called %s" % self.loudest
        self.PushFrame(frame,"OutBox")


