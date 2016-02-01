#!/usr/bin/env python
from icecube import icetray, dataclasses
import csv, time, sys
csv.register_dialect('spaces',delimiter=' ',skipinitialspace=True)

class fixsnow(icetray.I3Module):

    def __init__(self,context):

        icetray.I3Module.__init__(self, context)
        self.AddParameter('infile','CSV file with new snow depths',
                '/home/fmcnally/ShowerLLH/useful/snow_height.txt')
        self.AddParameter('runnumber','current run number', 0)
        self.AddOutBox("OutBox")

    def Configure(self):

        self.infile = self.GetParameter('infile')
        self.runnumber = self.GetParameter('runnumber')
        # If we don't have snow for this run, how far are we willing to go
        self.runthreshold = 1e3
        self.allsnow = []
        with open(self.infile) as snowfile:
            cols = snowfile.next().split()
            rows = csv.DictReader(snowfile, fieldnames=cols, dialect='spaces')
            for row in rows:
                self.allsnow.append(row)

    def Geometry(self,frame):

        if self.runnumber:
            self.UpdateSnows(frame)
        else:
            # print("No run number yet, we'll push this frame for now, and update the snow when a physics frame comes around.")
            self.PushFrame(frame)

    def UpdateSnows(self,frame):

        havesnow = False
        runsoutside = self.runthreshold

        for row in self.allsnow:

            first_run = int(row['first_run'])
            last_run  = int(row['last_run'])

            if self.runnumber >= first_run and self.runnumber <= last_run:
                self.newsnow = row
                havesnow = True

            if havesnow == False and min(abs(self.runnumber-first_run),abs(self.runnumber-last_run)) < runsoutside:
                runsoutside = min(abs(self.runnumber-first_run),abs(self.runnumber-last_run))
                self.newsnow = row

        if havesnow == False:
            print("No valid snow for run number: {0}".format(self.runnumber))
            if runsoutside < self.runthreshold:
                print("Using snow for runs: {0}-{1}".format(self.newsnow['first_run'],self.newsnow['last_run']))
            else:
                sys.exit("Could not find snow for any runs within {0} runs.".format(self.runthreshold))

        newgeo = icetray.I3Frame(icetray.I3Frame.Geometry)
        newgeo['I3Geometry'] = frame['I3Geometry']
        stgeo = newgeo['I3Geometry'].stationgeo
        for station,tanks in stgeo:
            tanks[0].snowheight = max(float(self.newsnow['{0:02d}A'.format(station)]),0.0)
            tanks[1].snowheight = max(float(self.newsnow['{0:02d}B'.format(station)]),0.0)
            stgeo[station]=tanks

        #frame['I3Geometry'].stationgeo = stgeo
        self.PushFrame(newgeo)
        if frame.Has('I3DetectorStatus'):
            newds = icetray.I3Frame(icetray.I3Frame.DetectorStatus)
            newds['I3DetectorStatus'] = frame['I3DetectorStatus']
            self.PushFrame(newds)

    def Physics(self,frame):
        if frame['I3EventHeader'].run_id != self.runnumber:
            self.runnumber = frame['I3EventHeader'].run_id
            self.UpdateSnows(frame)
        self.PushFrame(frame)



class printsnow(icetray.I3Module):

    def __init__(self, context):
        icetray.I3Module.__init__(self,context)
        self.AddOutBox("OutBox")

    def Configure(self):
        pass

    def Geometry(self, frame):
        station = 30
        print(station, frame['I3Geometry'].stationgeo[station][0].snowheight, frame['I3Geometry'].stationgeo[station][1].snowheight)
        self.PushFrame(frame)



if __name__ == '__main__':

    import sys

    if len(sys.argv) != 3:
        sys.exit("Usage: {0} gcdfile.i3.gz [runnumber] or [data.i3.gz]".format(sys.argv[0]))
 
    from I3Tray import *
    from icecube import icetray, dataclasses, dataio

    gcd = sys.argv[1]
    if sys.argv[2].find('.i3') != -1:
        files = [gcd, sys.argv[2]]
        runnumber = 0
    else:
        files = [gcd,]
        runnumber = int(sys.argv[2])

    tray = I3Tray()
    tray.AddModule("I3Reader","reader", FilenameList = files)
    tray.AddModule(printsnow,'ps1')
    tray.AddModule(fixsnow,'fs',runnumber = runnumber)
    tray.AddModule(printsnow,'ps2')
    tray.AddModule("TrashCan", "tc")
    tray.Execute()
    tray.Finish()
    del(tray)
