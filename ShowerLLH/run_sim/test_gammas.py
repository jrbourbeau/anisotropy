#!/usr/bin/env python 

################################################################################
# Runs the grid search on the given files and writes to hdf5 output
################################################################################
from icecube import icetray, dataio, MuonGun, dataclasses
from I3Tray import *
from icecube.phys_services import *
from icecube.tableio import I3TableWriter
from icecube.hdfwriter import I3HDFTableService
import numpy as np
import sys, time, glob

def maker(config, outFile, fileList):

    # Keys to write to frame
    keys  = []
    keys += ['I3EventHeader', 'MCPrimary', 'muon_energy', 'I3Triggers']

    t0 = time.time()

    tray = I3Tray()
    tray.AddModule('I3Reader', FileNameList=fileList)
    hdf = I3HDFTableService(outFile)

    tray.AddModule("I3NullSplitter", "fullevent")

    def icecube_muon_energies(frame):
        energies = []
        surface = MuonGun.Cylinder(1000,500)
        for p in MuonGun.muons_at_surface(frame, surface):
            energies.append(p.energy)

        frame['muon_energy'] =  dataclasses.I3VectorDouble(energies)
    '''
    # Write muons to a I3ParticleVect
    def get_muons(frame):
        tree      = frame['I3MCTree']
        particles = tree.get_daughters(frame['I3MCTree'][0])
        muons     = []

        for p in particles:
            if p.type_string == 'MuMinus' or p.type_string == 'MuPlus':
                muons.append(p)

        frame['muons'] = dataclasses.I3ParticleVect(muons) 
    '''
    tray.AddModule(icecube_muon_energies, 'muonEnergies')
    #tray.AddModule(get_muons, 'muons')
    tray.AddModule(I3TableWriter, tableservice=hdf, keys=keys, SubEventStreams = ['fullevent'])

    tray.Execute()
    tray.Finish()

    print "Time taken: ", time.time() - t0



if __name__ == "__main__":

    config   = sys.argv[1]
    outFile  = sys.argv[2]
    fileList = sys.argv[3:]
    maker(config, outFile, fileList)







