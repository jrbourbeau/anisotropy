#!/usr/bin/env python

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import glob, os, tabulate, subprocess
import mapFunctions
from mapfunctions.energyCuts import getEbins, getEnergyMaps

def plot(verbose=False):

    prefix = '/data/user/fmcnally/anisotropy/maps/merged'
    alpha = 1/20.
    deg = 180./np.pi
    colorDict = {0:'b'}

    # IceCube information
    files, meds, sigL, sigR = [],[],[],[]
    ebins = getEbins()
    files = getEnergyMaps('IC')
    print('ebins = {}'.format(ebins))
    print('files = {}'.format(files))
    dipoleParams = np.zeros((len(files), 6))
    for i, files_i in enumerate(files):
        #print('files_{}[0] = {}'.format(i,files_i[0]))
        data, bg, local = hp.read_map(files_i[0], range(3))
    opts = {'polar': False, 'swindow': 3, 'verbose': True, 'customOut': None,\
            'gplane': False, 'max': '1.0', 'batch': False, 'ramax': None,\
            'half': True, 'threshold': None, 'decmin': -90.0, 'mapName': 'fit',\
            'multi': None, 'scale': 3.0, 'output': False, 'title': False,\
            'ext': 'png', \
            'filepath': '/data/user/fmcnally/anisotropy/maps/merged/IC_24H_sid_4.25-4.5GeV.fits', \
            'smooth': 5.0, 'mask': False, 'llabel': False, 'rlabel': False,\
            'coords': 'C', 'min': '-1.0', 'decmax': -25.0, 'ramin': None,\
            'prelim': False, 'fix_multi': False, 'stype': 'tophat',\
            'outDir': '/home/jbourbeau/public_html/figures/'}
    #map = (data-bg)/bg
    print('filepath = {}'.format(opts['filepath']))
    map = mapFunctions.getMap(*[opts['filepath']],**opts)
    # Multiply by scale
    if opts['scale']:
        map[map!=hp.UNSEEN] *= (10**opts['scale'])

    # Setup the coordinate system
    rot = 180
    
    # Mollweide (default) or Cartesian projection
    proj = 'Mollweide'
    projaxis = hp.projaxes.HpxMollweideAxes
    
    # Find min and max for unmasked pixels
    unmasked = np.array([i for i in map if (i!=hp.UNSEEN and i!=np.inf)])
    min = float(opts['min']) if opts['min'] else unmasked.min()
    max = float(opts['max']) if opts['max'] else unmasked.max()
    
    colormap = plt.get_cmap('jet')
    
    pltParams = {'fig':1, 'rot':rot, 'min':min, 'max':max, \
            'cbar':False, 'notext':True, 'coord':opts['coords'], \
            'cmap':colormap}

    # Eliminate nan's and mask
    map[np.isnan(map)] = 0
    map = mapFunctions.maskMap(map, -90, 90)
    hp.mollview(map,**pltParams)
    plt.show()

if __name__ == "__main__":
    #plot(verbose=True)
    #subprocess.check_output(["echo", "Hello World!"])
    #subprocess.call('./plotFITS.py --half --mask -n fit -v --scale 3 --min -1.0 --max 1.0 /data/user/fmcnally/anisotropy/maps/merged/IC_24H_sid_4.25-4.5GeV.fits')
    subprocess.call(['./plotFITS.py','--half','--mask', '-v', '-b', '-nfit','-s 3',\
            '-m -1.0', '-M 1.0', '-f /data/user/fmcnally/anisotropy/maps/merged/IC_24H_sid_4.25-4.5GeV.fits'])
    files, meds, sigL, sigR = [],[],[],[]
    ebins = getEbins()
    files = getEnergyMaps('IC')
    print('ebins = {}'.format(ebins))
    print('files = {}'.format(files))
