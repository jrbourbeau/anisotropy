#!/usr/bin/env python

import re, glob

def getEbins():
    return ['4','4.25','4.5','4.75','5','5.25','5.5','6','6.5','100']

def getEnergyMaps(config, ebins=[]):

    if ebins == []:
        ebins = getEbins()

    prefix = '/data/user/fmcnally/anisotropy/maps/merged'
    mapFiles = glob.glob('%s/%s_24H_sid_*GeV.fits' % (prefix, config))
    mapFiles.sort()
    emins = [float(re.split('_|-|GeV', f)[-3]) for f in mapFiles]
    emaxs = [float(re.split('_|-|GeV', f)[-2]) for f in mapFiles]

    fileList = []
    for i in range(len(ebins)-1):
        fileList += [[file for j, file in enumerate(mapFiles) \
                if emins[j]>=float(ebins[i]) and emaxs[j]<=float(ebins[i+1])]]

    return fileList
