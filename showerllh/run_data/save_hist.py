#!/usr/bin/env python

###########################################################################
# Replacement for save_data.py takes a collection of hdf5 files, and      #
# builds desired histograms for rapid plotting                            #
###########################################################################

import numpy as np
import healpy as hp
import argparse, tables

from icecube import astro

import dataFunctions as df
from showerllh.analysis.skypos import getDecRA
from showerllh.analysis.llhtools import inPoly, getEbins
from showerllh.analysis.zfix import zfix


def hdf5extractor(config, file):

    ##=======================================================================
    ## Starting parameters

    rDict = {'proton':'p','helium':'h','oxygen':'o','iron':'f'}
    t1 = astro.Time()

    print 'Building arrays from %s...' % file
    t = tables.openFile(file)
    q = {}

    # Get reconstructed compositions from list of children in file
    children = []
    for node in t.walk_nodes('/'):
        try: children += [node.name]
        except tables.NoSuchNodeError:
            continue
    children = list(set(children))
    compList = [n.split('_')[-1] for n in children if 'ShowerLLH_' in n]

    # Get ShowerLLH cuts and info
    rrc = t.root.ShowerLLH_proton.col('exists').astype('bool')
    for value in ['zenith', 'azimuth']:
        q[value] = t.root.ShowerLLH_proton.col(value)
    for comp in compList:
        r = rDict[comp]
        for value in ['x','y','energy']:
            q[r+'ML_'+value] = t.getNode('/ShowerLLH_'+comp).col(value)
        q[r+'LLH'] = t.getNode('/ShowerLLHParams_'+comp).col('maxLLH')

    # Timing
    mjd_day = t.root.I3EventHeader.col('time_start_mjd_day')
    mjd_sec = t.root.I3EventHeader.col('time_start_mjd_sec')
    mjd_ns  = t.root.I3EventHeader.col('time_start_mjd_ns')
    q['mjd'] = np.zeros(len(mjd_day), dtype=np.float64)
    for i in range(len(mjd_day)):
        day = int(mjd_day[i])
        sec = int(mjd_sec[i])
        ns  = int(mjd_ns[i])
        t1.SetTime(day, sec, ns)
        q['mjd'][i] = t1.GetMJD()

    # Event ID
    run = t.root.I3EventHeader.col('Run')
    event = t.root.I3EventHeader.col('Event')
    subevent = t.root.I3EventHeader.col('SubEvent')
    eventIDs = []
    for i in range(len(run)):
        eventIDs += ['%s_%s_%s' % (run[i], event[i], subevent[i])]
    q['eventIDs'] = np.asarray(eventIDs)

    # Condition and prescale passed (filter[condition, prescale])
    # For notes on weights see bottom of file
    filtermask  = df.filter_mask(config)
    filternames = df.filter_names(config)
    f = {}
    for fname in filternames:
        f[fname] = t.getNode('/'+filtermask).col(fname)
        f[fname] = f[fname][:,0].astype(float)
    filterArray = np.array([f[fname] * df.it_weights(fname) 
            for fname in f.keys()])
    filterArray[filterArray == 0] = 100.
    q['weights'] = np.amin(filterArray, axis=0)

    # Other reconstruction info
    q['NStations'] = t.root.NStations.col('value')

    t.close()

    # Laputop values
    #for key in ['x','y','zenith','azimuth','s125','e_proton','e_iron','beta']:
    #    arrays += ['lap_'+key]

    # Get Laputop info
    #for value in ['x', 'y', 'zenith', 'azimuth']:
    #    q['lap_'+value] = t.root.Laputop.col(value)
    #for value in ['s125', 'e_proton', 'e_iron', 'beta']:
    #    q['lap_'+value] = t.root.LaputopParams.col(value)

    # Get most likely composition
    rList = [rDict[comp] for comp in compList]
    full_llhs = np.array([q[r+'LLH'] for r in rList])
    max_llh = np.amax(full_llhs, axis=0)
    q['llh_comp'] = np.array(['' for i in range(len(q['pLLH']))])
    for r in rList:
        q['llh_comp'][q[r+'LLH'] == max_llh] = r
    for key in ['x', 'y', 'energy']:
        q['ML_'+key] = np.array([q[r+'ML_'+key][i]
                for i, r in enumerate(q['llh_comp'])])

    # Check for multiple most-likely compositions (mark as bad)
    badVals = np.sum(full_llhs == max_llh, axis=0)
    badVals = (badVals-1).astype('bool')
    q['llh_comp'][badVals] = ''
    for key in ['x','y','energy']:
        q['ML_'+key][badVals] = np.nan

    # Calculate sky positions
    q['dec'], q['ra'] = getDecRA(q, verbose=False)

    # Containment cut
    it_geo = df.it_geo(config)
    q['cuts'] = {}
    q['cuts']['llh'] = inPoly(q['ML_x'], q['ML_y'], 0, config=it_geo)

    return q


def histWriter(config, file, outfile):

    # Bin values
    eList = ['p','h','o','f']
    decbins = ['0-12','12-24','24-40']
    rabins = ['0-60','60-120','120-180','180-240','240-300','300-360']

    # Build general list of key names to write
    keyList = []
    keyList += ['energy','energy_w','energy_z','energy_w_z']
    keyList += ['zenith','zenith_w','core','core_w']
    keyList += ['%s_err' % k for k in keyList]
    # Split by composition
    keyList = ['%s_%s' % (k, e) for k in keyList for e in eList]
    # Split spatially into 3 dec bins (~12 degrees each) and 6 ra bins
    keyList = ['%s_%s_%s' % (k, dec, ra) for k in keyList \
            for dec in decbins for ra in rabins]

    # Extract information from hdf5 file
    q = hdf5extractor(config, file)
    c0 = q['cuts']['llh']
    r = np.log10(q['ML_energy'])[c0]
    cosz = np.cos(q['zenith'])[c0]
    dist = np.sqrt(q['ML_x']**2 + q['ML_y']**2)[c0]
    fit = zfix(q['zenith'], bintype='logdist')[c0]
    w = q['weights'][c0]

    # Make cuts
    degree = np.pi / 180.
    for e in eList:
        q[e] = (q['llh_comp'] == e)[c0]
    for dec in decbins:
        decmin = (180 - float(dec.split('-')[1])) * degree
        decmax = (180 - float(dec.split('-')[0])) * degree
        q[dec] = ((q['dec'] >= decmin) * (q['dec'] < decmax))[c0]
    for ra in rabins:
        ramin = float(ra.split('-')[0]) * degree
        ramax = float(ra.split('-')[1]) * degree
        q[ra] = ((q['ra'] >= ramin) * (q['ra'] < ramax))[c0]

    # Method of intelligently producing histograms based on key names
    def smartHist(key, x, bins):
        tempx = x
        wts = None
        params = key.split('_')
        e, dec, ra = params[-3:]
        c1 = q[e] * q[dec] * q[ra]
        if 'z' in params:
            tempx = x - fit
        if 'w' in params:
            wts = w[c1]
            if 'err' in params:
                wts = (w[c1])**2
        h0 = np.histogram(tempx[c1], bins=bins, weights=wts)[0]
        return h0

    # Energy distribution
    h = {}
    print 'Calculating energy distributions...'
    bins = getEbins(reco=True)
    energyKeys = [k for k in keyList if 'energy' in k]
    for key in energyKeys:
        h[key] = smartHist(key, r, bins)

    # Zenith distribution
    print 'Calculating zenith distributions...'
    bins = np.linspace(0.8, 1, 81)
    zenithKeys = [k for k in keyList if 'zenith' in k]
    for key in zenithKeys:
        h[key] = smartHist(key, cosz, bins)

    # Core distribution
    print 'Calculating core position distributions...'
    bins = np.linspace(0, 700, 141)
    coreKeys = [k for k in keyList if 'core' in k]
    for key in coreKeys:
        h[key] = smartHist(key, dist, bins)

    print 'Saving...'
    np.save(outfile, h)


def skyWriter(config, file, outfile):

    nside = 64
    npix = hp.nside2npix(nside)

    # Binning for various parameters
    sbins = np.arange(npix+1, dtype=int)
    ebins = np.arange(5, 9.501, 0.05)
    dbins = np.linspace(0, 700, 141)
    lbins = np.linspace(-20, 20, 151)

    # Get desired information from hdf5 file
    d = hdf5extractor(config, file)
    c0 = d['cuts']['llh']
    r = np.log10(d['ML_energy'])[c0]
    fit = zfix(d['zenith'], bintype='logdist')[c0]
    w = d['weights'][c0]
    xy = np.sqrt(d['ML_x']**2 + d['ML_y']**2)[c0]
    dllh = (d['fLLH'] - d['pLLH'])[c0]

    # Bin in sky
    #zen = np.pi - d['zenith'][c0]
    #azi = d['azimuth'][c0]
    dec = d['dec'][c0]
    ra  = d['ra'][c0]
    x = hp.ang2pix(nside, dec, ra)

    # Energy cut
    ecut = (r >= 6.2)

    p = {'weights':w}
    q = {}
    q['energy_w'] = np.histogram2d(x, r-fit, bins=(sbins,ebins), **p)[0]
    q['dist_w'] = np.histogram2d(x, xy, bins=(sbins,dbins), **p)[0]
    q['llh_w'] = np.histogram2d(x, dllh, bins=(sbins,lbins), **p)[0]
    p['weights'] = w**2
    q['energy_err_w'] = np.histogram2d(x, r-fit, bins=(sbins,ebins), **p)[0]
    q['dist_err_w'] = np.histogram2d(x, xy, bins=(sbins,dbins), **p)[0]
    q['llh_err_w'] = np.histogram2d(x, dllh, bins=(sbins,lbins), **p)[0]
    # Energy cut versions
    q['llhcut_w'] = np.histogram2d(x[ecut], dllh[ecut], bins=(sbins,lbins),
            weights=w[ecut])[0]
    q['llhcut_err_w'] = np.histogram2d(x[ecut], dllh[ecut], bins=(sbins,lbins),
            weights=(w[ecut])**2)[0]

    q['energy'] = np.histogram2d(x, r-fit, bins=(sbins,ebins))[0]
    q['dist'] = np.histogram2d(x, xy, bins=(sbins,dbins))[0]
    q['llh'] = np.histogram2d(x, dllh, bins=(sbins,lbins))[0]
    q['energy_err'] = np.histogram2d(x, r-fit, bins=(sbins,ebins))[0]
    q['dist_err'] = np.histogram2d(x, xy, bins=(sbins,dbins))[0]
    q['llh_err'] = np.histogram2d(x, dllh, bins=(sbins,lbins))[0]
    # Energy cut versions
    q['llhcut'] = np.histogram2d(x[ecut], dllh[ecut], bins=(sbins,lbins))[0]
    q['llhcut_err'] = np.histogram2d(x[ecut], dllh[ecut], bins=(sbins,lbins))[0]

    np.save(outfile, q)



if __name__ == "__main__":

    p = argparse.ArgumentParser(description='Converts hdf5 file to npy dict')
    p.add_argument('-c', '--config', dest='config',
            help='Detector configuration [IT73 --> IT81-IV]')
    p.add_argument('-f', '--files', dest='files', nargs='*',
            help='Input files')
    p.add_argument('-o', '--outfiles', dest='outfiles', nargs='*',
            help='Output files')
    p.add_argument('--sky', dest='sky',
            default=False, action='store_true',
            help='Write the sky histograms')
    args = p.parse_args()

    for infile, outfile in zip(args.files, args.outfiles):
        if args.sky:
            skyWriter(args.config, infile, outfile)
        else:
            histWriter(args.config, infile, outfile)


###############################################################################
## Notes on weights

"""
 - Events that pass STA8 condition have a prescale and weight of 1.
 - Events that pass STA3ii condition have a 1/2 chance to pass the STA3ii
prescale. Those that fail have a 1/3 chance to pass the STA3 prescale. So, the 
total chance of an event passing is 1/2+(1/3*1/2) = 2/3. Weight = 1.5
 - Events that pass STA3 condition but not STA3ii condition have a prescale and
weight of 3
"""
