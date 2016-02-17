#!/usr/bin/env python

import numpy as np
import healpy as hp
import tables, glob

import myGlobals as my
from showerllh.run_data.save_hist import hdf5extractor

def histSpectrum(config):

    files = glob.glob('%s/%s_data/files/*.hdf5' % (my.llh_data, config))
    files.sort()

    nside = 64
    npix = hp.nside2npix(nside)
    sbins = np.arange(npix, dtype=int)
    ebins = np.arange(5, 9.501, 0.05)

    for file in files[:10]:

        print 'Working on %s...' % file

        d = hdf5extractor(config, file)
        c0 = d['cuts']['llh']
        r = np.log10(d['ML_energy'])[c0]
        fit = zfix(d['zenith'], bintype='logdist')[c0]
        w = d['weights'][c0]

        zen = d['zenith'][c0]
        azi = d['azimuth'][c0]
        pix = hp.ang2pix(nside, zen, azi)

        x = r - fit
        y = pix
        h0 = np.histogram2d(x, y, bins=(ebins,sbins), weights=w)[0]


if __name__ == "__main__":

    histSpectrum('IT81')
