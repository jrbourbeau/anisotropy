#!/usr/bin/env python

#=============================================================================
# File Name     : AnisotropyMaps.py
# Description   : Create a class with anisotropy analysis method functions
# Creation Date : 02-19-2016
# Last Modified : Sat 20 Feb 2016 09:27:12 PM CST
# Created By    : James Bourbeau
#=============================================================================

import numpy as np
import healpy as hp

from mapFunctions import getMap
from projFunctions import getHarmonicFitParams


class AnisotropyMaps(object):

    def __init__(self, infile, **kwargs):

        # Apply default values if arguments not in input dictionary
        defaults = {'multi': False, 'smooth': 5,
                    'stype': 'tophat', 'swindow': 3,
                    'verbose': False, 'mask': False, 'decmin': -90.,
                    'decmax': 90., 'ramin':0., 'ramax':360.,
                    'nRAbins':24, 'fix_multi': False}
        opts = {k: kwargs[k] for k in kwargs if k in defaults}
        opts.update({k: defaults[k] for k in defaults if k not in opts})
        self.Opts = opts
        self.MapFile = infile
        self.RelintMap = getMap(*[infile], mapName='relint', **self.Opts)
        self.RelerrMap = getMap(*[infile], mapName='relerr', **self.Opts)
        self.NPix = len(self.RelintMap)
        self.NSide = hp.npix2nside(self.NPix)

    def getNPix(self):
        return self.NPix

    def getRIRAProj(self):
        opts = self.Opts
        relint = self.RelintMap
        relerr = self.RelerrMap

        # Cut to desired dec range (equiv to healpy zenith range)
        deg2rad = np.pi / 180
        npix = len(relint)
        nside = hp.npix2nside(npix)
        theta, phi = hp.pix2ang(nside, range(npix))
        thetamax = (90 - opts['decmin']) * deg2rad
        thetamin = (90 - opts['decmax']) * deg2rad
        pass_dec_cut = (theta <= thetamax) * (theta >= thetamin)
        # Setup right-ascension bins
        ramin = opts['ramin'] * deg2rad
        ramax = opts['ramax'] * deg2rad
        rabins, rabinwidth = np.linspace(ramin, ramax,
            opts['nRAbins'] + 1, retstep=True)
        # Calculate phi for each pixel
        theta, phi = hp.pix2ang(nside, range(npix))
        # Bin in right ascension
        phiBins = np.digitize(phi, rabins) - 1
        # UNSEEN cut
        pass_UNSEEN_cut = (relint != hp.UNSEEN)
        pass_INF_cut = (relerr != np.inf)

        ri, ri_err = np.zeros((2, opts['nRAbins']))
        for i in range(opts['nRAbins']):
            pass_phiBin_cut = (phiBins == i)
            pass_all_cuts = pass_UNSEEN_cut * pass_phiBin_cut * pass_INF_cut * pass_dec_cut
            ri[i] = np.mean(relint[pass_all_cuts])
            ri_err[i] = np.sqrt(
                np.sum(relerr[pass_all_cuts]**2)) / pass_all_cuts.sum()
            #ri_err = np.sqrt(data * (bg + 1./20 * data)/bg**3)
        dx = rabinwidth / 2.
        ra = np.linspace(ramin + dx, ramax - dx, opts['nRAbins']) / deg2rad
        ra_err = dx * np.ones(opts['nRAbins']) / deg2rad

        return (ra, ri, ra_err, ri_err)

    def getProjFitParams(self, lmax, sigmay=None):
        ra, ri, ra_err, ri_err = self.getRIRAProj()
        popt, perr, chi2 = getHarmonicFitParams(ra, ri, lmax, sigmay=None)
        return popt, perr, chi2
