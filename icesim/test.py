#!/usr/bin/env python

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import argparse

from icecube import photospline

def test(s, htype='2d'):

    compDict = {'Proton':2212, 'Iron':1000260560}
    ecut = (s['energy_sp'] >= 4) * (s['energy_sp'] < 4.75)

    # Energy binning
    etrue = np.log10(s['energy'])
    ebins = np.linspace(etrue.min(), etrue.max(), 101)

    # Declination binning
    degree = np.pi / 180.
    dec = s['zen'] / degree - 90
    thetamax = 90.
    dbins = np.linspace(1, np.cos(thetamax*degree), 51)
    dbins = np.arccos(dbins) / degree - 90

    for comp in compDict:

        cut = ecut * (s['type'] == compDict[comp])

        fig, ax = plt.subplots()
        ax.set_title('%s Distribution' % comp)
        if htype == 'dec':
            h = ax.hist(dec[cut], bins=dbins, histtype='step')
        if htype == 'energy':
            h = ax.hist(etrue[cut], bins=ebins, histtype='step')
        if htype == '2d':
            h = ax.hist2d(etrue[cut], dec[cut], bins=[ebins, dbins])
            ax.set_xlabel(r'True Energy ($log_{10}(E/GeV)$)')
            ax.set_ylabel(r'Declination ($^{\circ}$)')
            fig.colorbar(h[3], ax=ax)
        plt.show()


if __name__ == "__main__":

    p = argparse.ArgumentParser()
    p.add_argument('--htype', dest='htype',
            default='2d',
            help='Histogram type [2d | energy | dec]')
    args = p.parse_args()

    simFile = '/data/user/fmcnally/anisotropy/sim/IC86_10649.npy'
    s = np.load(simFile)
    s = s.item()

    # Get reconstructed direction from DST
    theta, phi = hp.pix2ang(1024, s['reco2'].astype('int'))
    s['dstZenith']  = theta
    s['dstAzimuth'] = phi

    # Calculate splined energy values
    print 'Calculating energy values from spline tables...'
    splineFile = '/data/user/fmcnally/anisotropy/sim/IC86_10649_spline.fits'
    x = np.cos(s['dstZenith'])
    with np.errstate(divide='ignore'):
        y = np.log10(s['nchannel'])
    stable = photospline.I3SplineTable(splineFile)
    def spline_energy(x, y):
        if (x < 0.3):
            return 0
        return stable.eval([x, y])
    s['energy_sp'] = np.array(map(spline_energy, x, y))

    test(s, htype=args.htype)





