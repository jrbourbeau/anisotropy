#!/usr/bin/env python

#=============================================================================
# File Name     : realSH.py
# Description   : Script to create real, normalized spherical harmonic maps
#                 using scipy's sph_harm() function. These maps are stored
#                 for quick access later.
# Creation Date : 03-17-2016
# Last Modified : Thu 17 Mar 2016 05:12:40 PM CDT
# Created By    : James Bourbeau
#=============================================================================

import numpy as np
import healpy as hp
import argparse
from scipy.special import sph_harm

def realSH(l,m,theta,phi):
    if m < 0:
        real_SH = (1j/np.sqrt(2))*(sph_harm(m,l,phi,theta)
            -((-1)**m)*(sph_harm(-m,l,phi,theta)))
    elif m == 0:
        real_SH = sph_harm(m,l,phi,theta)
    else:
        real_SH = (1/np.sqrt(2))*(sph_harm(-m,l,phi,theta)
            +((-1)**m)*(sph_harm(m,l,phi,theta)))
    if type(real_SH)!=np.ndarray:
        real_SH = np.array([real_SH])
        return real_SH.real
    if np.isreal(real_SH).all():
        return real_SH.real
    elif (np.absolute([y for y in real_SH.imag if y!=0.])<1e-15).all():
        return real_SH.real
    else:
        raise SystemError('realSH is not returning a real value. Something has gone terribly wrong...')

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Calculate anisotropy')
    parser.add_argument('--lmax', dest='lmax', type=int, default=3, \
        help='Maximum l-value used to generate SH maps.')
    parser.add_argument('-o','--outdir', dest='outdir', \
        help='Output directory where .npy files will be stored.')

    args = parser.parse_args()
    opts = vars(args).copy()

    # Useful breakdown of l, m values to be used
    l = opts['lmax']
    nsph = sum([2*l_i+1 for l_i in range(l+1)])
    lvals = [[l_i]*(2*l_i+1) for l_i in range(l+1)]
    mvals = [[m for m in range(-l_i, l_i+1)] for l_i in range(l+1)]
    lvals = [item for sublist in lvals for item in sublist]
    mvals = [item for sublist in mvals for item in sublist]

    nside = 64
    npix = hp.nside2npix(nside)
    theta, phi = hp.pix2ang(nside, range(npix))
    fitparams = ['Y(%i,%i)' % (lvals[i], mvals[i]) for i in range(nsph)]
    SH = {}
    for i in range(len(lvals)):
        key = 'Y({},{})'.format(lvals[i], mvals[i])
        print('Generating {}...'.format(key))
        SH[key] = realSH(lvals[i],mvals[i],theta,phi)

    np.save(opts['outdir']+'normedSH.npy',SH)
    print('\nSphereical harmonics up to l={} saved\n'.format(opts['lmax']))
