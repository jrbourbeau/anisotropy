#!/usr/bin/env python

import numpy as np
import healpy as hp
import glob

def maskFile(file):

    m = hp.read_map(file)

    # Create mask
    nPix = len(m)
    nSide = hp.npix2nside(nPix)
    theta, phi = hp.pix2ang(nSide, range(nPix))
    theta = np.pi - theta
    mask = np.cos(theta)>=0.8
    masked_pix = np.logical_not(mask)

    m[masked_pix] = hp.UNSEEN

    hp.write_map(file, m)


if __name__ == "__main__":

    prefix = '/net/user/fmcnally/ShowerLLH/maps/smoothed/'
    fileList = glob.glob(prefix + '*.fits')
    for file in fileList:
        maskFile(file)
