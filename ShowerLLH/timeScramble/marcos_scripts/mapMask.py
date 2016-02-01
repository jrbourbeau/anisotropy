#!/usr/bin/env python

import healpy as H
import sys
from numpy import *

if __name__ == "__main__":

   degree = pi / 180.

   nside = int(sys.argv[1])
   thetai = float(sys.argv[2]) * degree
   thetaf = float(sys.argv[3]) * degree

   print "nside: ", nside
   print "low edge: ", thetai/degree, " hi edge: ", thetaf/degree

   npix = H.nside2npix(nside)
   omap = zeros(npix, dtype=double)

   for i in range(0,npix):
     (theta, phi) = H.pix2ang(nside,i)
     if theta >= thetai and theta <= thetaf:
        omap[i] = 1.

   H.write_map("mask.fits", omap)


