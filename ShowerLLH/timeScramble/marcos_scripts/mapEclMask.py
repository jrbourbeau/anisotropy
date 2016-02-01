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

   thetae = 23.4392794 * degree
   phie = 180 * degree
   epole = H.ang2pix(nside,thetae,phie)
   
   a = H.pix2vec(nside,epole)

   for i in range(0,npix):
     b = H.pix2vec(nside,i)
     delta = arccos(dot(a,b))

     if delta >= thetai and delta <= thetaf:
        omap[i] = 1.

   H.write_map("emask.fits", omap)


