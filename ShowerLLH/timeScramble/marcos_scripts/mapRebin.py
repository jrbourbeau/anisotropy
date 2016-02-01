#!/usr/bin/env python

import healpy as H
import sys
from numpy import *

if __name__ == "__main__":

   map = H.read_map(sys.argv[1])
   fnside = int(sys.argv[2])

   inside = H.get_nside(map)
   npix = H.nside2npix(inside)
   onpix = H.nside2npix(fnside)
   omap = zeros(onpix, dtype=double)

   for i in range(0,npix):
     (itheta, iphi) = H.pix2ang(inside,i)  
     opix = H.ang2pix(fnside,itheta,iphi)
     omap[opix] += map[i]

   H.write_map("udmap.fits", omap)


