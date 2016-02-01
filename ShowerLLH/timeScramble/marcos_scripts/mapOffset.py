#!/usr/bin/env python

import healpy as H
import sys
from numpy import *

if __name__ == "__main__":

   map1 = H.read_map(sys.argv[1])
   offset = float(sys.argv[2])

   nside1 = H.get_nside(map1)

   npix = H.nside2npix(nside1)
   omap = zeros(npix, dtype=double)

   for i in range(0,npix):
     omap[i] = map1[i] - offset

   H.write_map("moffset.fits", omap)


