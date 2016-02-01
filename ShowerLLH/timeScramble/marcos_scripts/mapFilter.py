#!/usr/bin/env python

import healpy as H
import sys
from numpy import *

if __name__ == "__main__":

   map = H.read_map(sys.argv[1])

   nside1 = H.get_nside(map)

   npix = H.nside2npix(nside1)
   omap = zeros(npix, dtype=double)
   
   for i in range(0,npix):
     if (map[i] > -20 and map[i] < 30):
       omap[i] = map[i]

   H.write_map("mfilter.fits", omap)
   

