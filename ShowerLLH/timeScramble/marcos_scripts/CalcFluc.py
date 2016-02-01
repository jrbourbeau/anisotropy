#!/usr/bin/env python

import healpy as H
import sys
from numpy import *

if __name__ == "__main__":

   data = H.read_map(sys.argv[1])
   bkg = H.read_map(sys.argv[2])

   nside = H.get_nside(data)
   npix = H.nside2npix(nside)
   map = zeros(npix, dtype=double)

   for i in range(0,npix):
     if (bkg[i] != 0) :
        map[i] = data[i] / bkg[i] - 1
        
   H.write_map("flucmap.fits", map)


