#!/usr/bin/env python

import healpy as H
import sys
from numpy import *
from pylab import * 
import matplotlib as mpl

if __name__ == "__main__":

   imap = H.read_map(sys.argv[1])

   nside = H.get_nside(imap)
   npix = H.nside2npix(nside)

   counts = zeros(4*nside, dtype=float)
   norm = zeros(4*nside, dtype=float)

   for i in range(0,npix):
     (x, y, z) = H.pix2vec(nside,i)
     ringId = H.ring_num(nside,z)

     counts[ringId] += imap[i]
     norm[ringId] += 1.

   omap = zeros(npix, dtype=float)
   mapmax = max(imap)
   
   for i in range(0,npix):
     (x, y, z) = H.pix2vec(nside,i)
     ringId = H.ring_num(nside,z)

     omap[i] = counts[ringId] / norm[ringId]  
      
   
   H.write_map("mavg_bg.fits",omap)
   H.write_map("mavg_data.fits",imap)
   
      
   




