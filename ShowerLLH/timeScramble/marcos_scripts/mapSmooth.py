#!/usr/bin/env python

import healpy as H
import sys
from numpy import *

def LiMaSignificance(Non, Nback, alpha):
   Noff = Nback / alpha;
   sign = 1
   if Non < alpha * Noff: 
      sign = -1
   return sign * sqrt(2*(Non *log(((1+alpha)*Non)/(alpha * (Non + Noff))) + Noff *log(((1+alpha)*Noff)/(Non+Noff))))    

if __name__ == "__main__":

   basename = sys.argv[1]
   smooth = int(sys.argv[2])

   bkgmap = H.read_map(basename + "_bg.fits")
   datmap = H.read_map(basename + "_data.fits")

   nside1 = H.get_nside(bkgmap)
   nside2 = H.get_nside(datmap)

   if nside1 != nside2: 
      print "ERROR: maps have different Nsides!"
      sys.exit()
      
   npix = H.nside2npix(nside1)
   relint = zeros(npix, dtype=double)
   signif = zeros(npix, dtype=double)

   print "smoothing"

   for i in range(0,npix):
     v = H.pix2vec(nside1, i)
     pixlist = H.query_disc(nside1, v, smooth)

     datval = [datmap[j] for j in pixlist]
     bkgval = [bkgmap[j] for j in pixlist]

     sumbkg = sum(bkgval)
     sumdat = sum(datval)

     relint[i] = sumdat / sumbkg - 1
     signif[i] = LiMaSignificance(sumdat, sumbkg, 1/20.)


   sigmafile = basename + "_sigma_" + smooth + "deg.fits"
   relintfile = basename + "_ri_" + smooth + "deg.fits"

   H.write_map(sigmafile, signif)
   H.write_map(relintfile, relint)


