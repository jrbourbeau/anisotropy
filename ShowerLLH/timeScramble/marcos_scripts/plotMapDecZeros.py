#!/usr/bin/env python

import healpy as H
import sys
from numpy import *
from pylab import * 
from mpl_toolkits.axes_grid.parasite_axes import SubplotHost
import matplotlib as mpl

if __name__ == "__main__":

   imap = H.read_map(sys.argv[1])

   nside = H.get_nside(imap)
   npix = H.nside2npix(nside)

   decs = zeros(npix, dtype=float)
   nzmap = zeros(npix, dtype=float)

   n = 0

   for i in range(0,npix):
     (theta, phi) = H.pix2ang(nside,i)
     if theta > pi/2 and imap[i] != 0:
        decs[n] = 180 - theta * 180 / pi 
        nzmap[n] = imap[i]
        n = n+1

   fig = figure(figsize=(12,6))
   ax = SubplotHost(fig, 111)
   fig.add_subplot(ax)

   print size(decs)

   ax.semilogy(decs,nzmap,'o')
   #ax.scatter(decs,nzmap)
   xrange(0,65)
   xlabel("Zenith angle [deg]")
   ylabel("Rel intensity")

   show() 


