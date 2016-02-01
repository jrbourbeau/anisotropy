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

   for i in range(0,npix):
     (theta, phi) = H.pix2ang(nside,i)
     decs[i] = theta * 180. / pi


   fig = figure(figsize=(12,6))
   ax = SubplotHost(fig, 111)
   fig.add_subplot(ax)

   ax.scatter(decs,imap)

   show() 


