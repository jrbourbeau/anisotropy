#!/usr/bin/python

from pylab import *
import numpy as n
import pyfits
from optparse import OptionParser
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib as mpl
import sys
import glob
import matplotlib.font_manager as font_manager
import healpy as H

#rc('xtick', labelsize=18)
#rc('ytick', labelsize=18)

#-------------------------------------------------------

def containsAll(str, set):
   return 0 not in [c in str for c in set]

def sigmaRI(data, bkg, alpha):
   return sqrt(data * (bkg + alpha * data)/bkg**3)

#-------------------------------------------------------
#-------------------------------------------------------

if __name__ == "__main__":

   # Set up command line options
   usage = "usage: %prog [options] bkg.fits data.fits"
   parser = OptionParser(usage)
   parser.add_option("-r", "--ramin", dest="ramin", type=float,
                     default=0, help="minimum RA")
   parser.add_option("-R", "--ramax", dest="ramax", type=float,
                     default=360, help="maximum RA")
   parser.add_option("-D", "--decmax", dest="decmax", type=float,
                     default=-35, help="maximum Dec")
   parser.add_option("-d", "--decmin", dest="decmin", type=float,
                     default=-75, help="minimum Dec")
   parser.add_option("-a", "--alpha", dest="alpha", type=float,
                     default=20, help="1/alpha")
   parser.add_option("-n", "--nbins", dest="nbins", type=int,
                     default=24, help="number of bins")
   #parser.add_option("-I", "--", dest="decmin", type=float,
   #                  help="minimum Dec")
   #parser.add_option("-i", "--decmin", dest="decmin", type=float,
   #                  help="minimum Dec")
   parser.add_option("-z","--zero", action="store_true", dest="zeroline",
                     default=False, help="Draw zero line")
   parser.add_option("-f","--flipra", action="store_true", dest="flipra",
                     default=False, help="Flips RA in x axis")

   (options, args) = parser.parse_args()

   print options

   degree = n.pi / 180
   
   bkgmap = H.read_map(args[0])
   datamap = H.read_map(args[1])

   ramin = options.ramin * degree
   ramax = options.ramax * degree
   decmax = options.decmax 
   decmin = options.decmin
   nbins = options.nbins
   
   thetamin = (90 - decmax) * degree
   thetamax = (90 - decmin) * degree

   bins = n.linspace(ramin, ramax, nbins+1)
   data = n.zeros(nbins)
   bkg  = n.zeros(nbins)

   nside = H.get_nside(bkgmap)
   npix = H.nside2npix(nside)
   
   for i in range(0,npix):
     (theta, phi) = H.pix2ang(nside,i)
     if theta >= thetamin and theta <= thetamax:
        data[n.digitize([phi],bins)-1] += datamap[i]
        bkg[n.digitize([phi],bins)-1] += bkgmap[i]
   
   ri = (data - bkg) / bkg     
   sigmay = sigmaRI(data, bkg, 1/options.alpha)
   dx = (ramax - ramin)/(2*nbins)
   ra = n.linspace(ramin+dx, ramax-dx, nbins) / degree 
   sigmax = dx * n.ones(nbins) / degree

   fig = figure(figsize=(8,6))
   ax = fig.add_subplot(111)

   ax.errorbar(ra,ri,xerr=[sigmax,sigmax],yerr=[sigmay,sigmay],fmt="r.",capsize=0,label="",linewidth=1.5,markersize=4,color='blue',marker='o',mec='blue')
   #ax.fill_between(ra,ri-sigmay,ri+sigmay,facecolor='#66FF00',linewidth=0)
   
   xlabel(r"Right Ascension $[^{\circ}]$",fontsize=14)
   ylabel(r"$\Delta N/\langle N \rangle$",fontsize=14)
   
   if options.zeroline:
      xzero = arange(0, 360, 1)
      yzero = 0 * xzero
      ax.plot(xzero,yzero,linewidth=1.5,linestyle='--',color='black')

   #ax.set_ylim(-0.00105,0.00105)
   if options.flipra:
      ax.set_xlim(options.ramin,options.ramax)
   else:   
      ax.set_xlim(options.ramax,options.ramin)

   grid()

   draw()
   show()

#-------------------------------------------------------
#-------------------------------------------------------

