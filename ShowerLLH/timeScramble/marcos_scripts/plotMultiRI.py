#!/usr/bin/env python

from pylab import *
import numpy as n
import pyfits
from optparse import OptionParser
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib as mpl
from matplotlib import cm
import sys
import glob
import matplotlib.font_manager as font_manager
import healpy as H

params = {'legend.fontsize': 10 }
rcParams.update(params)

#rc('xtick', labelsize=18)
#rc('ytick', labelsize=18)

#-------------------------------------------------------

def containsAll(str, set):
   return 0 not in [c in str for c in set]

def sigmaRI(data, bkg, alpha):
   return sqrt(data * (bkg + alpha * data)/bkg**3)

def plotFit(ax, fitfile, options):
   fit = H.read_map(fitfile)
   
   ramin = options.ramin * degree
   ramax = options.ramax * degree
   decmax = options.decmax 
   decmin = options.decmin
   nbins = options.nbins
   
   thetamin = (90 - decmax) * degree
   thetamax = (90 - decmin) * degree

   bins = n.linspace(ramin, ramax, nbins+1)
   data = n.zeros(nbins)
   N = n.zeros(nbins)
   
   nside = H.get_nside(fit)
   npix = H.nside2npix(nside)

   for i in range(0,npix):
     (theta, phi) = H.pix2ang(nside,i)
     if theta >= thetamin and theta <= thetamax:
        data[n.digitize([phi],bins)-1] += fit[i]
        N[n.digitize([phi],bins)-1] += 1

   ri = data / N
   dx = (ramax - ramin)/(2*nbins)
   ra = n.linspace(ramin+dx, ramax-dx, nbins) / degree

   ax.plot(ra,ri,linewidth=1.5)

def getFiles(basename): 
   bkgname = basename + "_bg.fits"
   dataname = basename + "_data.fits"
   
   bkgmap = H.read_map(bkgname)
   datamap = H.read_map(dataname)
   return (bkgmap, datamap)
   

def calcRI(ax, bkgmap, datamap, basename, options): 
   #bkgname = basename + "_bg.fits"
   #dataname = basename + "_data.fits"

   #bkgmap = H.read_map(bkgname)
   #datamap = H.read_map(dataname)

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

   ax.errorbar(ra,1e3*ri,xerr=[sigmax,sigmax],yerr=[1e3*sigmay,1e3*sigmay],capsize=0,label=basename,linewidth=2,markersize=5,marker='o',fmt='.',
   mew=1, mec='none'
   )
   #ax.fill_between(ra,ri-sigmay,ri+sigmay,linewidth=0,label=basename,alpha=0.5)


def calcAvgBkg(ax, basename, options):
   dataname = basename + "_data.fits"
   datamap = H.read_map(dataname)

   nside = H.get_nside(datamap)
   npix = H.nside2npix(nside)

   counts = zeros(4*nside, dtype=float)
   norm = zeros(4*nside, dtype=float)

   for i in range(0,npix):
     (x, y, z) = H.pix2vec(nside,i)
     ringId = H.ring_num(nside,z)

     counts[ringId] += datamap[i]
     norm[ringId] += 1.

   bkgmap = zeros(npix, dtype=float)
   mapmax = max(datamap)

   for i in range(0,npix):
     (x, y, z) = H.pix2vec(nside,i)
     ringId = H.ring_num(nside,z)

     bkgmap[i] = counts[ringId] / norm[ringId]
   
   calcRI(ax, bkgmap, datamap, basename + " (avg)", options) 


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
   parser.add_option("-t", "--fit", dest="fit",
                     default=None, help="Fit file")
   parser.add_option("-I", "--rimax", dest="rimax", type=float,
                     default=None, help="max RI")
   parser.add_option("-i", "--rimin", dest="rimin", type=float,
                     default=None, help="min RI")
   parser.add_option("-z","--zero", action="store_true", dest="zeroline",
                     default=False, help="Draw zero line")
   parser.add_option("-f","--flipra", action="store_true", dest="flipra",
                     default=False, help="Flips RA in x axis")
   parser.add_option("-v","--avg", action="store_true", dest="avgbkg",
                     default=False, help="Bkg from avg data")
   parser.add_option("-o", "--output", dest="output", default=None,
                      help="Output image file")
   parser.add_option("-b","--batchmode", action="store_true", dest="batchMode",
                      default=False, help="Execute without interaction")
   parser.add_option("-p", "--paper", action="store_true", dest="paper",
                     default=False, help="paper style")

   (options, args) = parser.parse_args()


   if options.paper:
        rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
        rc('text', usetex=True)
        rc('xtick', labelsize=16)
        rc('ytick', labelsize=16)

   print options

   degree = n.pi / 180

   fig = figure(figsize=(8,6))
   ax = fig.add_subplot(111)

   for str in args:
     (bgmap, datamap) = getFiles(str)
     calcRI(ax, bgmap, datamap, str, options)
     if options.avgbkg:
        calcAvgBkg(ax, str, options)

   if options.fit != None:
      plotFit(ax, options.fit, options)
   
   xlabel(r"Right Ascension $[^{\circ}]$",fontsize=14)
   ylabel(r"$\Delta N/\langle N \rangle$",fontsize=14)

   if options.zeroline:
      xzero = arange(0, 360, 1)
      yzero = 0 * xzero
      ax.plot(xzero,yzero,linewidth=1.5,linestyle='--',color='black')

   if options.flipra:
      ax.set_xlim(options.ramin,options.ramax)
   else:   
      ax.set_xlim(options.ramax,options.ramin)
      
   ax.set_ylim(options.rimin,options.rimax)

   #grid()
   legend(loc=2, numpoints=1)

   if options.paper:
      majorFormatterX = FormatStrFormatter('%d')
      majorFormatterY = FormatStrFormatter('%0.1f')

      ax.yaxis.set_major_formatter(majorFormatterY)
      ax.xaxis.set_major_formatter(majorFormatterX)

   draw()

   if options.output:
     if options.output.find("jpg") > 0:
           fig.savefig(options.output, dpi=300)
     else:
           savefig(options.output)
           jpgOutput = options.output.replace("png", "jpg")
           fig.savefig(jpgOutput, dpi=300)

   if not options.batchMode:
      show()

#-------------------------------------------------------
#-------------------------------------------------------

