#!/usr/bin/python

from pylab import *
import numpy as n
import pyfits
from optparse import OptionParser
from mpl_toolkits.axes_grid.parasite_axes import SubplotHost
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib as mpl
import sys
import glob
import matplotlib.font_manager as font_manager

rc('xtick', labelsize=18)
rc('ytick', labelsize=18)

#-------------------------------------------------------

def containsAll(str, set):
   return 0 not in [c in str for c in set]

#-------------------------------------------------------

def AddErrorBars(ax, filename, col=0):
   file = n.genfromtxt(filename)
   x = file[:,0]
   y = file[:,1]
   sigma = file[:,2]

   ax.errorbar(x,y,yerr=[sigma,sigma],fmt="r.",capsize=0,label="",linewidth=2,markersize=7,color='black',marker='o',mec='black')

#-------------------------------------------------------

def ReadFile(str, ax):
   l = []
   Cl = []
   if str.startswith('-') or str.isdigit():
      return
   elif str.endswith('.dat'):
      print "Reading dat file ", str
      AddErrorBars(ax,str)
      return

   else:
      print "Uh? The file ",str," doesn't seem to contain relative intensity data"
      return

#-------------------------------------------------------
#-------------------------------------------------------

if __name__ == "__main__":
   # Set up command line options
   usage = "usage: %prog [options] [list of Cl or FITS files with spectra]"
   parser = OptionParser(usage)

   majorFormatterX = FormatStrFormatter('%d')
   #majorFormatterY = FormatStrFormatter('%s')

   parser.add_option("-X", "--xmax", dest="xmax", type=float, help="Maximum RA", default = 360)
   parser.add_option("-x", "--xmin", dest="xmin", type=float, help="Minimum RA", default = 0)
   parser.add_option("-Y", "--ymax", dest="ymax", type=float, help="Maximum RI", default = 1e-3)
   parser.add_option("-y", "--ymin", dest="ymin", type=float, help="Minimum RI", default = -1e-3)

   (options, args) = parser.parse_args()
   if len(args) < 1:
      parser.error("incorrect number of arguments")
   
   fig = figure(figsize=(9,8))
   #fig.subplots_adjust(0.1,0.1,0.97)
   ax = SubplotHost(fig, 111)
   fig.add_subplot(ax)

   for str in sys.argv[1:]:
      if containsAll(str,'*'):
         files = glob.glob(str)
         for file in files:
            ReadFile(file,ax)
      else:
         ReadFile(str,ax)

   xlabel(r"Right Ascension $[^{\circ}]$",fontsize=18)
   ylabel(r"$\Delta N/\langle N \rangle$",fontsize=60)

   #ax.axis["bottom"].label.set_size(18)
   #ax.axis["bottom"].major_ticklabels.set_fontsize(18)
   #ax.axis["left"].label.set_size(18)
   #ax.axis["left"].major_ticklabels.set_fontsize(18)
   
   xzero = arange(0, 360, 1)
   yzero = 0 * xzero
   ax.plot(xzero,yzero,linewidth=2)

   ax.set_ylim(options.ymin,options.ymax)
   ax.set_xlim(options.xmin,options.xmax)



   #ax.grid()

   #draw()
   show()

#-------------------------------------------------------
#-------------------------------------------------------

