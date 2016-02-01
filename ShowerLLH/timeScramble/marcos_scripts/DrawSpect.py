#!/usr/bin/env python

from pylab import *
import numpy as n
import pyfits
from optparse import OptionParser
from mpl_toolkits.axes_grid.parasite_axes import SubplotHost
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib as mpl
import sys
import glob

#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#rc('text', usetex=True)
rc('xtick', labelsize=18)
rc('ytick', labelsize=18)
rcParams.update({ "font.family" : "serif" })

#-------------------------------------------------------

def containsAll(str, set):
   return 0 not in [c in str for c in set]

#-------------------------------------------------------

SigColor = 0
Color = 0

def AddSigmaBands(ax, file, plotstyle='slog', l2=False,col=0):
   clfile = n.genfromtxt(file)
   l = clfile[:,0]
   Cl = clfile[:,1]
   l2sigma = clfile[:,2]
   u2sigma = clfile[:,3]
   l1sigma = clfile[:,4]
   u1sigma = clfile[:,5]

   colors = [ \
      {'median': '#009900', '1sig': '#66FF00', '2sig': '#CCFF33'},  \
      {'median': '#990000', '1sig': '#FF6600', '2sig': '#FFCC33'}, \
      {'median': '#000099', '1sig': '#0066FF', '2sig': '#33CCFF'}, \
   ]

   z1 = [x*(x+1)/(2.*n.pi) for x in l]

   if l2:
      Cl = z1*Cl
      l2sigma = z1*l2sigma
      l1sigma = z1*l1sigma
      u1sigma = z1*u1sigma
      u2sigma = z1*u2sigma

   if plotstyle == 'loglog':
      ax.loglog(l,Cl,linewidth=4,color=colors[col]['median'],label=file)
   elif plotstyle == 'slog':
      ax.semilogy(l,Cl,linewidth=4,color=colors[col]['median'],label=file)
   else:
      ax.plot(l,Cl,linewidth=4,color=colors[col]['median'],label=file)

   ax.fill_between(l,l2sigma,u2sigma,facecolor=colors[col]['2sig'],linewidth=0)
   ax.fill_between(l,l1sigma,u1sigma,facecolor=colors[col]['1sig'],linewidth=0)

#-------------------------------------------------------

def AddRandomBands(ax,plotstyle='slog',b=0):
   ''' 
   Not really random bands, but 1 and 2 sigma bands from the power spectrum 
   of random 10000 maps
   '''
   
   clfile = n.genfromtxt('/Users/santander/icecube/anisotropy/tools/etc/bands3.cl')

   if b == 1:
      clfile = n.genfromtxt('/Users/santander/icecube/anisotropy/tools/etc/bands2.cl')


   l = clfile[:,0]
   Cl = clfile[:,1]
   l2sigma = clfile[:,2]
   u2sigma = clfile[:,3]
   l1sigma = clfile[:,4]
   u1sigma = clfile[:,5]

   # Blues

   if b == 1:
     ax.fill_between(l,l2sigma,u2sigma,facecolor='#33CCFF',linewidth=0,alpha=0.8)
     ax.fill_between(l,l1sigma,u1sigma,facecolor='#0066FF',linewidth=0,alpha=0.8)
   
   # Oranges
   else:
     ax.fill_between(l,l2sigma,u2sigma,facecolor='#FFCC33',linewidth=0)
     ax.fill_between(l,l1sigma,u1sigma,facecolor='#FF6600',linewidth=0)
   
   # Grays
   #ax.fill_between(l,l2sigma,u2sigma,facecolor='#CCCCCC',linewidth=0)
   #ax.fill_between(l,l1sigma,u1sigma,facecolor='#999999',linewidth=0)
   
   if plotstyle == 'loglog':
      ax.loglog(l,Cl,linewidth=4,color='#000099')
   elif plotstyle == 'slog':  
      if b == 1:
         ax.semilogy(l,Cl,linewidth=3,color='#000099')
      else:
         ax.semilogy(l,Cl,linewidth=3,color='#990000')
      #ax.semilogy(l,Cl,linewidth=4,color='#990000')
      #ax.semilogy(l,Cl,linewidth=4,color='#000000',nonposy='mask')
   else:   
      ax.plot(l,Cl,linewidth=4,color='#000099')

#-------------------------------------------------------

def MultipoleToAngle(l):
    """
    Convert multipole moment to angular scale
    """
    return (180. / l)

#-------------------------------------------------------

def TickToLabel(tick):
    """Convert a floating point tick value to a string label
    """
    if abs(tick - int(tick)) > 0:
        return ("%.1f" % tick)
    else:
        return ("%d" % int(tick))

#-------------------------------------------------------

def SetupAngleAxis(ax, lticks=None):
    """Create the angle axis
    """
    # Convert ticks (multipoles) to degrees
    if not lticks:
        lticks = [1.8, 3.6, 7.2, 18., 36., 90, 180.]
    dticks = [MultipoleToAngle(l) for l in lticks]
    tickLabels = [TickToLabel(d) for d in dticks]

    ax2 = ax.twin()
    ax2.set_xticks(lticks)
    ax2.set_xticklabels(tickLabels)
    ax2.axis["top"].set_label("Angular scale [$\circ$]")
    ax2.axis["top"].label.set_visible(True)
    ax2.axis["top"].label.set_size(16)
    ax2.axis["top"].major_ticklabels.set_fontsize(16)
    ax2.axis["right"].major_ticks.set_visible(False)
    ax2.axis["right"].major_ticklabels.set_visible(False)

#-------------------------------------------------------

def ReadSpectrum(str, ax, plotstyle, l2=False):
   l = []
   Cl = []
   if str.startswith('-') or str.isdigit():
      return
   elif str.endswith('.cl'):
      print "Reading Cl file ", str
      clfile = n.genfromtxt(str)
      l = clfile[:,0]
      Cl = clfile[:,1]
      if len(clfile[0]) > 2:
         global SigColor
         AddSigmaBands(ax,str,plotstyle,l2,col=SigColor)
         SigColor += 1
         return
   
   elif str.endswith('.fits'):      
      print "Reading FITS file ", str
      fitsfile = pyfits.open(str)
      if len(fitsfile[1].columns.names) < 2 and not containsAll(fitsfile[1].columns.names[0],'C_l'):
         print "ERROR, the FITS file given (",str,") doesn't seem to be an spectrum file"
         return
      tbdata = fitsfile[1].data
      Cl = tbdata.field('Temperature C_l')
      l = range(0,len(Cl),1)

      if not containsAll(fitsfile[1].columns.names[0],'C_l'):
         l = tbdata.field('l')
      else:
         print "No l column found, building l indices from Cl values"
   else:
      print "Uh? The file ",str," doesn't seem to be a Cl nor FITS file containing a spectrum"
      return

   z1 = [x*(x+1)/(2.*n.pi) for x in l]   
  
   if l2:
      Cl = z1*Cl
   
   # If I want to define some colors
   #colors=['black','magenta','red']

   if plotstyle == 'loglog':
      ax.loglog(l,Cl,'-',label=str,linewidth=4,markersize=12)
   elif plotstyle == 'slog':
      global Color
      #ax.semilogy(l,Cl,'o-',label=str,linewidth=3,color=colors[Color],markersize=8)
      ax.semilogy(l,Cl,'-',label=str,linewidth=3,markersize=12)
      #leg=str.split('_')[1]
      #ax.semilogy(l,Cl,label=leg,linewidth=4)
      Color += 1
   else:
      ax.plot(l,Cl,'-',label=str,linewidth=4,markersize=12)


#-------------------------------------------------------
#-------------------------------------------------------

if __name__ == "__main__":
   # Set up command line options
   usage = "usage: %prog [options] [list of Cl or FITS files with spectra]"
   parser = OptionParser(usage)

   #majorFormatterX = FormatStrFormatter('%d')
   #majorFormatterY = FormatStrFormatter('%f')

   parser.add_option("-L", "--lmax", dest="lmax", type=float, help="Maximum multipole", default = 40)
   parser.add_option("-l", "--lmin", dest="lmin", type=float, help="Minimum multipole", default = 0)
   parser.add_option("-p", "--pstyle", dest="pstyle", type=int, help="Plot style [0 - Semilogy, 1 - Loglog, 2 - Linear]", default = 0)
   parser.add_option("-s", "--scale", dest="l2", help="l(l+1) scale", action="store_true", default = False)
   parser.add_option("-b", "--bands", dest="band", help="1 and 2 sigma bands", action="store_true", default = False)
   parser.add_option("-c", "--colors", dest="colors", help="spectrum colors", action="store_true", default = False)
   parser.add_option("-M", "--maxCl", dest="Clmax", type=float, help="Maximum Cl", default = -1)
   parser.add_option("-m", "--minCl", dest="Clmin", type=float, help="Minimum Cl", default = -1)

   pstyles = ['slog','loglog','lin']
    
   (options, args) = parser.parse_args()
   if len(args) < 1:
      parser.error("incorrect number of arguments")
   
   if options.pstyle > 2:
      parser.error("plot style error, possible values are: 0 - Semilogy, 1 - Loglog, 2 - Linear")

   fig = figure(figsize=(12,6))
   ax = SubplotHost(fig, 111)
   fig.add_subplot(ax)

   if options.band:
      AddRandomBands(ax,pstyles[options.pstyle],1)
      #AddRandomBands(ax,pstyles[options.pstyle],1)
  
   if options.colors:
      nColors = len([x for x in sys.argv[1:] if containsAll(x,'.cl') or containsAll(x,'.fits')])
      print nColors
      cmap = mpl.cm.get_cmap(name='jet')
      mycolors = [cmap(i) for i in n.linspace(0, 0.9, nColors)]
      ax.set_color_cycle(mycolors)

   for str in sys.argv[1:]:
      if containsAll(str,'*'):
         files = glob.glob(str)
         for file in files:
            ReadSpectrum(file,ax,pstyles[options.pstyle],options.l2)
      else:
         ReadSpectrum(str,ax,pstyles[options.pstyle],options.l2)

   ax.legend()

   xlabel(r"$\ell$",fontsize=16)
   if options.l2:
      ylabel(r"$l(l+1)/ 2 \pi \; C_{\ell}$",fontsize=36)
   else:
      ylabel(r"$C_{\ell}$",fontsize=16)

   xticks = [1, 3, 5, 10, 15, 20, 30, 100]

   if (options.lmax < 100):  
     ax.set_xticks(xticks)
     ax.set_xticklabels(xticks)
     SetupAngleAxis(ax, xticks)

   ax.axis["bottom"].label.set_size(16)
   ax.axis["bottom"].major_ticklabels.set_fontsize(16)
   ax.axis["left"].label.set_size(16)
   ax.axis["left"].major_ticklabels.set_fontsize(16)

   #ax.yaxis.set_major_formatter(majorFormatterY)
   #ax.xaxis.set_major_formatter(majorFormatterX)

   if (options.Clmin != -1 and options.Clmax != -1):  
     ax.set_ylim(options.Clmin,options.Clmax)

   ax.set_xlim(options.lmin,options.lmax)


   ax.grid()

   #draw()
   show()

#-------------------------------------------------------
#-------------------------------------------------------

