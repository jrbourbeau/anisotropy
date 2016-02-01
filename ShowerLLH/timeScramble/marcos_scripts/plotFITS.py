#!/usr/bin/env python

from pylab import *
from optparse import OptionParser
import healpy as H

degree = pi / 180.0

def IsInIC40Mask(alpha, delta):
    """Return true if coordinates are within the IC40 mask region 
    """
    #return delta >= (-85*degree) and delta <= (-5*degree)
    return delta <= (-25*degree)
    #return delta <= (0*degree)

def IsInMilagroMask(alpha, delta):
    """Return true if coordinates are within the Milagro mask region
    """
    #return delta >= (-6*degree) and delta <= (63*degree)
    return delta >= (5*degree) and delta <= (63*degree)

def IsInCombinedMask(alpha, delta):
    return (delta >= (5*degree) and delta <= (63*degree)) or delta <= (-25*degree)

def TickToLabel(tick):
    """Convert a floating point tick value to a string label
    """
    if abs(tick - int(tick)) > 0:
        return ("%.1f" % tick)
    else:
        return ("%d" % int(tick))

def SetupColorBar(fig, title=None, ticks=None, coord="C", geo=False):
    """Create the color bar for a HEALPix figure
    """
    fig = figure(1)
    for ax in fig.get_axes():
        if type(ax) is H.projaxes.HpxMollweideAxes:
            tickLabels = []
            if ticks:
                tickLabels = [TickToLabel(t) for t in ticks]
            cb = fig.colorbar(ax.get_images()[0], ax=ax,
                              orientation="horizontal",
                              #orientation="vertical",
                              shrink=0.8, aspect=35,
                              pad=0.05, fraction=0.1,
                              ticks=ticks)
            cb.ax.set_xticklabels(tickLabels)
            cb.set_label(title)
            if coord == "C" and geo == False:
                ax.annotate("0h", xy=(1.8, -0.625))
                ax.annotate("24h", xy=(-1.9, -0.625))
            elif coord != "C" and geo == False:
                ax.annotate("-180$^\circ$", xy=(1.7, -0.635))
                ax.annotate("180$^\circ$", xy=(-1.9, -0.635))
            elif geo == True:
                ax.annotate("$360^{\circ}$", xy=(1.8, -0.625))
                #ax.annotate("$90^{\circ}$", xy=(0, -0.9))
                ax.annotate("$0^{\circ}$", xy=(-1.9, -0.625))


if __name__ == "__main__":
    # Set up command line options
    usage = "usage: %prog [options] INPUT.fits"
    parser = OptionParser(usage)
    parser.add_option("-T", "--title", dest="title", default="",
                      help="Plot title")
    parser.add_option("-L", "--label", dest="label", default="",
                      help="Color bar label")
    parser.add_option("-m", "--min", dest="min", type=float,
                      help="Plot minimum value")
    parser.add_option("-M", "--max", dest="max", type=float,
                      help="Plot maximum value")
    parser.add_option("-t", "--ticks", dest="ticks", default=None,
                      help="Ticks to use in plot")
    parser.add_option("-s", "--scale", dest="scale", type=float,
                      help="Scale the map after input")
    parser.add_option("-b","--batchmode", action="store_true", dest="batchMode",
                      default=False, help="Execute without interaction")
    parser.add_option("--mask", action="store_true", dest="mask",
                      default=False, help="Mask part of the sky")
    parser.add_option("-o", "--output", dest="output", default=None,
                      help="Output image file")
    parser.add_option("-c", "--coords", dest="coords", default=None,
                      help="C=equatorial, G=galactic, E=ecliptic")
    parser.add_option("-z","--nonzero", action="store_true", dest="nonZero",
                      default=False, help="Set minimum to lowest, non-zero value")
    parser.add_option("-g","--geo", action="store_true", dest="geo",
                      default=False, help="Flip to geo coordinates")
                      
    (options, args) = parser.parse_args()
    if len(args) != 1:
        parser.error("incorrect number of arguments")

    # Open the FITS file
    map = H.read_map(args[0])
    if options.scale:
        map *= options.scale
    dMin, dMax = (map[argmin(map)], map[argmax(map)])

    if options.min:
        dMin = options.min
    if options.max:
        dMax = options.max
    # Set minimum to lowest non-zero value
    if options.nonZero:
        mapmin = [x for x in map if x != 0]
        dMin = mapmin[argmin(mapmin)]
    
    print "Min: %f\nMax: %f" % (dMin, dMax)

    # Set up the tick axis
    ticks = []
    if options.ticks == None:
        ticks = [dMin, dMax]
    else:
        tlist = options.ticks.split(",")
        ticks = [float(t.strip()) for t in tlist]

    # Set up the coordinate system
    coords = "C"
    rotation = 180
    #rotation = 0
    if options.coords == "G":
        coords="CG"
        rotation = 0
    elif options.coords == "E":
        coords="CE"
        rotation = 0

    flip="astro"

    if options.geo == True:
        flip = "geo"

    # Mask out unused pixels
    if options.mask:
        npix = map.size
        nside = H.npix2nside(npix)
        for i in range(0, npix):
            delta, alpha = H.pix2ang(nside, i)
            delta = 90*degree - delta
            if not IsInIC40Mask(alpha, delta):
            #if not IsInCombinedMask(alpha, delta):
            #if not IsInMilagroMask(alpha, delta):
                if options.geo:
                   map[i] = inf
                else:    
                   map[i] = H.UNSEEN
            
    # Plot the local coordinate map
    H.mollview(map,
               fig=1, coord=coords, rot=rotation,
               title=options.title,
               min=dMin, max=dMax,
               cbar=False,flip = flip,
               margins=(0.1,0.9,0.1,0.9),
               notext=True)
    H.graticule()

    # Set up the color bar
    fig = figure(1)

    SetupColorBar(fig,
        title=options.label, 
        ticks=ticks,
        coord=coords,
        geo=options.geo)

    if options.output:
        savefig(options.output)

    if not options.batchMode:
        show()

