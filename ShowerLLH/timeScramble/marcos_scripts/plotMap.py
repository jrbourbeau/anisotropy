#!/usr/bin/env python

from pylab import *
from optparse import OptionParser
import healpy as H
from matplotlib.ticker import FormatStrFormatter
from matplotlib.colors import LinearSegmentedColormap
from scipy import interpolate

degree = pi / 180.0

def IsInIC40Mask(alpha, delta):
    """ Return true if coordinates are within the IC40 mask region 
    """
    return delta >= (-85*degree) and delta <= (-5*degree)

def IsInIceTopMask(alpha, delta):
    """ Return true if coordinates are within the IC40 mask region 
    """
    return delta >= (-85*degree) and delta <= (-35*degree)

def IsInIC59Mask(alpha, delta):
    """ Return true if coordinates are within the IC59 mask region 
    """
    return delta >= (-88*degree) and delta <= (-25*degree)

def IsInMilagroMask(alpha, delta):
    """ Return true if coordinates are within the Milagro mask region
    """
    return delta >= (-6*degree) and delta <= (63*degree)

def IsInCombinedMask(alpha, delta):
    """ Return true if coordinates are within the Milagro mask region
    """
    return (delta >= (-6*degree) and delta <= (63*degree)) or (delta >= -88 * degree and delta <= -25 * degree)

def SetupColorBar(fig, title=None, ticks=None, coord="C", flip="astro"):
    """ Create the color bar for a HEALPix figure
    """
    fig = figure(1)
    for ax in fig.get_axes():
        if type(ax) is H.projaxes.HpxMollweideAxes:
            cb = fig.colorbar(ax.get_images()[0], ax=ax,
                              orientation="horizontal",
                              shrink=0.8, aspect=50,
                              pad=0.03, fraction=0.1,
                              ticks=ticks,
                              format=FormatStrFormatter("%g"))
            for label in cb.ax.get_xticklabels():
                #label.set_fontsize("large")
                label.set_fontsize("xx-large")
            cb.set_label(title, size="xx-large")
            #cb.set_label(title, size="large")
            if coord == "C":
                #ax.annotate("0h", xy=(1.8, 0.625), size="large")
                #ax.annotate("24h", xy=(-1.95, 0.625), size="large")
                if flip == "astro":
                   #ax.annotate("0$^\circ$", xy=(1.8, 0.625), size="x-large")
                   #ax.annotate("360$^\circ$", xy=(-1.95, 0.625), size="x-large")
                   ax.annotate("0$^\circ$", xy=(1.8, -0.75), size="xx-large")
                   ax.annotate("360$^\circ$", xy=(-1.95, -0.75), size="xx-large")
                else:
                   ax.annotate("360$^\circ$", xy=(1.8, 0.625), size="x-large")
                   ax.annotate("0$^\circ$", xy=(-1.95, 0.625), size="x-large")
                #ax.annotate("0$^\circ$", xy=(1.8, -0.625))
                #ax.annotate("360$^\circ$", xy=(-1.95, -0.625))
            else:
                if flip == "astro":
                   ax.annotate("-180$^\circ$", xy=(1.65, 0.625), size="large")
                   ax.annotate("180$^\circ$", xy=(-1.9, 0.625), size="large")
                else:
                   ax.annotate("180$^\circ$", xy=(1.65, 0.625), size="large")
                   ax.annotate("-180$^\circ$", xy=(-1.9, 0.625), size="large")


def cmap_discretize(cmap, N):
    """Return a discrete colormap from the continuous colormap cmap.
    
        cmap: colormap instance, eg. cm.jet. 
        N: Number of colors.
    
    Example
        x = resize(arange(100), (5,100))
        djet = cmap_discretize(cm.jet, 5)
        imshow(x, cmap=djet)
    """

    cdict = cmap._segmentdata.copy()
    # N colors
    colors_i = linspace(0,1.,N)
    # N+1 indices
    indices = linspace(0,1.,N+1)
    for key in ('red','green','blue'):
        # Find the N colors
        D = array(cdict[key])
        I = interpolate.interp1d(D[:,0], D[:,1])
        colors = I(colors_i)
        # Place these colors at the correct indices.
        A = zeros((N+1,3), float)
        A[:,0] = indices
        A[1:,1] = colors
        A[:-1,2] = colors
        # Create a tuple for the dictionary.
        L = []
        for l in A:
            L.append(tuple(l))
        cdict[key] = tuple(L)
    # Return colormap object.
    return matplotlib.colors.LinearSegmentedColormap('colormap',cdict,1024)


def SetupThresholdColormap(amin, amax, threshold):
    """ Create a color map that draws all values below the threshold in
        grayscale, and everything above in the usual "jet" rainbow RGB scale
    """
    thresh = (threshold - amin) / (amax - amin)
    if threshold <= amin or threshold >= amax:
        thresh = 0.
#    dthresh = 1 - thresh
#    threshDict = { "blue"  : ((0.0, 1.0, 1.0), (thresh, 0.4, 0.5),
#                             (thresh+0.11*dthresh,  1, 1),
#                             (thresh+0.34*dthresh, 1, 1),
#                             (thresh+0.65*dthresh, 0, 0),
#                             (1, 0, 0)),
#                   "green" : ((0.0, 1.0, 1.0), (thresh, 0.4, 0.0),
#                              (thresh+0.125*dthresh, 0, 0),
#                              (thresh+0.375*dthresh, 1, 1),
#                              (thresh+0.64*dthresh, 1, 1),
#                              (thresh+0.91*dthresh, 0, 0),
#                              (1, 0, 0)),
#                   "red"   : ((0.0, 1.0, 1.0), (thresh, 0.4, 0.0),
#                              (thresh+0.35*dthresh, 0, 0),
#                              (thresh+0.66*dthresh, 1, 1),
#                              (thresh+0.89*dthresh, 1, 1),
#                              (1, 0.5, 0.5)) }
    threshDict = { "red"   : ((0.0, 1.0, 1.0),
                              (thresh, 0.4, 0.7),
                              (1.0, 1.0, 1.0)),
                   "green" : ((0.0, 1.0, 1.0),
                              (thresh, 0.4, 0.0),
                              (1.0, 1.0, 1.0)),
                   "blue"  : ((0.0, 1.0, 1.0),
                              (thresh, 0.4, 0.0),
                              (1.0, 0.5, 1.0)) }
    return LinearSegmentedColormap("thresholdColormap", threshDict, 256)

def SetupAbsThresholdColormap(amin, amax, threshold):
    """ Create a color map for "two-sided" thresholds.  Below the threshold,
        the map is a cool green-blue palette.  Between the lower and upper
        threshold, the map is gray-white-gray.  Above the upper threshold,
        the map is a warm red-yellow palette.
    """
    x1 = (-threshold - amin) / (amax - amin)
    x3 = (amax - threshold) / (amax - amin)
    x2 = 1. - x1 - x3
    gvl = 0.5
    thrDict = {
        "red"    : ((0.0, 1.0, 0.5), (x1, 0.0, gvl), (x1 + 0.5*x2, 1.0, 1.0),
                    (x1 + x2, gvl, 0.7), (1.0, 1.0, 1.0)),
        "green"  : ((0.0, 1.0, 1.0), (x1, 0.0, gvl), (x1 + 0.5*x2, 1.0, 1.0),
                    (x1 + x2, gvl, 0.0), (1.0, 1.0, 1.0)),
        "blue"   : ((0.0, 1.0, 1.0), (x1, 0.7, gvl), (x1 + 0.5*x2, 1.0, 1.0),
                    (x1 + x2, gvl, 0.0), (1.0, 0.5, 1.0)) }
    return LinearSegmentedColormap("thresholdColormap", thrDict, 256)

def SetupBVRYColormap():
    cdict = { 'red':   ((0.0,  0.0, 0.0),
                    (0.25, 0.0, 0.0),
                    (0.57, 1.0, 1.0),
                    (1.0,  1.0, 1.0)),

          'green': ((0.0,  0.0, 0.0),
                    (0.42, 0.0, 0.0),
                    (0.92, 1.0, 1.0),
                    (1.0,  1.0, 1.0)),

          'blue':  ((0.0,  0.0, 0.0),
                    (0.25, 1.0, 1.0),
                    (0.42, 1.0, 1.0),
                    (0.92, 0.0, 0.0),
                    (1.0,  1.0, 1.0))
        }
    return LinearSegmentedColormap('BVRY', cdict)

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
    parser.add_option("--ncols", dest="ncols", type=float,
                      help="Number of colors in palette")
    parser.add_option("-l", "--logz", action="store_true", dest="logz",
                      default=False, help="Plot z-axis on a log scale")
    parser.add_option("-b","--batchmode", action="store_true", dest="batchMode",
                      default=False, help="Execute without interaction")
    #parser.add_option("--mask", action="store_true", dest="mask",
    #                  default=False, help="Mask part of the sky")
    parser.add_option("--mask", dest="mask",
                      default="X", help="I=IceCube, M=Milagro, C=combined")
    parser.add_option("-S", "--smooth", dest="sigma", type=float,
                      help="Smooth with Gaussian of given sigma, in degrees")
    parser.add_option("-o", "--output", dest="output", default=None,
                      help="Output image file")
    parser.add_option("-c", "--coords", dest="coords", default=None,
                      help="C=equatorial, G=galactic, E=ecliptic")
    parser.add_option("--grid", dest="grid", default='C',
                      help="C=equatorial, G=galactic, E=ecliptic")
    parser.add_option("-x", "--threshold", dest="threshold", type=float,
                      help="Threshold value for plotting data")
    parser.add_option("-y", "--bvry", action="store_true",dest="bvry", default=False,
                      help="BVRY color palette")
    parser.add_option("-z","--nonzero", action="store_true", dest="nonZero",
                      default=False, help="Set minimum to lowest, non-zero value")
    parser.add_option("-w","--whitebkg", action="store_true", dest="whiteBkg",
                      default=False, help="Set bkg color to white")
    parser.add_option("-g","--geo", action="store_true", dest="geo",
                      default=False, help="Flip coordinates to geo")
    parser.add_option("--half", action="store_true", dest="half",
                      default=False, help="Show only bottom half of the sky")
                      
    (options, args) = parser.parse_args()
    if len(args) != 1:
        parser.error("incorrect number of arguments")

    rcParams.update({ "font.family" : "serif" })

    # Open the FITS file
    map = H.read_map(args[0])
    if options.scale:
        map *= options.scale

    # Plot in log scale
    if options.logz:
        npix = map.size
        for i in range(0, npix):
            if map[i] > 0:
                map[i] = log10(map[i])
            else:
                if options.whiteBkg:
                   map[i] = inf
                else: 
                   map[i] = H.UNSEEN

    # Mask out unused pixels
    if options.mask != "X":
        npix = map.size
        nside = H.npix2nside(npix)
        for i in range(0, npix):
            delta, alpha = H.pix2ang(nside, i)
            delta = 90*degree - delta
            if options.mask == "I":
               if not IsInIC59Mask(alpha, delta):
                  if options.whiteBkg:
                     map[i] = inf
                  else:
                     map[i] = H.UNSEEN

            if options.mask == "M":
               if not IsInMilagroMask(alpha, delta):
                  if options.whiteBkg:
                     map[i] = inf
                  else:
                     map[i] = H.UNSEEN

            if options.mask == "C":
               if not IsInCombinedMask(alpha, delta):
                  if options.whiteBkg:
                     map[i] = inf
                  else:
                     map[i] = H.UNSEEN
            
            if options.mask == "T":
               if not IsInIceTopMask(alpha, delta):
                  if options.whiteBkg:
                     map[i] = inf
                  else:
                     map[i] = H.UNSEEN

    # Find the minimum and maximum of the map for pixels not masked out
    notMasked = [x for x in map if (x > H.UNSEEN and x != inf)]
    dMin, dMax = (notMasked[argmin(notMasked)], notMasked[argmax(notMasked)])

#    if options.threshold:
#        #map[abs(map) < options.threshold] = H.UNSEEN
#        map[map < options.threshold] = H.UNSEEN
#        #dMin = options.threshold

    # Smooth if desired
    if options.sigma:
        map = H.smoothing(map, sigma=options.sigma, degree=True)

    if options.min:
        dMin = options.min
    if options.max:
        dMax = options.max
    # Set minimum to lowest non-zero value
    if options.nonZero:
        mapmin = [x for x in map if x != 0]
        dMin = mapmin[argmin(mapmin)]
    
    print "Min: %f\nMax: %f" % (dMin, dMax)

    colormap = cm.jet
    #get_cmap("gist_heat")

    if options.ncols != None:
       colormap = cmap_discretize(colormap,options.ncols)
    if options.threshold:
        colormap = SetupAbsThresholdColormap(dMin, dMax, options.threshold)
    if options.bvry:
        colormap = SetupBVRYColormap()

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
    if options.coords == "G":
        coords="CG"
        rotation = 0
    elif options.coords == "E":
        coords="CE"
        rotation = 250

    # Flip the coordinate system to geo if wanted

    flip = "astro"
    if options.geo:
       flip = "geo"

    # Setting figure sizes for plot
    w = 9 
    h = 6
    extent = (0.02,0.05,0.96,0.95)

    if options.half:
       h /= 2.2
       extent = (0.02,0.11,0.96,1.1)

    fig = figure(1,figsize=(w,h))  

    # Plot the local coordinate map
    H.mollview(map, fig=1, coord=coords, rot=rotation,
               title=options.title,
               min=dMin, max=dMax,
               cbar=False,
               notext=True,
               flip = flip,
               cmap=colormap,
               extent=extent
               #cmap=cm.RdYlBu_r,
               )
    
    # Draw a grid
    H.graticule(coord=options.grid)
    #H.graticule(coord='G',dpar=360,dmer=360)

    if options.half:
       for ax in fig.get_axes():
          if type(ax) is H.projaxes.HpxMollweideAxes:
             ax.set_ylim(-1,0.005)

    # Set up the color bar
    SetupColorBar(fig,
        title=options.label, 
        ticks=ticks,
        coord=coords,
        flip = flip)

    defaultSize = fig.get_size_inches()
    fig.set_size_inches(defaultSize[0], defaultSize[1] + 1.05)
    #fig.savefig("image.png")

    if options.output:
        if options.output.find("jpg") > 0:
            fig.savefig(options.output, dpi=300)
        else:
            savefig(options.output)
            jpgOutput = options.output.replace("png", "jpg")
            fig.savefig(jpgOutput, dpi=300)

    if not options.batchMode:
        show()

