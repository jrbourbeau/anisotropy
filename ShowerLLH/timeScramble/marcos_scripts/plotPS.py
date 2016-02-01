#!/usr/bin/env python

from numpy import *
from pylab import *
from optparse import OptionParser
from mpl_toolkits.axes_grid.parasite_axes import SubplotHost
import healpy as H

degree = pi / 180.0

def MultipoleToAngle(l):
    """
    Convert multipole moment to angular scale
    """
    return (180. / l)

def IsInIC40Mask(alpha, delta):
    """Return true if coordinates are within the IC40 mask region
    """
    return delta >= (-85*degree) and delta <= (-5*degree)

def TickToLabel(tick):
    """Convert a floating point tick value to a string label
    """
    if abs(tick - int(tick)) > 0:
        return ("%.1f" % tick)
    else:
        return ("%d" % int(tick))

def PSInIC40(map):
    """Hand-calculate the power spectrum of a map in just the IC40 mask
    """
    npix = len(map)
    nside = H.npix2nside(npix)
    avg = 0
    count = 0
    maskMap = zeros(npix, dtype=float)
    for i in range(0, npix):
        delta, alpha = H.pix2ang(nside, i)
        delta = 90*degree - delta
        if not IsInIC40Mask(alpha, delta):
            avg += map[i]
            count += 1
    avg /= count
    for i in range(0, npix):
        delta, alpha = H.pix2ang(nside, i)
        delta = 90*degree - delta
        if not IsInIC40Mask(alpha, delta):
            map[i] -= avg

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
    ax2.axis["top"].set_label("angular scale [$\circ$]")
    ax2.axis["top"].label.set_visible(True)
    ax2.axis["top"].label.set_size(16)
    ax2.axis["top"].major_ticklabels.set_fontsize(14)
    ax2.axis["right"].major_ticks.set_visible(False)
    ax2.axis["right"].major_ticklabels.set_visible(False)

if __name__ == "__main__":
    # Set up command line options
    usage = "usage: %prog [options] INPUT.fits"
    parser = OptionParser(usage)
    parser.add_option("-l", "--lmin", dest="lmin", type=float,
                      help="Minimum multipole")
    parser.add_option("-L", "--lmax", dest="lmax", type=float,
                      help="Maximum multipole")
    parser.add_option("-c", "--clrange", dest="clRange", default=None,
                      help="Power spectrum range")
    parser.add_option("-s", "--scale", dest="scale", type=float,
                      help="Scale the map after input")
    parser.add_option("-b","--batchmode", action="store_true", dest="batchMode",
                      default=False, help="Execute without interaction")
    parser.add_option("-o", "--output", dest="output", default=None,
                      help="Output image file")
    (options, args) = parser.parse_args()
    if len(args) != 1:
        parser.error("incorrect number of arguments")
    if options.clRange:
        if not len(options.clRange.split(",")) == 2:
            parser.error("clrange: 2 values required")

    # Read in the FITS data
    map = H.read_map(args[0])
    dMin, dMax = (map[argmin(map)], map[argmax(map)])
    avg = average(map)
    map -= avg
    print "Min: %f\nMax: %f\nAvg: %f" % (dMin, dMax, avg)

    # Calculate and plot the power spectrum
    if options.scale:
        map *= options.scale
    cl = H.anafast(map, iter=20, use_weights=True)
    l = arange(1, len(cl), dtype=float)
    cl = cl[1:]
    lMin, lMax = (l[argmin(l)], l[argmax(l)])
    clMin, clMax = (cl[argmin(cl)], cl[argmax(cl)])
    print "Calculated spectrum up to l=%f" % lMax
    if options.lmin:
        lMin = options.lmin
    if options.lmax:
        lMax = options.lmax
    if options.clRange:
        clist = options.clRange.split(",")
        clMin, clMax = [float(c.strip()) for c in clist]

    npix = map.size
    nside = H.npix2nside(npix)
    for i in range(0, npix):
        delta, alpha = H.pix2ang(nside, i)
        delta = 90*degree - delta
        if not IsInIC40Mask(alpha, delta):
            map[i] = 0.
    clMask = H.anafast(map, iter=20, use_weights=True)
    clMask = clMask[1:]

    # Create the figure
    fig = figure(figsize=(12,6))
    ax = SubplotHost(fig, 111)
    fig.add_subplot(ax)

    #ax.loglog(l, cl, "r", l, clMask, "b", linewidth=2)
    ax.loglog(l, cl, "r", linewidth=2, label="All Sky")
    ax.loglog(l, clMask, "b", linewidth=2, label="IC40 Mask")
    ax.set_xlim(lMin, lMax)
    ax.set_ylim(clMin, clMax)
    ax.set_xticks([1, 5, 20, 100])
    ax.set_xticklabels([1, 5, 20, 100])
    ax.axis["bottom"].set_label("multipole $\ell$")
    ax.axis["bottom"].label.set_size(16)
    ax.axis["bottom"].major_ticklabels.set_fontsize(14)
    ax.axis["left"].set_label("$\hat{C}_\ell$ [arb. units]")
    ax.axis["left"].label.set_size(16)
    ax.axis["left"].major_ticklabels.set_fontsize(14)
    ax.legend()

    SetupAngleAxis(ax)

    if options.output:
        savefig(options.output)

    if not options.batchMode:
        show()

