#!/usr/bin/env python

from icecube import icetray, dataio, dataclasses
from numpy import *
import sys
import matplotlib.pyplot as plt
sys.path.append('/home/zgriffith/ShowerLLH/analysis')
from useful import getVertex

# Read a gcd file and return the tank positions
def getTanks(gcd):

    # Extract Geometry and DetectorStatus
    i3file = dataio.I3File(gcd)
    finished, geo_check, status_check = False, False, False
    while not finished:
        frame = i3file.pop_frame()
        if frame == None:
            print 'Did not find both Geometry and DetectorStatus in gcd file'
            sys.exit(1)
        if frame.Has('I3Geometry') and not geo_check:
            geo = frame.Get('I3Geometry')
            geo_check = True
        if frame.Has('I3DetectorStatus') and not status_check:
            stat = frame.Get('I3DetectorStatus')
            status_check = True
        if geo_check and status_check:
            finished = True

    # add tanks that are active to tankpositions
    tankpositions = []
    for number, station in geo.stationgeo:
        for tank in station:
            om1, om2 = tank.omkey_list
            dst = stat.dom_status
            if om1 not in dst.keys() or om2 not in dst.keys():
                continue
            elif dst[om1].pmt_hv > 5e-7 and dst[om2].pmt_hv > 5e-7:
                tankpositions.append(tank.position)

    return tankpositions


# Return new positions given a step of length l 
def mover(X, Y, stepX, stepY, theta, l):

    if not (stepX%2 == stepY%2).all():
        print "movements don't agree"
    dx = l*cos(theta + pi/3)
    dy = l*sin(theta + pi/3)
    x_moved = X + stepX*dx
    y_moved = Y + stepY*dy
    return x_moved, y_moved


# Return a hexagon with nside tanks per side separated by length l
def getHex(l, nside, theta, x0, y0):

    # create a baseline
    xx = array([x0 + l*cos(theta)*n for n in range(-nside+1, nside)])
    yy = array([y0 + l*sin(theta)*n for n in range(-nside+1, nside)])
    # expand into a parallelogram
    xx_bins, yy_bins = [],[]
    for i in range(2*nside-1):
        stepX, stepY = arange(-nside+1, nside), arange(-nside+1, nside)
        x1, y1 = mover(xx[i], yy[i], stepX, stepY, theta, l)
        xx_bins.append(x1)
        yy_bins.append(y1)
    xx_bins = asarray(xx_bins)
    yy_bins = asarray(yy_bins)
    # refine shape to hexagon
    nmid = 2*nside - 1
    spots = array([[True for i in range(nmid)] for j in range(nmid)])
    for i in range(nmid):
        for j in range(nmid):
            k = i+j
            if k < nside-1 or k > 2*nmid-nside-1:
                spots[i][j] = False

    return xx_bins[spots], yy_bins[spots]


# Create grid and save in resource directory
def makeGrids(resourcedir, plot=False):

    gcd, binSpots = {},{}
    gcd['IT59'] = '/data/sim/sim-new/downloads/GCD_20_04_10/' + \
                        'GeoCalibDetectorStatus_IC59.55040_official.i3.gz'
    gcd['IT73'] = '/data/sim/sim-new/downloads/GCD_31_08_11/' + \
                        'GeoCalibDetectorStatus_IC79.55380_L2a.i3.gz'
    gcd['IT81'] = '/data/sim/sim-new/downloads/GCD/' + \
                        'GeoCalibDetectorStatus_IC86.55697_V2.i3.gz'

    #for config in ['IT59', 'IT73']:
    for config in ['IT81']:

        l = 125.   # grid spacing
        Z = 1947   # uniform plane for depth

        # Find the tank positions
        print 'Loading', gcd[config]
        tankpositions = getTanks(gcd[config])
        xx = array([pos.x for pos in tankpositions])
        yy = array([pos.y for pos in tankpositions])

        # pick out the positions of the tanks on the bottom edge and fit
        bot = yy<-390
        yy_bot = yy[bot]
        xx_bot = xx[bot]
        m, b = polyfit(xx_bot, yy_bot, 1)
        theta = arctan(m)

        # create a baseline
        xbase = array([xx_bot[0]+l*cos(theta)*n for n in range(-7,7)])
        ybase = array([yy_bot[0]+l*sin(theta)*n for n in range(-7,7)])
        # build a parallelogram off the baseline
        xx_bins, yy_bins = [],[]
        for i in range(len(xbase)):
            stepX, stepY = arange(-2,12), arange(-2,12)
            x1, y1 = mover(xbase[i], ybase[i], stepX, stepY, theta, l)
            xx_bins.append(x1)
            yy_bins.append(y1)
        xx_bins = asarray(xx_bins)
        yy_bins = asarray(yy_bins)

        # Refine grid to mimic detector shape
        spots = array([[True for i in range(14)] for j in range(14)])
        for i in range(len(spots)):
            for j in range(len(spots[i])):
                dmin = sqrt((xx-xx_bins[i][j])**2 + (yy-yy_bins[i][j])**2).min()
                if dmin > 2.33*l:
                    spots[i][j] = False

        xx_final = xx_bins[spots]
        yy_final = yy_bins[spots]

        # formula for number of points in a hex grid with side n
        hexCounts = lambda n: 3*n*(n-1)+1
        # Declare the number of points per side and spacing
        l_mid, n_mid = l/6, 6
        l_fine, n_fine = 5, 6

        # build coarse, middle, and fine grids
        coarse_grid = zeros((len(xx_final), 3))
        middle_grid = zeros((len(xx_final), hexCounts(n_mid), 3))
        fine_grid = zeros((len(xx_final),hexCounts(n_mid),hexCounts(n_fine),3))

        for i in range(len(xx_final)):
            coarse_grid[i][0] = xx_final[i]
            coarse_grid[i][1] = yy_final[i]
            x_mid, y_mid = getHex(l_mid, n_mid, theta, xx_final[i], yy_final[i])
            for j in range(len(x_mid)):
                middle_grid[i][j][0] = x_mid[j]
                middle_grid[i][j][1] = y_mid[j]
                x_fine, y_fine = getHex(l_fine,n_fine,theta,x_mid[j],y_mid[j])
                for k in range(len(x_fine)):
                    fine_grid[i][j][k][0] = x_fine[k]
                    fine_grid[i][j][k][1] = y_fine[k]

        # Give all grids a uniform z position
        coarse_grid[:,2] = Z
        middle_grid[:,2] = Z
        fine_grid[:,2]   = Z

        binSpots['coarse'] = coarse_grid
        binSpots['middle'] = middle_grid
        binSpots['fine'] = fine_grid

        # Vertices (for plotting)
        vertx, verty = getVertex(-50, config=config, tank_x=xx, tank_y=yy)
        vertx = append(vertx, vertx[0])
        verty = append(verty, verty[0])

        # plotting
        if plot:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(xx_final, yy_final, 'k.')
            ax.plot(xx, yy, 'ks')
            ax.plot(vertx, verty)
            for i in range(len(xx_final)):
                x_temp, y_temp = getHex(l/6, 6, theta, xx_final[i], yy_final[i])
                ax.plot(x_temp, y_temp, '.')
            ax.plot(0, 0, 'o')
            plt.show()

        outFile = '%s/%s_BinSpots.npy' % (resourcedir, config)
        save(outFile, binSpots)


if __name__ == "__main__":

    if len(sys.argv) != 2:
        print 'Usage: %s [resourcedir]' % sys.argv[0]
        sys.exit()

    resourcedir = sys.argv[1]
    makeGrids(resourcedir)

