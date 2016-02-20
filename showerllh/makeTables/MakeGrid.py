#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import argparse

import myGlobals as my
import simFunctions_IT as simFunctions
#from showerllh.analysis.llhTools import getVertex
from showerllh.analysis.llhtools import getVertex

# Read a gcd file and return the tank positions
def getTanks(outFile):

    from icecube import dataio, dataclasses
    tank_xy = {}

    for config in ['IT59','IT73','IT81']:

        gcd = simFunctions.getGCD(config)

        # Extract Geometry and DetectorStatus
        print 'Loading {}...'.format(gcd)
        i3file = dataio.I3File(gcd)
        finished, geo_check, status_check = False, False, False
        while not finished:
            frame = i3file.pop_frame()
            if frame == None:
                print 'Did not find both Geometry and DetectorStatus in file'
                raise
            if frame.Has('I3Geometry') and not geo_check:
                geo = frame.Get('I3Geometry')
                geo_check = True
            if frame.Has('I3DetectorStatus') and not status_check:
                stat = frame.Get('I3DetectorStatus')
                status_check = True
            if geo_check and status_check:
                finished = True

        # add tanks that are active to tankpositions
        tank_xy[config] = []
        for number, station in geo.stationgeo:
            for tank in station:
                om1, om2 = tank.omkey_list
                dst = stat.dom_status
                if om1 not in dst.keys() or om2 not in dst.keys():
                    continue
                elif dst[om1].pmt_hv > 5e-7 and dst[om2].pmt_hv > 5e-7:
                    tank_xy[config].append((tank.position.x, tank.position.y))

    np.save(outFile, tank_xy)


# Create hex with spacing l, nside spots per side, tilt theta, and center xy0
def getHex(l, nside, theta, xy0):

    # Create a baseline
    x0, y0 = xy0
    xx = np.array([x0 + l*np.cos(theta)*n for n in range(-nside+1, nside)])
    yy = np.array([y0 + l*np.sin(theta)*n for n in range(-nside+1, nside)])

    # Expand into a parallelogram
    dx = l * np.cos(theta + np.pi/3)
    dy = l * np.sin(theta + np.pi/3)
    xx_bins = np.array([xx + i*dx for i in range(-nside+1, nside)])
    yy_bins = np.array([yy + i*dy for i in range(-nside+1, nside)])

    # Refine shape to hexagon
    nmid = 2*nside - 1
    k = np.array([[i+j for j in range(nmid)] for i in range(nmid)])
    spots = ((k >= nside-1) * (k <= 2*nmid-nside-1))

    return zip(xx_bins[spots], yy_bins[spots])


# Create grid for ShowerLLH iterative grid search
def makeGrids(**kwargs):

    steps  = kwargs['steps']
    nhex   = kwargs['nhex']

    l = steps[0]    # grid spacing
    Z = 1947        # uniform plane for depth

    # Find the tank positions
    tankpos = np.load(kwargs['tankFile'])
    #print('tankpos = {}'.format(tankpos))
    tankpos = tankpos.item()
    #print('tankpos = {}'.format(tankpos))
    tank_xy = tankpos[kwargs['config']]
    #print('tank_xy = {}'.format(tank_xy))
    tank_x, tank_y = np.transpose(tank_xy)
    #print('tank_x = {}'.format(tank_x))
    #print('tank_y = {}'.format(tank_y))

    # pick out the positions of the tanks on the bottom edge and fit
    bot = (tank_y < -390)
    x_bot, y_bot = tank_x[bot], tank_y[bot]
    m, b = np.polyfit(x_bot, y_bot, 1)
    theta = np.arctan(m)

    # Create oversized hexagon centered at (0,0)
    big_nside = int((len(x_bot)/2 * 125) / l) * 2
    grid_xy = getHex(l, big_nside, theta, [0,0])

    # Reduce to shape of detector
    grid = {}
    grid[0] = []
    for x, y in grid_xy:
        dmin = np.sqrt((tank_x-x)**2 + (tank_y-y)**2).min()
        if dmin <= kwargs['boundaryDist']:
            grid[0].append((x, y))

    # Recursive function for building nested grids
    def fillGrid(coord, depth, grid):
        hexCoords = getHex(steps[depth], nhex[depth-1], theta, coord)
        try: grid[depth][coord] = hexCoords
        except KeyError:
            grid[depth] = {coord:hexCoords}
        if depth != len(steps)-1:
            for newCoord in hexCoords:
                fillGrid(newCoord, depth+1, grid)

    # Build nested grids
    if len(steps) > 1:
        for coord in grid[0]:
            fillGrid(coord, 1, grid)

    # Vertices (for plotting)
    vertx, verty = getVertex(-50, config=kwargs['config'], tank_xy=tank_xy)
    vertx = np.append(vertx, vertx[0])
    verty = np.append(verty, verty[0])

    # Plotting
    g0x, g0y = np.transpose(grid[0])
    if kwargs['plot']:
        fig, ax = plt.subplots()
        ax.plot(g0x, g0y, 'k.')
        ax.plot(tank_x, tank_y, 'ks')
        ax.plot(vertx, verty)
        #for x, y in grid[0]:
        #    xy_temp = getHex(steps[1], nhex[0], theta, (x,y))
        #    x_temp, y_temp = np.transpose(xy_temp)
        #    ax.plot(x_temp, y_temp, '.')
        ax.plot(0, 0, 'go')
        plt.show()

    np.save(kwargs['outFile'], grid)


if __name__ == "__main__":

    # Global variables setup for path names
    my.setupShowerLLH(verbose=False)
    resourcedir = my.llh_resource

    p = argparse.ArgumentParser(
            description='Makes grid points for use in iterative grid search')

    p.add_argument('-c', '--config', dest='config',
            default='IT81',
            help='Detector configuration (determines grid shape)')
    p.add_argument('-s', '--steps', dest='steps', nargs='*', type=float,
            default=[125, 20, 5],
            help='Step sizes for each iteration of the grid')
    p.add_argument('-b', '--boundaryDist', dest='boundaryDist', type=float,
            help='Distance outside tanks to search (default = 2.33 x l)')
    p.add_argument('-n', '--nhex', dest='nhex', type=int,
            help='Number of points per side for smaller grid hexes')
    p.add_argument('--newGCD', dest='newGCD', action='store_true',
            default=False,
            help='Option to recalculate tank positions from gcd files')
    p.add_argument('-o', '--outFile', dest='outFile',
            help='Output filename')
    p.add_argument('--tankFile', dest='tankFile',
            default=resourcedir+'/tankpos.npy',
            help='File containing tank coordinates')
    p.add_argument('-p', '--plot', dest='plot', action='store_true',
            default=False,
            help='Option to plot grid')

    args = p.parse_args()

    # Set default values if not given
    if not args.nhex:
        args.nhex = [6 for i in range(len(args.steps)-1)]
    if not args.boundaryDist:
        args.boundaryDist = 2.33 * args.steps[0]
    if not args.outFile:
        outDir = resourcedir
        args.outFile = '%s/%s_grid.npy' % (outDir, args.config)

    # Check that step  and nhex sizes are valid
    if not all(i > j for i, j in zip(args.steps, args.steps[1:])):
        print 'Grid steps must be in decreasing size'
        raise
    if len(args.nhex) != len(args.steps)-1:
        print 'Incorrect number of arguments for either nhex or steps'
        raise

    if args.newGCD:
        getTanks(args.tankFile)

    opts = vars(args).copy()
    makeGrids(**opts)



