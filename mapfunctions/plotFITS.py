#!/usr/bin/env python

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.colors import LinearSegmentedColormap
import argparse
from mapFunctions import getMap
import os, re, sys
import colormaps as cmaps
# from plotFunctions import cmap_discretize
from plotFunctions import diverge_map


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


def SetupColorBar(label=None, min=None, max=None, ticks=None, coord="C",
                    fontsize='large', projaxis=hp.projaxes.HpxMollweideAxes):
    """ Create the color bar for a HEALPix figure
    """
    fig = plt.figure(1)
    ax = [axis for axis in fig.get_axes() if type(axis) is projaxis][0]
    shrink = 0.7
    orientation = 'horizontal'
    if projaxis == hp.projaxes.HpxCartesianAxes:
        shrink = 0.65
        orientation = 'vertical'
    #    pad = 0.03
    if projaxis == hp.projaxes.HpxOrthographicAxes:
        shrink = 0.6

    cb = fig.colorbar(ax.get_images()[0], ax=ax,
                orientation=orientation,
                shrink=shrink, aspect=50,
                pad=0.01, fraction=0.1,
                ticks=ticks,
                format=FormatStrFormatter("%g"))

    for l in cb.ax.get_xticklabels():
        l.set_fontsize(fontsize)

    ax0 = cb.ax.yaxis if orientation=='vertical' else cb.ax.xaxis
    ls  = ax0.get_ticklabels()
    ls1 = [l.get_text() for l in ls]
    if float(ls1[0]) != float(min):
        ls1[0] = '%s' % min
    if float(ls1[-1]) != float(max):
        ls1[-1] = '%s' % max
    diffs = [float(ls1[i+1]) - float(ls1[i]) for i in range(len(ls1) - 1)]
    if diffs[0] < diffs[1]/2.:
        ls1[1] = ''
    if diffs[-1] < diffs[-2]/2.:
        ls1[-2] = ''
    ax0.set_ticklabels(ls1)

    cb.set_label(label, size=fontsize)

    if projaxis == hp.projaxes.HpxMollweideAxes:
        if coord == "C":
            ax.annotate("0$^\circ$", xy=(1.8, -0.75), size=fontsize)
            ax.annotate("360$^\circ$", xy=(-1.99, -0.75), size=fontsize)
        else:
            ax.annotate("-180$^\circ$", xy=(1.65, 0.625), size=fontsize)
            ax.annotate("180$^\circ$", xy=(-1.9, 0.625), size=fontsize)


def makeTitle(args, options):
    files = [os.path.splitext(os.path.basename(f))[0] for f in args]
    params = np.transpose([f.split('_') for f in files])
    newparams = ['' for i in params]
    for i in range(len(params)):
        values = sorted(list(set(params[i])))
        if len(values) == 1:
            newparams[i] = values[0]
            continue
        if values[0][:2] in ['IC','IT']:
            newparams[i] = '-'.join(values)
        if 'GeV' in values[0]:
            evals = np.array([re.split('-|GeV', v) for v in values]).flatten()
            evals = sorted([j for j in evals if j != ''])
            eflts = [float(j) for j in evals]
            evals = [j for (k,j) in sorted(zip(eflts,evals))]
            newparams[i] = '%s-%sGeV' % (evals[0], evals[-1])
    title = ' '.join(newparams)
    title += ' %s %02ddeg' % (options.mapName, options.smooth)
    title = title.replace('.', 'pt')
    if options.multi:
        title += ' l%ssub' % options.multi
    if options.polar:
        title += ' polar'
    return title


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
            description='description: Plotting program for healpix maps.')

    parser.add_argument('-f','--filepath', help='Path to healpix fit file to be analyzed')
    parser.add_argument('-m', '--min', dest='min', help='Plot minimum value')
    parser.add_argument('-M', '--max', dest='max', help='Plot maximum value')
    parser.add_argument('-d', '--decmin', dest='decmin', type=float,
            default=-90., help='Minimum declination value (90->-90)')
    parser.add_argument('-D', '--decmax', dest='decmax', type=float,
            default=90., help='Maximum declination value (90->-90)')
    parser.add_argument('-r', '--ramin', dest='ramin', type=float,
            help='Minimum RA value')
    parser.add_argument('-R', '--ramax', dest='ramax', type=float,
            help='Maximum RA value')
    parser.add_argument('--mask', dest='mask', default=False,
            action='store_true', help='Intelligent masking')
    parser.add_argument('-b', '--batchmode', action='store_true', dest='batch',
            default=False, help='Execute without interaction')
    parser.add_argument('-o', '--output', dest='output', action='store_true',
            default=False, help='Output image file')
    parser.add_argument('-n', '--mapName', dest='mapName',
            help='Map type desired (signal, relint, relint_err, data, bg)')
    parser.add_argument('-s', '--scale', dest='scale', type=float,
            help='Scale the map after input')
    parser.add_argument('-S', '--smooth', dest='smooth', type=float, default=0,
            help='Desired smoothing radius (in degrees)')
    parser.add_argument('--stype', dest='stype', default='tophat',
            help='Option for smoothing type [tophat|gauss]')
    parser.add_argument('--swindow', dest='swindow', type=float, default=3,
            help='Option for smoothing window')
    parser.add_argument('--multi', dest='multi', type=int,
            help='Use Multipole subtraction')
    parser.add_argument('--fix_multi', dest='fix_multi',
            default=False, action='store_true',
            help='Fix multipole subtraction to values from cumulative map')
    parser.add_argument('-x', '--threshold', dest='threshold', type=float,
            help='Threshold value for plotting data')
    parser.add_argument('-c', '--coords', dest='coords', default='C',
            help='C=equatorial, G=galactic, E=ecliptic')
    parser.add_argument('--gplane', dest='gplane',
            default=False, action='store_true',
            help='Show the galactic plane')
    parser.add_argument('--half', action="store_true", dest='half',
            default=False, help='Show only bottom half of the sky')
    parser.add_argument('--title', action="store_true", dest="title",
            default=False, help='Show the title on the plot')
    parser.add_argument('--outDir', dest='outDir',
            default='/home/jbourbeau/public_html/figures/',
            help='Option for changing output directory')
    parser.add_argument('--prelim', action='store_true', dest='prelim',
            default=False, help='Indicate plot is preliminary')
    parser.add_argument('--llabel', dest='llabel',
            default=False, help='Optional left label overlay on map')
    parser.add_argument('--rlabel', dest='rlabel',
            default=False, help='Optional right label overlay on map')
    parser.add_argument('--polar', dest='polar',
            default=False, action='store_true',
            help='Polar gnomonic view of map')
    parser.add_argument('--customOut', dest='customOut',
            help='Option for custom output file name')
    parser.add_argument('--ext', dest='ext',
            default='png', help='Output file extension')
    parser.add_argument('-v', '--verbose', action='store_true', dest='verbose',
            default=False, help='Verbose output')
    parser.add_argument('--symmetric', action='store_true', dest='symmetric',
            default=True, help='Using diverging colomap with symmetric bounds (i.e. 0=white)')


    args = parser.parse_args()
    #print(args)

    if len(sys.argv) < 1:
        parser.error('Incorrect number of arguments')

    if not args.mapName:
        args.mapName = raw_input('Choose map type [data|bg|signal|relint]: ')
        print

    # detector = os.path.basename(args[0])[:2]
    detector = os.path.basename(args.filepath)[:2]
    if args.mask:
        args.decmax = -25.
        if detector == 'IT':
            args.decmax = -35.
            # if any([k in re.split('_|\.',args[0]) for k in ['p','h','o','f']]):
            if any([k in re.split('_|\.',args.filepath) for k in ['p','h','o','f']]):
                args.decmax = 37.9 - 90.
        #args.mask = False

    opts = vars(args).copy()
    # map = getMap(*args, **opts)
    map = getMap(*[args.filepath], **opts)

    # Multiply by scale
    if args.scale:
        map[map!=hp.UNSEEN] *= (10**args.scale)

    # Setup the coordinate system
    rot = 180
    if args.coords == 'G':
        args.coords = 'CG'
        rot = 0
    if args.coords == 'E':
        args.coords = 'CE'
        rot = 250

    # Mollweide (default) or Cartesian projection
    proj = 'Mollweide'
    projaxis = hp.projaxes.HpxMollweideAxes
    if args.polar:
        proj='Orthographic'
        projaxis = hp.projaxes.HpxOrthographicAxes
    if args.ramin!=None and args.ramax!=None:
        proj = 'Cartesian'
        projaxis = hp.projaxes.HpxCartesianAxes

    # Find min and max for unmasked pixels
    unmasked = np.array([i for i in map if (i!=hp.UNSEEN and i!=np.inf)])
    min = float(args.min) if args.min else unmasked.min()
    max = float(args.max) if args.max else unmasked.max()
    # Make min-max symmetric for relative intensity maps
    # if args.mapName in ['relint','relint_err','signal','fit']:
    #     if (min < 0.) and (max >= 0.) and (abs(min) <= abs(max)):
    #         min = -max
    #     if (min < 0.) and (max >= 0.) and (abs(min) >= abs(max)):
    #         max = -min
    if opts['symmetric']==True:
        if (min < 0.) and (max >= 0.) and (abs(min) <= abs(max)):
            min = -max
        if (min < 0.) and (max >= 0.) and (abs(min) >= abs(max)):
            max = -min
        # colormap = cmap_discretize(plt.get_cmap('seismic'),np.linspace(min,max,30))
        # colormap = plt.get_cmap('jet')
        # colormap = plt.get_cmap('seismic')
        colormap = plt.get_cmap('RdBu_r')
        # colormap = diverge_map(high=('#A00000'),low=('#3F54C0'))
        # colormap = cmaps.viridis
    else:
        colormap = cmaps.viridis
        # colormap = cmap_discretize(cmaps.viridis,np.linspace(min,max,20))
    if not args.min:
        args.min = '%.2f' % min
    if not args.max:
        args.max = '%.2f' % max

    # Setup colormap with option for threshold
    #colormap = plt.get_cmap('jet')
    # if args.mapName in ['relint','relint_err','signal','fit','single']:
    #     #colormap = plt.get_cmap('seismic')
    #     colormap = cmap_discretize(plt.get_cmap('seismic'),np.linspace(min,max,30))
    #     #colormap = cmap_discretize(plt.get_cmap('coolwarm'),np.linspace(min,max,20))
    # else:
    #     colormap = cmap_discretize(cmaps.viridis,np.linspace(min,max,20))
    #colormap = plt.get_cmap('jet')
    #colormap = plt.get_cmap('coolwarm')
    #colormap = cmaps.viridis
    if args.threshold:
        colormap = SetupAbsThresholdColormap(min, max, args.threshold)
    colormap.set_under('white')
    colormap.set_bad('gray')

    # Automatically generated title
    title = makeTitle([args.filepath], args)

    fontsize='large'
    if args.output:
        args.title = False
        fontsize = 'small'
        if proj == 'Cartesian':
            fontsize = 'medium'

    fig = plt.figure(1)
    pltParams = {'fig':1, 'rot':rot, 'title':title, 'min':min, 'max':max, \
            'cbar':False, 'notext':True, 'coord':args.coords, \
            'cmap':colormap}

    if proj == 'Cartesian':
        lonra = [args.ramin, args.ramax]
        latra = [args.decmin, args.decmax]
        hp.cartview(map, lonra=lonra, latra=latra, **pltParams)
        ax = fig.axes[0]
        ax.axis('on')
        ax.set_xlabel(r'Right Ascension [$^{\circ}$]')
        ax.set_ylabel(r'Declination [$^{\circ}$]')
        xlabels = ['%i' % (rot-xtick) for xtick in ax.get_xticks()]
        ax.set_xticklabels(xlabels)
    elif proj == 'Orthographic':

        pltParams['rot'] = [0,-90,rot]
        hp.orthview(map, half_sky=True, **pltParams)
    else:
        hp.mollview(map, **pltParams)
    hp.graticule(verbose=False)

    if args.gplane:
        theta = np.pi/2 * np.ones(100)
        phi = np.arange(0, 360, 3.6) * np.pi/180.
        #hp.projscatter(theta, phi, coord='G')
        hp.projplot(theta, phi, 'k--', coord='G')
        hp.projplot(np.pi/2, 0, 'k^', coord='G')
        #hp.projscatter(np.pi/2, 0, coord='G')
        #hp.projscatter(np.pi/2, np.pi/2., coord='G')
        #hp.projscatter(np.pi/2, np.pi, coord='G')
        #hp.projscatter(np.pi/2, 3/2.*np.pi, coord='G')

    # Set up the color bar
    labelDict = {'relint':'Relative Intensity'}
    labelDict.update({'signal':'Significance [$\sigma$]'})
    labelDict.update({'data':'Data','bg':'Background'})
    labelDict.update({'relint_err':'Relative Intensity Error'})
    labelDict.update({'fit':'Multipole Fit (Relative Intensity)'})
    labelDict.update({'single':'Multipole Fit (Relative Intensity)'})
    label = labelDict[args.mapName]
    if args.scale:
        label += ' [x 10$^{-%d}$]' % args.scale
    cbarParams = {'label':label, 'min':args.min, 'max':args.max, \
            'coord':args.coords, 'fontsize':fontsize, 'projaxis':projaxis}
    SetupColorBar(**cbarParams)

    # Options for half size map and labels
    w, h = fig.get_size_inches()
    lParams = {'size':'large','color':'white','family':'sans-serif'}
    ax = [axis for axis in fig.get_axes() if type(axis) is projaxis][0]
    if args.polar:
        lParams['color'] = 'black'
        if args.llabel:
            lbl = '%s %s' % (args.llabel[:-3], args.llabel[-3:])
            ax.annotate(lbl, xy=(-.9,.9), **lParams)
        if args.rlabel:
            ax.annotate(args.rlabel, xy=(.6,.9), **lParams)
    if not args.title:
        ax.set_title(title, visible=False)
    if args.half:
        #ax.set_ylim(-1, 0.5)
        ax.set_ylim(-1, 0.005)
        if args.llabel:
            lbl = '%s %s' % (args.llabel[:-3], args.llabel[-3:])
            ax.annotate(lbl, xy=(-1.85,-0.24), **lParams)
        if args.rlabel:
            ax.annotate(args.rlabel, xy=(1.3,-0.24), **lParams)
        #fig.set_size_inches(w, h/2.8, forward=True)
        fig.set_size_inches(w, h/2.8, forward=True)
    if args.prelim:
        ax = [axis for axis in fig.get_axes() if type(axis) is projaxis][0]
        IT = 'IceTop' if detector=='IT' else 'IceCube'
        ax.annotate(IT+' Preliminary', xy=(0.2,-0.24), **lParams)

    if args.output:
        outFile  = args.outDir + title.replace(' ', '_')
        if args.customOut:
            outFile = args.outDir + args.customOut
        plt.savefig(outFile+'.'+args.ext, dpi=300, bbox_inches='tight')

    if not args.batch:
        plt.show()
