#!/usr/bin/env python

import subprocess, glob, math, os, argparse

import myGlobals as my
from anisotropy.icesim.analysis import readDist as readDist_IC
from anisotropy.topsim.analysis import readDist as readDist_IT
from anisotropy.mapfunctions.energyCuts import getEbins, getEnergyMaps

def medianLabel(config, emin, emax):

    if config[:2] == 'IC':
        median, sigL, sigR = readDist_IC(config, emin, emax)
    if config[:2] == 'IT':
        median, sigL, sigR = readDist_IT(config, emin, emax)
    median = 10**(9+median)
    median = round(median, -int(math.floor(math.log10(median))) + 1)
    if median >= 1e15:
        return '%.1fPeV' % (median / 1e15)
    if median >= 1e12:
        if median >= 1e13:
            return '%iTeV' % (median / 1e12)
        return '%.1fTeV' % (median / 1e12)


if __name__ == "__main__":

    # Global variables setup for path names
    my.setupAnisotropy(verbose=False)
    mapPrefix = my.ani_maps

    p = argparse.ArgumentParser(
            description='Makes all plots for anisotropy paper')
    p.add_argument('-b', '--batch', dest='batch',
            default=False, action='store_true',
            help='Option to not show interactive plots')
    p.add_argument('--largesmall', dest='largesmall',
            default=False, action='store_true',
            help='Large and small scale structure maps for IceCube')
    p.add_argument('--unsmoothed', dest='unsmoothed',
            default=False, action='store_true',
            help='Unsmoothed relative intensity map')
    p.add_argument('--ic59', dest='ic59',
            default=False, action='store_true',
            help='IC59 20deg smoothed map for comparison')
    p.add_argument('--it', dest='it',
            default=False, action='store_true',
            help='IceTop relint & sig maps')
    p.add_argument('--square', dest='square',
            default=False, action='store_true',
            help='Square subset maps for comparison of small-scale features')
    p.add_argument('--ebins', dest='ebins',
            default=False, action='store_true',
            help='IceCube maps binned in energy')
    p.add_argument('--polar', dest='polar',
            default=False, action='store_true',
            help='IceCube maps binned in energy (polar view)')
    p.add_argument('--powerspec', dest='powerspec',
            default=False, action='store_true',
            help='Power spectrum plot')
    p.add_argument('--solar', dest='solar',
            default=False, action='store_true',
            help='Solar dipole maps')
    p.add_argument('--dipole', dest='dipole',
            default=False, action='store_true',
            help='Plot dipole phase as a function of energy')
    p.add_argument('--minima', dest='minima',
            default=False, action='store_true',
            help='Plot phase of absolute minimum as a function of energy')
    p.add_argument('--reco', dest='reco',
            default=False, action='store_true',
            help='Plot table used for IC energy reconstructions')
    p.add_argument('--edist', dest='edist',
            default=False, action='store_true',
            help='Plot true energy distributions for reco energy bins')
    p.add_argument('--edist2', dest='edist2',
            default=False, action='store_true',
            help='Plot true energy distributions for reco energy bins')
    p.add_argument('--proj', dest='proj',
            default=False, action='store_true',
            help='Plot 1D projection for sid + solar')
    p.add_argument('--projerr', dest='projerr',
            default=False, action='store_true',
            help='Plot 1D projection for anti + ext')
    p.add_argument('--projcomp', dest='projcomp',
            default=False, action='store_true',
            help='Plot 1D projection for each year')
    p.add_argument('--smallscale', dest='smallscale',
            default=False, action='store_true',
            help='Relative intensity for each small-scale region over time')
    p.add_argument('--all', dest='all',
            default=False, action='store_true',
            help='Make all plots')
    args = p.parse_args()
    opts = vars(args)
    if args.all:
        opts = {key:True for key in opts}

    argList  = []
    outDir = '/home/jbourbeau/anisotropy/paperplots/'
    if not os.path.isdir(outDir):
        outDir = '/home/jbourbeau/'
    ext = 'png'
    batch = args.batch

    ##=========================================================================
    ## Skymaps

    cmd = '{}/mapfunctions/plotFITS.py'.format(my.ani_home)
    defArgs  = '-o --mask --outDir {} --ext {}'.format(outDir, ext)
    if batch:
        defArgs += ' -b'

    # Large- and small-scale structure
    mapFile = mapPrefix + '/IC_24H_sid.fits'
    #a  = [mapFile+' -n relint -S 5 -s 3 -m -1.5 -M 1.5 --half']
    a  = [mapFile+' -n relint -S 5 -s 3 -m -1.5 -M 1.5 --half --gplane']
    a += [mapFile+' -n relint -S 5 -s 4 -m -5 -M 5 --multi 2 --half --gplane']
    #a += [mapFile+' -n signal -S 5 -m -45 -M 30 --half']
    #a += [mapFile+' -n signal -S 5 -m -9 -M 9 --multi 2 --half']
    a += [mapFile+' -n signal -S 5 -m -45 -M 45 --half --gplane']
    a += [mapFile+' -n signal -S 5 -m -12 -M 12 --multi 2 --half --gplane']
    if opts['largesmall']:
        argList += a

    # Unsmoothed relative intensity
    mapFile = mapPrefix + '/IC_24H_sid.fits'
    a =  ['-f '+mapFile+' -n relint -s 3 --half -m -2 -M 2']
    a += ['--customOut IC_relint_unsmoothed']
    if opts['unsmoothed']:
        argList += a

    # IC59 map
    mapFile  = mapPrefix + '/IC59_24H_sid.fits'
    a = ['-f '+mapFile+' -n relint -S 20 -s 4 -m -2 -M 2 --multi 2 --half']
    if opts['ic59']:
        argList += a

    # IceTop maps
    mapFile = mapPrefix + '/IT_24H_sid_STA8.fits'
    #a  = [mapFile+' -n relint -S 20 -s 3 -m -3 -M 3']
    #a += [mapFile+' -n signal -S 20 -m -9 -M 9']
    label = medianLabel('IT', 8, 100)
    a  = '%s -n relint -S 20 -s 3 -m -3 -M 3 --llabel %s' % (mapFile, label)
    a += ' --rlabel IceTop --half'
    if opts['it']:
        argList += [a]

    # Anisotropy as a function of energy
    a = []
    eBins = getEbins()
    eMins, eMaxs = eBins[:-1], eBins[1:]
    mapFiles = getEnergyMaps('IC')
    for i, f in enumerate(mapFiles):
        mapFile = ' '.join(f)
        label = medianLabel('IC', eMins[i], eMaxs[i])
        tempArgs = '%s -n relint -S 20 -s 3 --llabel %s' % (mapFile, label)
        tempArgs += ' --half'
        minmax = 3 if float(eMins[i]) >= 6.5 else 1
        tempArgs += ' -m -%i -M %i' % (minmax, minmax)
        # Custom naming for merged file
        if len(f) > 1:
            tempArgs += ' --customOut IC_24H_sid_5pt5-6GeV_relint_20deg'
        a += [tempArgs]
    if opts['ebins']:
        argList += a

    # Polar plots
    a = []
    for i, f in enumerate(mapFiles):
        mapFile = ' '.join(f)
        label = medianLabel('IC', eMins[i], eMaxs[i])
        tempArgs = '%s -n relint -S 20 -s 3 --polar' % mapFile
        tempArgs += ' --llabel %s' % label
        minmax = 3 if float(eMins[i]) >= 6.5 else 1
        tempArgs += ' -m -%i -M %i' % (minmax, minmax)
        if len(f) > 1:
            tempArgs += ' --customOut IC_24H_sid_5pt5-6GeV_relint_20deg_polar'
        a += [tempArgs]
    label = medianLabel('IT', 8, 100)
    mapFile = '%s/IT_24H_sid_STA8.fits' % mapPrefix
    a += ['%s -n relint -S 20 -s 3 -m -3 -M 3 --polar --llabel %s --rlabel %s' % \
            (mapFile, label, 'IceTop')]
    if opts['polar']:
        argList += a

    # Solar dipole
    mapFile = mapPrefix + '/IC_24H_solar.fits'
    a  = [mapFile+' -n relint -S 20 -s 4 -m -3 -M 3 --half']
    a += [mapFile+' -n signal -S 20 -m -33 -M 33 --half']
    mapFile  = mapPrefix + '/IT_24H_solar_NotSTA8.fits'
    mapFile += ' %s/IT_24H_solar_STA8.fits' % mapPrefix
    newArg  = mapFile+' -n relint -S 20 -s 4 -m -3 -M 3 --half'
    newArg += ' --customOut IT_24H_solar_relint'
    a += [newArg]
    newArg  = mapFile+' -n signal -S 20 -m -3.5 -M 3.5 --half'
    newArg += ' --customOut IT_24H_solar_signal'
    a += [newArg]
    if opts['solar']:
        argList += a

    # Apply default arguments and plot command function
    argList = ['{} {} {}'.format(cmd, a, defArgs) for a in argList]
    for a in argList:
        a = a.split(' ')
        proc = subprocess.Popen(a)


    ##=========================================================================
    ## Other plots

    # Square comparison maps
    if opts['square']:
        cmd = '%s/mapFunctions/plotFITS.py' % my.ani_home
        mapFile = mapPrefix + '/IC_24H_sid.fits'
        a  = '%s %s -n signal -S 5 -m -12 -M 12 --multi 2' % (cmd, mapFile)
        a += ' --ramin -90 --ramax 0 --decmin -80 --decmax -35'
        a += ' -o --mask --outDir %s --ext %s' % (outDir, ext)
        a += ' --customOut IC_24H_square'
        if batch:
            a += ' -b'
        proc = subprocess.Popen(a.split(' '))
        # Recreate IC59 original map
        b = a.replace('IC_', 'IC59_')
        b = b.replace('-S 5', '-S 20')
        proc = subprocess.Popen(b.split(' '))

    # Power spectrum
    out = outDir + 'IC_Power_Spectrum.' + ext
    if opts['powerspec']:
        mapFile = '%s/maps/merged/IC_24H_sid.fits' % my.ani_data
        cmd = 'python %s/polspice/spice.py' % my.ani_home
        a = '%s %s -a 1 -A 130 -t 140 --multi 2 -o %s' % (cmd, mapFile, out)
        #a += ' --nofull'
        a = a.split(' ')
        if batch:
            a += ['-b']
        proc = subprocess.Popen(a)

    # Dipole phase
    #out = outDir + 'IC_Dipole_Phase.' + ext
    #if opts['dipole']:
    #    from anisotropy.mapFunctions.dipolePlot import energyPlot
    #    offset = 90
    #    energyPlot(offset=offset, out=out, batch=batch)
    #out = outDir + 'IC_Minimum_Phase.' + ext
    #if opts['minima']:
    #    from anisotropy.mapFunctions.dipolePlot import minimumPlot
    #    offset=90
    #    minimumPlot(offset=offset, out=out, batch=batch)
    if opts['dipole']:
        out = '{}/IC_Dipole.{}'.format(outDir, ext)
        cmd = 'python {}/mapfunctions/phasePlot.py'.format(my.ani_home)
        a = '{} -f energy -l 3 -o {} --offset 90 -n 72'.format(cmd, out)
        #a = '%s -f energy -l 4 -o %s --offset 90' % (cmd, out)
        a = a.split(' ')
        if batch:
            a += ['-b']
        proc = subprocess.Popen(a)


    # Reconstructed energy table
    config = 'IC59'
    out = outDir + config+'_Median_Energy.' + ext
    if opts['reco']:
        from anisotropy.icesim.plots import reco_energy_plot
        reco_energy_plot(config, out=out, batch=batch)

    # Energy distribution plot
    config = 'IC86'
    out = outDir + config+'_Energy_Distributions.' + ext
    if opts['edist']:
        from anisotropy.icesim.plots import eres
        eres(config, out=out, batch=batch)

    config = 'IC59'
    out = outDir + config+'_Energy_Dist2.' + ext
    if opts['edist2']:
        from anisotropy.icesim.plots import eres2
        eres2(config, out=out, batch=batch)

    # One-dimensional projection
    out = outDir + 'IC_proj1d.' + ext
    mapTypes = ['sid','solar']
    #mapTypes = ['anti','ext']
    if opts['proj']:
        mapFiles = ['%s/IC_24H_%s.fits' % (mapPrefix, i) \
                for i in mapTypes]
        a  = ['%s/mapFunctions/proj1d.py' % my.ani_home]
        a += mapFiles
        a += ['-z','-o',out,'--syserr','--labels','method']
        #a += ['-z','-o',out]
        if batch:
            a += ['-b']
        proc = subprocess.Popen(a)

    # One-dimensional projection
    out = outDir + 'IC_proj1d_err.' + ext
    mapTypes = ['anti','ext']
    if opts['projerr']:
        mapFiles = ['%s/IC_24H_%s.fits' % (mapPrefix, i) \
                for i in mapTypes]
        a  = ['%s/mapFunctions/proj1d.py' % my.ani_home]
        a += mapFiles
        a += ['-z','-o',out,'--labels','method']
        if batch:
            a += ['-b']
        proc = subprocess.Popen(a)

    # One-dimensional projection for each year
    if opts['projcomp']:
        cmd = '%s/mapFunctions/proj1d.py' % my.ani_home
        configs = ['IC59','IC79','IC86','IC86-II','IC86-III','IC86-IV']
        files = ['%s/%s_24H_sid.fits' % (mapPrefix, cfg) for cfg in configs]
        files = ' '.join(files)
        out = outDir + 'IC_proj1d_comp.' + ext
        a = '%s %s -z -o %s --syserr --labels configs' % (cmd, files, out)
        if batch:
            a += ' -b'
        proc = subprocess.Popen(a.split(' '))

    if opts['smallscale']:
        cmd = 'python %s/mapFunctions/testsmall.py' % my.ani_home
        out = outDir + 'IC_SmallScale.' + ext
        a = '%s --paper -o %s' % (cmd, out)
        if batch:
            a += ' -b'
        proc = subprocess.Popen(a.split(' '))



