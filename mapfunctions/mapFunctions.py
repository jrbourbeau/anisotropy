#!/usr/bin/env python

import healpy as hp
import numpy as np
import os, glob, ROOT, time
from numpy import sqrt, pi
from numpy.polynomial import Legendre
import matplotlib.pyplot as plt
from projFunctions import getProjDipole


def maskMap(map, decmin, decmax):

    degree = pi / 180.
    npix  = len(map)
    nside = hp.npix2nside(npix)
    theta, phi = hp.pix2ang(nside, range(npix))

    thetaMin = (90 - decmin) * degree
    thetaMax = (90 - decmax) * degree
    thetaCut = (theta <= thetaMin) * (theta >= thetaMax)

    new_map = np.copy(map)
    new_map[np.logical_not(thetaCut)] = hp.UNSEEN

    return new_map


# Li Ma Significance
def LMSignificance(nData, nBG, alpha, nData_wsq=None, nBG_wsq=None):

    with np.errstate(invalid='ignore', divide='ignore'):

        # Allow for scaling term if weighted square maps necessary
        scale = 1.
        if (nData_wsq != None) or (nBG_wsq != None):
            scale = (nData + nBG/alpha) / (nData_wsq + nBG_wsq/alpha)
        Non  = nData * scale
        Noff = nBG/alpha * scale

        sn = np.sign(nData - nBG)
        sigma = sn * sqrt(2*(Non*np.log(((1+alpha)*Non) / (alpha*(Non+Noff)))
            + Noff * np.log(((1+alpha)*Noff) / (Non+Noff))))

    return sigma


# Return the real Cartesian spherical harmonic for a given l, m
def norm_sphharm(l, m, vx, vy, vz):

    if l==0 and m==0:
        return 1/2. * sqrt(1/pi)
    if l==1 and m==-1:
        return vy * sqrt(3/(4*pi))
    if l==1 and m==0:
        return vz * sqrt(3/(4*pi))
    if l==1 and m==1:
        return vx * sqrt(3/(4*pi))
    if l==2 and m==-2:
        return vx*vy * 1/2. * sqrt(15/pi)
    if l==2 and m==-1:
        return vy*vz * 1/2. * sqrt(15/pi)
    if l==2 and m==0:
        return (3*vz**2 - 1.) * 1/4. * sqrt(5/pi)
    if l==2 and m==1:
        return vx*vz * 1/2. * sqrt(15/pi)
    if l==2 and m==2:
        return (vx**2 - vy**2) * 1/4. * sqrt(15/pi)
    if l==3 and m==-3:
        return (3*vx**2 - vy*vx) * vy * 1/4. * sqrt(35/(2*pi))
    if l==3 and m==-2:
        return vx*vy*vz * 1/2. * sqrt(105/pi)
    if l==3 and m==-1:
        return (5*vz**2 - 1) * vy * 1/4. * sqrt(21/(2*pi))
    if l==3 and m==0:
        return (vz**2 - 1) * vz * 1/4. * sqrt(7/pi)
    if l==3 and m==1:
        return (5*vz**2 - 1) * vx * 1/4. * sqrt(21/(2*pi))
    if l==3 and m==2:
        return (vx**2 - vy**2) * vz * 1/4. * sqrt(105/pi)
    if l==3 and m==3:
        return (vx**2 - 3*vy**2) * vx * 1/4. * sqrt(35/(2*pi))
    raise

# Return the real Cartesian spherical harmonic for a given l, m
def real_sphharm(l, m, vx, vy, vz):

    if l==0 and m==0:
        return 1.
    if l==1 and m==-1:
        return vy
    if l==1 and m==0:
        return vz
    if l==1 and m==1:
        return vx
    if l==2 and m==-2:
        return 2*vx*vy
    if l==2 and m==-1:
        return 2*vy*vz
    if l==2 and m==0:
        return sqrt(1/3.)*(3*vz**2 - 1.)
    if l==2 and m==1:
        return 2*vx*vz
    if l==2 and m==2:
        return (vx**2 - vy**2)
    if l==3 and m==-3:
        return (3*vx**2 - vy*vx) * vy
    if l==3 and m==-2:
        return sqrt(8/3.)*vx*vy*vz
    if l==3 and m==-1:
        return sqrt(3/5.)*(5*vz**2 - 1) * vy
    if l==3 and m==0:
        return sqrt(2/5.)*(vz**2 - 1) * vz
    if l==3 and m==1:
        return sqrt(3/5.)*(5*vz**2 - 1) * vx
    if l==3 and m==2:
        return sqrt(2/3.)*(vx**2 - vy**2) * vz
    if l==3 and m==3:
        return (vx**2 - 3*vy**2) * vx


## Creates dipole and quadrupole fit map
def multifit(l, data, bg, alpha, **kwargs):

    defaults = {'params':False, 'out':False, 'verbose':False, \
                'decmax':-25, 'decmin':-90}
    opts = {k:kwargs[k] for k in kwargs if k in defaults}
    opts.update({k:defaults[k] for k in defaults if k not in opts})

    npix = len(data)
    nside = hp.npix2nside(npix)
    degree = pi / 180.
    minZ = opts['decmin'] * degree
    maxZ = opts['decmax'] * degree
    print('skymap decmin = {}'.format(opts['decmin']))
    print('skymap decmax = {}'.format(opts['decmax']))

    # Useful breakdown of l, m values to be used
    nsph = sum([2*l_i+1 for l_i in range(l+1)])
    lvals = [[l_i]*(2*l_i+1) for l_i in range(l+1)]
    mvals = [[m for m in range(-l_i, l_i+1)] for l_i in range(l+1)]
    lvals = [item for sublist in lvals for item in sublist]
    mvals = [item for sublist in mvals for item in sublist]

    # Calculate relative intensity and variance
    with np.errstate(invalid='ignore', divide='ignore'):
        skymap = (data - bg) / bg
        skymapVar = data * (bg + alpha*data) / (bg**3)
    skymap[np.isnan(skymap)] = 0
    skymapVar[np.isnan(skymapVar)] = np.inf

    # Restrict range to desired zenith angles
    vx, vy, vz = hp.pix2vec(nside, [i for i in range(npix)])
    zcut = (vz >= minZ) * (vz <= maxZ)
    skymap[np.logical_not(zcut)] = 0
    skymapVar[np.logical_not(zcut)] = np.inf
    ndata = zcut.sum()

    # Original chi-squared function
    # Chi-squared function to minimize for fit
    def chi2(npar, derivatives, f, par, internal_flag):
        fit = np.zeros(len(vx))
        for i in range(len(lvals)):
            fit += par[i] * norm_sphharm(lvals[i], mvals[i], vx, vy, vz)
        df = skymap - fit
        f[0] = (df**2 / skymapVar).sum()

    # Chi-squared function to minimize with Tibet dipole
    def chi2Tibet(npar, derivatives, f, par, internal_flag):
        fit = np.zeros(len(vx))
        for i in range(len(lvals)):
            fit += par[i] * norm_sphharm(lvals[i], mvals[i], vx, vy, vz)
        df = skymap - fit

        opts['ramin'] = 0.
        opts['ramax'] = 360.
        opts['nbins'] = 24
        opts['lmax'] = l
        opts['decmin'] = -30.
        opts['decmax'] = 90.
        opts['plot'] = False

        deg2rad = np.pi/180.
        amp, amp_err, phase, phase_err = getProjDipole(fit, **opts)
        amp_tibet = 8.5e-4
        amp_err_tibet = 0.2e-4
        phase_tibet = 21.9*deg2rad
        phase_err_tibet = 1.6*deg2rad

        d1 = amp*np.cos(phase*deg2rad)
        d2 = amp*np.sin(phase*deg2rad)
        d1_tibet = amp_tibet*np.cos(phase_tibet)
        d2_tibet = amp_tibet*np.sin(phase_tibet)
        d1_var = (np.cos(phase_tibet)*amp_err_tibet)**2+(amp_tibet*np.sin(phase_tibet)*phase_err_tibet)**2
        d2_var = (np.sin(phase_tibet)*amp_err_tibet)**2+(amp_tibet*np.cos(phase_tibet)*phase_err_tibet)**2
        f[0] = (df**2 / skymapVar).sum()
        #f[0] = (df**2 / skymapVar).sum()+((amp-8.5e-4)/(0.2e-4))**2 +((phase-21.9)/(1.6))**2
        #f[0] = (df**2 / skymapVar).sum()+(d1-d1_tibet)**2/(d1_var) +(d2-d2_tibet)**2/(d2_var)

    # Setup minimizer
    minimizer = ROOT.TMinuit(nsph)
    minimizer.SetFCN(chi2Tibet)
    #minimizer.SetFCN(chi2)
    error_code = ROOT.Long(0)
    #minimizer.mnexcm("SET PRINTOUT", np.array([-1]), 1, error_code)
    minimizer.mnexcm("SET PRINTOUT", np.array([3]), 1, error_code)

    fitparams = ['Y(%i,%i)' % (lvals[i], mvals[i]) for i in range(nsph)]
    for i, p in enumerate(fitparams):
        #minimizer.mnparm(i, p, 0., 1e-6, -1e6, 1e6, error_code)
        minimizer.mnparm(i, p, 1e-6, 1e-6, -1e6, 1e6, error_code)

    # Iterate MIGRAD up to 1000 times
    minimizer.mnexcm("SET STR", np.array([2]), 1, error_code)
    minimizer.mnexcm("MIGRAD", np.array([1000]), 1, error_code)

    # Extract best fit parameters
    p = getFitParams(minimizer, fitparams, ndata)

    if opts['verbose']:
        outputFit(p, fitparams, 1e4)

    if opts['params']:
        return p

    # Create map with dipole/quadrupole fit
    fitmap = np.zeros(len(vx))
    for i in range(nsph):
        #fitmap += p[fitparams[i]] * real_sphharm(lvals[i],mvals[i],vx,vy,vz)
        fitmap += p[fitparams[i]] * norm_sphharm(lvals[i],mvals[i],vx,vy,vz)

    # Plot relative intensity vs right ascension
    # comment decmin/decmax in this section for Tibet
    #opts['decmin'] = -90.
    #opts['decmax'] = -25.
    #opts['plot'] = True
    deg2rad = np.pi/180.
    amp, amp_err, phase, phase_err = getProjDipole(fitmap, **opts)
    print('\nDipole = {}'.format(amp*10**4))
    print('Phase = {} \n'.format(phase))
    amp_tibet = 8.5e-4
    amp_err_tibet = 0.2e-4
    phase_tibet = 21.9*deg2rad
    phase_err_tibet = 1.6*deg2rad

    d1 = amp*np.cos(phase*deg2rad)
    d2 = amp*np.sin(phase*deg2rad)
    d1_tibet = amp_tibet*np.cos(phase_tibet)
    d2_tibet = amp_tibet*np.sin(phase_tibet)
    d1_var = (np.cos(phase_tibet)*amp_err_tibet)**2+(amp_tibet*np.sin(phase_tibet)*phase_err_tibet)**2
    d2_var = (np.sin(phase_tibet)*amp_err_tibet)**2+(amp_tibet*np.cos(phase_tibet)*phase_err_tibet)**2

    df = skymap - fitmap
    chi2_IC_term = (df**2 / skymapVar).sum()
    chi2_tibet_term = (d1-d1_tibet)**2/(d1_var) +(d2-d2_tibet)**2/(d2_var)
    print('\nchi2_IC_term = {}'.format(chi2_IC_term))
    print('chi2_tibet_term = {} \n'.format(chi2_tibet_term))

    if opts['out']:
        hp.write_map(out, fitmap)
        return

    return fitmap


## Get dipole/quadrupole fit parameters in dictionary
def getFitParams(minimizer, fitparams, ndata):

    p = {}
    amin, edm, errdef = [ROOT.Double(0) for i in range(3)]
    nvpar, nparx, ierr = [ROOT.Long(0) for i in range(3)]
    minimizer.mnstat(amin, edm, errdef, nvpar, nparx, ierr)
    ndof = ROOT.Long(ndata - nvpar)
    prob = ROOT.TMath.Prob(amin, ndof)
    p['chi2'] = amin
    p['ndof'] = ndof
    p['prob'] = prob

    for i, par in enumerate(fitparams):
        p[par], p['d'+par] = ROOT.Double(0), ROOT.Double(0)
        minimizer.GetParameter(i, p[par], p['d'+par])

    return p


## Ouput parameters from dipole/quadrupole fit
def outputFit(p, fitparams, scale):

    print '\nchi2/ndf = %.1f / %d = %.2e' % (p['chi2'], p['ndof'], p['prob'])
    print 'Fit values (x %d):' % scale
    for par in fitparams:
        print ' %s = %.3f +/- %.3f' % (par, scale*p[par], scale*p['d'+par])
    #for par in fitparams:
    #    print ' %.3f %.3f' % (scale*p[par], scale*p['d'+par])


# Takes in data, bg & dipole-quadrupole maps and creates false data map
def multi_subtraction(l, data, bg, alpha, **kwargs):

    defaults = {'decmin':-25., 'decmax':-90., 'verbose':False,
                'fix_multi':False, 'fix_data':None, 'fix_bg':None}
    opts = {k:kwargs[k] for k in kwargs if k in defaults}
    opts.update({k:defaults[k] for k in defaults if k not in opts})

    # Create residual map
    with np.errstate(invalid='ignore'):
        relint = (data - bg) / bg
    relint[np.isnan(relint)] = 0.

    if opts['fix_multi']:
        fit = multifit(l, opts['fix_data'], opts['fix_bg'], alpha, **opts)
    else:
        fit = multifit(l, data, bg, alpha, **opts)
    residual = relint - fit

    # Create false data map and normalize
    new_data = bg * (residual + 1)
    norm = data.sum() / float(new_data.sum())
    new_data = new_data * norm

    return new_data


# "Top hat" smoothing over a given angle in degrees
def smoothMap(map, wtsqr=False, **opts):

    if wtsqr==True and opts['stype']=='tophat':
        return None

    npix  = len(map)
    nside = hp.npix2nside(npix)
    smooth_rad = opts['smooth'] * pi/180.
    smooth_map = np.zeros(npix)

    if opts['stype'] == 'tophat':
        vec = np.transpose(hp.pix2vec(nside, np.arange(npix)))
        for i in range(npix):
            neighbors = hp.query_disc(nside, vec[i], smooth_rad)
            smooth_map[i] += map[neighbors].sum()

    if opts['stype'] == 'gauss':
        k = 1. / (smooth_rad)**2
        c3 = k / (4*pi * np.sinh(k))
        checkrad = opts['swindow'] * smooth_rad
        vec = np.transpose(hp.pix2vec(nside, np.arange(npix)))
        for i in range(npix):
            neighbors = hp.query_disc(nside, vec[i], checkrad)
            gaussBeam = c3 * np.exp(k * np.dot(vec[i], vec[neighbors].T))
            if wtsqr:
                gaussBeam = gaussBeam**2
            smooth_map[i] = (map[neighbors] * gaussBeam).sum()

    if opts['stype'] == 'double':

        # Fixed values for gaussian parameters predetermined
        degree = pi / 180.
        f1 = 0.25
        sig1 = 0.46
        k1 = 1. / (sig1*degree)**2
        sig2 = 1.47
        k2 = 1. / (sig2*degree)**2
        checkrad = opts['swindow'] * smooth_rad
        vec = np.transpose(hp.pix2vec(nside, np.arange(npix)))
        for i in range(npix):
            neighbors = hp.query_disc(nside, vec[i], checkrad)
            gaussBeam = f1 * k1/(2*pi) * \
                    np.exp(k1*(np.dot(vec[i], vec[neighbors].T) - 1))
            gaussBeam += (1-f1) * k2/(2*pi) * \
                    np.exp(k2*(np.dot(vec[i], vec[neighbors].T) - 1))
            if wtsqr:
                gaussBeam = gaussBeam**2
            smooth_map[i] = (map[neighbors] * gaussBeam).sum()

    if opts['stype'] == 'gauss2':

        degree = pi / 180.
        #lmax = 3*nside - 1
        lmax = 2*nside
        m2 = map.copy()
        #theta, phi = hp.pix2ang(nside, range(npix))
        #m2[theta < pi/2 + 35*degree] = hp.UNSEEN
        alm = hp.map2alm(m2, lmax)
        k = 1. / (smooth_rad)**2
        bl = np.zeros(lmax+1, dtype='double')
        if not wtsqr:
            bl[0] = 1
            bl[1] = 1/np.tanh(k) - 1/k
            for l in range(1, lmax):
                bl[l+1] = -(2*l+1)/k * bl[l] + bl[l-1]
                # NEED TO ASK FIO WHAT'S HAPPENING HERE
                if bl[l+1] <= 0:
                    bl[l+1] = 0
                    bl[l] = 0

        if wtsqr:
            bl[0] = k / (4*pi * np.tanh(k))
            bl[1] = (2*k*np.cosh(2*k) - np.sinh(2*k)) / (16*pi*np.sinh(k)**2)
            for l in range(1, lmax):
                bl[l+1] = -(2*l+1)/(2*k) * bl[l] + bl[l-1]
                if bl[l+1] <= 0:
                    bl[l+1] = 0
                    bl[l] = 0
            bl *= 4*pi / npix

        a = hp.Alm()
        l_idx, m_idx = a.getlm(lmax)
        for l in range(l_idx.max()):
            alm[l_idx==l] *= bl[l]

        smooth_map = hp.alm2map(alm, nside, lmax, verbose=False)

    return smooth_map

# Takes a file consisting of data, bg, and local maps. Returns desired
# final type (relint, signal, etc) with smoothing, masking, and fitting options.
def getMap(*inFiles, **kwargs):

    # Apply default values if arguments not in input dictionary
    defaults = {'mapName':None, 'multi':False, 'smooth':0, \
                'stype':'tophat', 'swindow':3,\
                'verbose':False, 'mask':False, 'decmin':-90., 'decmax':90.,
                'fix_multi':False}
    opts = {k:kwargs[k] for k in kwargs if k in defaults}
    opts.update({k:defaults[k] for k in defaults if k not in opts})
    alpha = 1/20.

    if not os.path.basename(inFiles[0])[:2] in ['IT','IC']:
        print('Based on filename, detector might not be IC\
                or IT. Please double check RA and Dec cuts.')

    # Require mapName input
    if not opts['mapName']:
        raise SystemExit('mapName parameter not given!')

    # Intelligent masking
    if opts['mask']:
        opts['decmax'] = -25.
        opts['decmin'] = -90.
        if os.path.basename(inFiles[0])[:2] == 'IT':
            opts['decmax'] = -35.

    # Warn for multipole fit
    if opts['multi'] and (opts['decmax'] > -25):
        print 'WARN: multipole fitter applied with mask at %s' % opts['decmax']
        print 'Suggested mask value: -25'

    # Option for verbose mode
    if opts['verbose']:
        print 'Input parameters:'
        for key in sorted(opts.keys()):
            print ' --%s: %s' % (key, opts[key])

    if opts['mapName'] == 'single':
        single = hp.read_map(inFiles[0], verbose=False)
    else:
        # Read in (multiple) input files
        data, bg, local = np.sum([hp.read_map(f, range(3), verbose=False)
                for f in inFiles], axis=0)

    # Option for multipole subtraction
    if opts['fix_multi']:
        ftot = '/data/user/fmcnally/anisotropy/maps/merged/IC_24H_sid.fits'
        opts['fix_data'],opts['fix_bg'] = hp.read_map(ftot, range(2), verbose=False)
    if opts['multi']:
        data = multi_subtraction(opts['multi'], data, bg, alpha, **opts)

    # Option for top-hat smoothing radius
    if opts['smooth'] != 0:
        d0 = data.copy()
        data_wsq = smoothMap(data, wtsqr=True, **opts)
        bg_wsq = smoothMap(bg, wtsqr=True, **opts)
        data = smoothMap(data, **opts)
        bg = smoothMap(bg, **opts)

    # Return desired map type
    if opts['mapName'] == 'data':
        map = data
    elif opts['mapName'] == 'bg':
        map = bg
    elif opts['mapName'] == 'signal':
        map = LMSignificance(data, bg, alpha, data_wsq, bg_wsq)
    elif opts['mapName'] == 'relint':
        with np.errstate(invalid='ignore', divide='ignore'):
            map = (data-bg) / bg
    elif opts['mapName'] == 'relerr':
        with np.errstate(invalid='ignore', divide='ignore'):
            map = (data/bg) * sqrt(1/data + alpha/bg)
    elif opts['mapName'] == 'fit':
        #map = multifit(3, data, bg, alpha, **opts)
        #map = multifit(2, data, bg, alpha, **opts)
        map = multifit(1, data, bg, alpha, **opts)
    elif opts['mapName'] == 'single':
        map = single
    else:
        raise SystemExit('Unrecognized mapName: %s' % opts['mapName'])

    # Eliminate nan's and mask
    map[np.isnan(map)] = 0
    # For desplaying full skyp fit maps comment out next line
    #map = maskMap(map, opts['decmin'], opts['decmax'])

    return map
