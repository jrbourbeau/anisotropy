#!/usr/bin/env python

#=========================================================================
# File Name     : multipolefit.py
# Description   : Perform spherical harmonic multipole fitting
# Creation Date : 05-06-2016
# Last Modified : Fri 06 May 2016 12:16:40 PM CDT
# Created By    : James Bourbeau
#=========================================================================

import os
import glob
import ROOT
import time
import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
from numpy import sqrt, pi
from projFunctions import getProjDipole, getRIRAProj, getHarmonicFitParams

# fit SH to skymap
def multifit(l, data, bg, alpha, **kwargs):

    defaults = {'params': False, 'out': False, 'verbose': False,
                'decmax': -25, 'decmin': -90, 'chi2': None,
                'init': 1e-15, 'step': 1e-06}
    # opts = {k:kwargs[k] for k in kwargs if k in defaults}
    opts = kwargs.copy()
    opts.update({k: defaults[k] for k in defaults if k not in opts})
    print('opts = {}'.format(opts))

    # Useful breakdown of l, m values to be used
    nsph = sum([2 * l_i + 1 for l_i in range(l + 1)])
    lvals = [[l_i] * (2 * l_i + 1) for l_i in range(l + 1)]
    mvals = [[m for m in range(-l_i, l_i + 1)] for l_i in range(l + 1)]
    lvals = [item for sublist in lvals for item in sublist]
    mvals = [item for sublist in mvals for item in sublist]

    # Calculate relative intensity and variance
    with np.errstate(invalid='ignore', divide='ignore'):
        skymap = (data - bg) / bg
        skymapVar = data * (bg + alpha * data) / (bg**3)
    skymap[np.isnan(skymap)] = 0
    skymapVar[np.isnan(skymapVar)] = np.inf

    # Cut to desired dec range (equiv to healpy zenith range)
    deg2rad = np.pi / 180
    npix = len(data)
    nside = hp.npix2nside(npix)
    theta, phi = hp.pix2ang(nside, range(npix))
    thetamax = (90 - opts['decmin']) * deg2rad
    thetamin = (90 - opts['decmax']) * deg2rad
    pass_dec_cut = (theta <= thetamax) * (theta >= thetamin)
    skymap[np.logical_not(pass_dec_cut)] = 0
    skymapVar[np.logical_not(pass_dec_cut)] = np.inf
    ndata = pass_dec_cut.sum()

    vx, vy, vz = hp.pix2vec(nside, [i for i in range(npix)])
    normedSH = np.load('normedSH.npy')
    normedSH = normedSH.item()
    fitparams = ['Y(%i,%i)' % (lvals[i], mvals[i]) for i in range(nsph)]
    fiterrparams = ['dY(%i,%i)' % (lvals[i], mvals[i]) for i in range(nsph)]

    def chi2(npar, derivatives, f, par, internal_flag):
        fit = np.zeros(len(vx))
        for i in range(len(lvals)):
            #fit += par[i] * norm_sphharm(lvals[i], mvals[i], vx, vy, vz)
            fit += par[i] * normedSH[fitparams[i]]
        df = skymap - fit
        f[0] = (df**2 / skymapVar).sum()

    # Setup minimizer
    minimizer = ROOT.TMinuit(nsph)
    minimizer.SetFCN(chi2)
    error_code = ROOT.Long(0)
    minimizer.mnexcm("SET PRINTOUT", np.array([-1]), 1, error_code)
    # minimizer.mnexcm("SET PRINTOUT", np.array([3]), 1, error_code)

    step = opts['step']
    init = opts['init']
    for i, p in enumerate(fitparams):
        minimizer.mnparm(i, p, init, step, -1e6, 1e6, error_code)
        # minimizer.mnparm(i, p, 0., 1e-6, -1e6, 1e6, error_code)
        # minimizer.mnparm(i, p, 1e-12, 1e-6, -1e6, 1e6, error_code)

    # Iterate MIGRAD up to 1000 times
    minimizer.mnexcm("SET STR", np.array([2]), 1, error_code)
    minimizer.mnexcm("MIGRAD", np.array([2000]), 1, error_code)

    # Extract best fit parameters
    p = getFitParams(minimizer, fitparams, ndata)
    #np.save('/home/jbourbeau/anisotropy/dev/typedSHparms_lmax_2.npy', p)

    if opts['verbose']:
        outputFit(p, fitparams, 1e4)

    if opts['params']:
        return p

    if opts['chi2'] != None:
        return fitparams, fiterrparams, p

    # Create map with dipole/quadrupole fit
    fitmap = np.zeros(len(vx))
    fiterrmap = np.zeros(len(vx))
    for i in range(nsph):
        fitmap += p[fitparams[i]] * normedSH[fitparams[i]]
        fiterrmap += p[fiterrparams[i]] * normedSH[fitparams[i]]

    if opts['out']:
        hp.write_map(out, fitmap)
        return

    return fitmap

# # fit SH to skymap
# def multifit(l, data, bg, alpha, **kwargs):
#
#     defaults = {'params':False, 'out':False, 'verbose':False, \
#                 'decmax':-25, 'decmin':-90, 'chi2': None}
#     # opts = {k:kwargs[k] for k in kwargs if k in defaults}
#     opts = kwargs.copy()
#     opts.update({k:defaults[k] for k in defaults if k not in opts})
#     print('opts = {}'.format(opts))
#
#     # Useful breakdown of l, m values to be used
#     nsph = sum([2*l_i+1 for l_i in range(l+1)])
#     lvals = [[l_i]*(2*l_i+1) for l_i in range(l+1)]
#     mvals = [[m for m in range(-l_i, l_i+1)] for l_i in range(l+1)]
#     lvals = [item for sublist in lvals for item in sublist]
#     mvals = [item for sublist in mvals for item in sublist]
#
#     # Calculate relative intensity and variance
#     with np.errstate(invalid='ignore', divide='ignore'):
#         skymap = (data - bg) / bg
#         skymapVar = data * (bg + alpha*data) / (bg**3)
#     skymap[np.isnan(skymap)] = 0
#     skymapVar[np.isnan(skymapVar)] = np.inf
#
#     # Cut to desired dec range (equiv to healpy zenith range)
#     deg2rad = np.pi / 180
#     npix  = len(data)
#     nside = hp.npix2nside(npix)
#     theta, phi = hp.pix2ang(nside, range(npix))
#     thetamax = (90-opts['decmin']) * deg2rad
#     thetamin = (90-opts['decmax']) * deg2rad
#     pass_dec_cut = (theta <= thetamax) * (theta >= thetamin)
#     skymap[np.logical_not(pass_dec_cut)] = 0
#     skymapVar[np.logical_not(pass_dec_cut)] = np.inf
#     ndata = pass_dec_cut.sum()
#
#     # Chi-squared function to minimize with Tibet dipole
#     # Tibet measurements for the amp and phase
#     amp_tibet = 8.5e-4
#     amperr_tibet = 0.2e-4
#     phase_tibet = 21.9*deg2rad
#     phase_err_tibet = 1.6*deg2rad
#     d1_tibet = amp_tibet*np.cos(phase_tibet)
#     d2_tibet = amp_tibet*np.sin(phase_tibet)
#     d1_var = (np.cos(phase_tibet)*amperr_tibet)**2+(amp_tibet*np.sin(phase_tibet)*phase_err_tibet)**2
#     d2_var = (np.sin(phase_tibet)*amperr_tibet)**2+(amp_tibet*np.cos(phase_tibet)*phase_err_tibet)**2
#     # Tibet measurements for the proj relint from Markus
#     tibetdata = np.loadtxt('Tibet_data.dat')
#     tibetRA = tibetdata[:,0]
#     tibetRAerr = [10.]*18
#     tibetRI = tibetdata[:,1]-1.
#     tibetRIerr = np.array([2.8e-5]*18)
#
#     opts.update({'ramin':0., 'ramax':360., 'nbins':18, 'lmax':3,
#         'decmin':-30., 'decmax':90., 'plot':False})
#
#     vx, vy, vz = hp.pix2vec(nside, [i for i in range(npix)])
#     normedSH = np.load('normedSH.npy')
#     normedSH = normedSH.item()
#     fitparams = ['Y(%i,%i)' % (lvals[i], mvals[i]) for i in range(nsph)]
#     fiterrparams = ['dY(%i,%i)' % (lvals[i], mvals[i]) for i in range(nsph)]
#     def chi2(npar, derivatives, f, par, internal_flag):
#         fit = np.zeros(len(vx))
#         for i in range(len(lvals)):
#             #fit += par[i] * norm_sphharm(lvals[i], mvals[i], vx, vy, vz)
#             fit += par[i] * normedSH[fitparams[i]]
#         df = skymap - fit
#
#         if opts['chi2']=='standard' or opts['chi2']==None:
#             f[0] = (df**2 / skymapVar).sum()
#         if opts['chi2']=='d1d2':
#             amp, amp_err, phase, phase_err = getProjDipole(fit, **opts)
#             d1 = amp*np.cos(phase*deg2rad)
#             d2 = amp*np.sin(phase*deg2rad)
#             f[0] = (df**2 / skymapVar).sum()+(d1-d1_tibet)**2/(d1_var) +(d2-d2_tibet)**2/(d2_var)
#         if opts['chi2']=='RI':
#             RA, RI, RAerr = getRIRAProj(fit,**opts)
#             f[0] = (df**2 / skymapVar).sum()+(((RI-tibetRI)/tibetRIerr)**2).sum()
#         if opts['chi2']=='tibetonly':
#             RA, RI, RAerr = getRIRAProj(fit,**opts)
#             f[0] = (((RI-tibetRI)/tibetRIerr)**2).sum()
#
#     # Setup minimizer
#     minimizer = ROOT.TMinuit(nsph)
#     minimizer.SetFCN(chi2)
#     error_code = ROOT.Long(0)
#     minimizer.mnexcm("SET PRINTOUT", np.array([-1]), 1, error_code)
#     # minimizer.mnexcm("SET PRINTOUT", np.array([3]), 1, error_code)
#
#     step = opts['step']
#     init = opts['init']
#     for i, p in enumerate(fitparams):
#         minimizer.mnparm(i, p, init, step, -1e6, 1e6, error_code)
#         # minimizer.mnparm(i, p, 0., 1e-6, -1e6, 1e6, error_code)
#         # minimizer.mnparm(i, p, 1e-12, 1e-6, -1e6, 1e6, error_code)
#
#     # Iterate MIGRAD up to 1000 times
#     minimizer.mnexcm("SET STR", np.array([2]), 1, error_code)
#     minimizer.mnexcm("MIGRAD", np.array([2000]), 1, error_code)
#
#     # Extract best fit parameters
#     if opts['chi2']=='RI':
#         ndata += 18
#     if opts['chi2']=='tibetonly':
#         ndata = 18
#     p = getFitParams(minimizer, fitparams, ndata)
#     #print('p = {}'.format(p))
#     #np.save('/home/jbourbeau/anisotropy/dev/typedSHparms_lmax_2.npy', p)
#
#     if opts['verbose']:
#         outputFit(p, fitparams, 1e4)
#
#     if opts['params']:
#         return p
#
#     if opts['chi2']!=None:
#         return fitparams, fiterrparams, p
#
#     # Create map with dipole/quadrupole fit
#     fitmap = np.zeros(len(vx))
#     fiterrmap = np.zeros(len(vx))
#     for i in range(nsph):
#         #fitmap += p[fitparams[i]] * norm_sphharm(lvals[i],mvals[i],vx,vy,vz)
#         fitmap += p[fitparams[i]] * normedSH[fitparams[i]]
#         fiterrmap += p[fiterrparams[i]] * normedSH[fitparams[i]]
#
#     if opts['out']:
#         hp.write_map(out, fitmap)
#         return
#
#     return fitmap


# Get dipole/quadrupole fit parameters in dictionary
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
    p['nvpar'] = nvpar
    p['ierr'] = ierr

    for i, par in enumerate(fitparams):
        p[par], p['d' + par] = ROOT.Double(0), ROOT.Double(0)
        minimizer.GetParameter(i, p[par], p['d' + par])

    # For some reason there was an issue using numpy.save() when the
    # values in p were ROOT.Double(0). A quick fix is to convert
    # all the values in p to floats
    for key in p:
        p[key] = float(p[key])

    return p


# Ouput parameters from dipole/quadrupole fit
def outputFit(p, fitparams, scale):

    print '\nchi2/ndf = %.1f / %d = %.4f' % (p['chi2'], p['ndof'], p['chi2'] / float(p['ndof']))
    print '\nprob = %.2e' % (p['prob'])
    print '\nierr = {}'.format(p['ierr'])
    print '\nFit values (x %d):' % scale
    for par in fitparams:
        print ' %s = %.3f +/- %.3f' % (par, scale * p[par], scale * p['d' + par])
