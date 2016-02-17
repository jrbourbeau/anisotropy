#!/usr/bin/env python

###############################################################################
## An attempt to update the TimeScramble.cc code to run with ShowerLLH       ##
###############################################################################

from icecube import ShowerLLH
from icecube import astro
from numpy import *
import healpy as hp


# Produce time scrambled maps for data given start and end date
def timeScramble(data, nInt, method='equatorial', out=False):

    # Constants
    nSide = 64
    nBGResample = 20            # Number of times to resample events for bg
    alpha = 1./nBGResample
    if 24%nInt != 0:
        print 'Integration time does not work cleanly into 24 hours!'
        sys.exit(1)
    dt = nInt / 24.             # Integration time (days)

    # Take care of weighting by creating a multiplier for each event
    keyList = ['ShowerPlane_azimuth', 'ShowerPlane_zenith', 'mjd']
    nEvents = len(data['mjd'])
    data['mult'] = zeros(nEvents, dtype=int)
    data['mult'] += 3 * data['w3']
    data['mult'] += 1 * data['w1']
    data['mult'] += random.randint(1,3,nEvents) * data['w1.5']

    # Setup maps
    nPixels = hp.nside2npix(nSide)
    maps = {}
    #for key in ['Local', 'Data', 'BG', 'Signal']:
    for key in ['Data', 'BG']:
        maps[key] = zeros(nPixels)

    temp, temp_wt = {},{}
    mjd_start = data['mjd'][0]
    finished = False

    # Loop through integration time periods
    while not finished:

        # Limit time range
        temp_cut = (data['mjd']>mjd_start) * (data['mjd']<mjd_start+dt)
        if temp_cut.sum() == 0:
            finished = True
            continue
        for key in data.keys():
            temp[key] = data[key][temp_cut]

        # Create new, weighted arrays
        for key in keyList:
            temp_wt[key] = []
        for i in range(len(temp['mjd'])):
            for j in range(temp['mult'][i]):
                for key in keyList:
                    temp_wt[key].append(temp[key][i])
        for key in keyList:
            temp_wt[key] = asarray(temp_wt[key])

        # Update maps for this time interval
        temp_maps = cscrambler(temp_wt, nSide, nBGResample, method)
        for key in temp_maps.keys():
            maps[key] += temp_maps[key]

        mjd_start += dt

    maps['BG'] *= alpha
    #maps['Signal'] = LMSignificance(maps['Data'], maps['BG'], alpha)

    # Save to file
    if out:
        for key in maps.keys():
            hp.write_map(out+'_'+key.lower()+'.fits', maps[key])

    return maps


# Calculate maps and perform time scrambling
def cscrambler(data, nSide, nBGResample, method):

    # Setup maps
    test = ShowerLLH.TimeScramble()
    temp_maps = {}
    nPixels = hp.nside2npix(nSide)
    #for key in ['Local', 'Data', 'BG']:
    for key in ['Data', 'BG']:
        temp_maps[key] = zeros(nPixels)

    # Begin looping through data
    mjd = data['mjd']
    theta = data['ShowerPlane_zenith']
    phi = data['ShowerPlane_azimuth']

    # Fill local and data maps
    local_pix = hp.ang2pix(nSide, theta, phi)
    dec, ra = array(map(l2e, mjd, theta, phi)).T
    data_pix = hp.ang2pix(nSide, dec, ra)
    for i in range(len(local_pix)):
        #temp_maps['Local'][local_pix[i]] += 1.0
        temp_maps['Data'][data_pix[i]] += 1.0

    # Fill background map
    bg = test.cscramble(mjd, theta, phi, nSide, nBGResample, method)
    temp_maps['BG'] = asarray(bg)

    return temp_maps


# Li Ma Significance
def LMSignificance(nData, nBG, alpha):

    inv_err = geterr()['invalid']
    div_err = geterr()['divide']
    seterr(divide  = 'ignore')
    seterr(invalid = 'ignore')

    Non  = nData
    Noff = nBG / alpha
    sn = sign(Non - nBG)
    sigma = sn * sqrt(2 * (Non * log(((1+alpha)*Non) / (alpha*(Non+Noff))) + 
            Noff * log(((1+alpha)*Noff) / (Non+Noff))))

    # Check for nan's
    nanCatch = (sigma!=sigma)
    sigma[nanCatch] = 0

    seterr(divide  = div_err)
    seterr(invalid = inv_err)

    return sigma


# General local to equatorial conversion
def l2e(mjd, theta, phi, method='equatorial'):

    ice = astro.IceCubeDetector()
    local = astro.LocalCoord()
    t = astro.Time()

    t.SetTime(mjd)
    local.SetLocalRad(theta, phi)

    if (method == 'equatorial'):
        eqApparent = ice.LocalToEquatorial(local, t)
        phi = eqApparent.GetRaRad()
    if (method == 'anti_sid'):
        eqApparent = ice.LocalToEquatorial_FromAntiSid(local, t)
        phi = eqApparent.GetRaRad()
    if (method == 'ext_sid'):
        eqApparent = ice.LocalToEquatorial_FromExtSid(local, t)
        phi = eqApparent.GetRaRad()
    if (method == 'solar'):
        eqApparent = ice.LocalToEquatorial(local, t)
        eqSolar = ice.PlanetToEquatorial(0, t)
        phi = eqApparent.GetRaRad() - eqSolar.GetRaRad()

    theta = pi/2 - eqApparent.GetDecRad()
    if phi < 0:
        phi += 2*pi

    return theta, phi


""" Pure python implementation of scrambling code - ~10x slower """
def scrambler(data, nSide, nBGResample, method):

    # Setup maps
    temp_maps = {}
    nPixels = hp.nside2npix(nSide)
    for key in ['Local', 'Data', 'BG']:
        temp_maps[key] = zeros(nPixels)

    # Begin looping through data
    mjd = data['mjd']
    theta = data['ShowerPlane_zenith']
    phi = data['ShowerPlane_azimuth']

    # Fill local and data maps
    local_pix = hp.ang2pix(nSide, theta, phi)
    dec, ra = array(map(l2e, mjd, theta, phi)).T
    data_pix = hp.ang2pix(nSide, dec, ra)
    for i in range(len(local_pix)):
        temp_maps['Local'][local_pix[i]] += 1.0
        temp_maps['Data'][data_pix[i]] += 1.0

    # Function for getting nBGResample background counts with random times
    def getBG(th, ph):
        rndMJD = random.choice(mjd, nBGResample)
        dec, ra = array([l2e(mjd_one, th, ph) for mjd_one in rndMJD]).T
        pixelID = hp.ang2pix(nSide, dec, ra)
        return pixelID

    # Fill background map
    bg_pix = array(map(getBG, theta, phi)).flatten()
    for pix in bg_pix:
        temp_maps['BG'][pix] += 1.0

    return temp_maps



