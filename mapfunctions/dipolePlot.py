#!/usr/bin/env python

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import glob, os, tabulate
import mapFunctions

#from anisotropy.icesim.analysis import readDist as readDist_IC
#from anisotropy.topsim.analysis import readDist as readDist_IT
#from anisotropy.mapFunctions.energyCuts import getEbins, getEnergyMaps

from frankAnisotropy.icesim.analysis import readDist as readDist_IC
from frankAnisotropy.topsim.analysis import readDist as readDist_IT
from frankAnisotropy.mapFunctions.energyCuts import getEbins, getEnergyMaps

# Return dipole strength, error, theta, and phi
def getDipole(p):

    py  = p['Y(1,-1)']
    pz  = p['Y(1,0)']
    px  = p['Y(1,1)']
    dpy = p['dY(1,-1)']
    dpz = p['dY(1,0)']
    dpx = p['dY(1,1)']

    # Calculate dipole strength
    #p0  = np.sqrt(3.*(px**2 + py**2 + pz**2)/(4.*np.pi))
    p0  = np.sqrt(px**2 + py**2 + pz**2)
    dp0 = np.sqrt(((px*dpx)**2 + (py*dpy)**2 + (px*dpz)**2) / p0**2)

    # Calculate dipole phase
    theta, phi = hp.vec2ang(np.array([px, py, pz]))

    # phi = arctan(y/x)
    dphidy =  px / (px**2 + py**2)
    dphidx = -py / (px**2 + py**2)
    dPhi = np.sqrt((dphidx * dpx)**2 + (dphidy * dpy)**2)

    # theta = arccos(z / np.sqrt(x**2 + y**2 + z**2))
    dthetadx = px * pz / (np.sqrt((px**2 + py**2)/(px**2+py**2+pz**2)) * \
            np.sqrt((px**2 + py**2 + pz**2)**(3)))
    dthetady = py * pz / (np.sqrt((px**2 + py**2)/(px**2+py**2+pz**2)) * \
            np.sqrt((px**2 + py**2 + pz**2)**(3)))
    dthetadz = -np.sqrt((px**2 + py**2)/(px**2+py**2+pz**2)) / \
            np.sqrt(px**2 + py**2 + pz**2)
    dTheta = np.sqrt((dthetadx*dpx)**2 + (dthetady*dpy)**2 + (dthetadz*dpz)**2)

    return p0, dp0, theta, dTheta, phi, dPhi


def configPlot(scale=1e4, verbose=False):

    prefix = '/data/user/fmcnally/anisotropy/maps/merged/'
    fileList = np.array(glob.glob(prefix + 'IC?*_24H_sid.fits'))
    alpha = 1/20.
    deg = 180./np.pi

    # Sort fileList according to configs
    configs = np.array([os.path.basename(f).split('_')[0] for f in fileList])
    csort = configs.argsort()
    configs, fileList = configs[csort], fileList[csort]

    dipoleParams = np.zeros((len(fileList), 6))
    for i, f in enumerate(fileList):

        data = mapFunctions.getMap(f, mapName='data')
        bg = mapFunctions.getMap(f, mapName='bg')
        p = mapFunctions.multifit(1, data, bg, alpha, params=True)
        if verbose:
            print configs[i]
            ptable = [[key, p[key], p['d'+key]] for key in sorted(p.keys())]
            print tabulate(ptable, headers=['Ylm','Value','Error'])
        dipoleParams[i] = getDipole(p)

    p0, dp0, theta, dTheta, phi, dPhi = dipoleParams.T
    if verbose:
        for i, config in enumerate(configs):
            print config, ':'
            print 'Theta: %.3f +/- %.3f' % (theta[i]*deg, dTheta[i]*deg)
            print 'Phi: %.3f +/- %.3f' % (phi[i]*deg, dPhi[i]*deg)

    fig = plt.figure(figsize=(17,6))
    x = range(len(configs))
    pltParams = {'fmt':'.','lw':2,'ms':14}

    ax = fig.add_subplot(121)
    ax.set_title('Dipole Strength vs Energy')
    ax.set_xlabel(r'log$_{10}$(E/GeV)')
    ax.set_ylabel('Dipole Strength (x %s)' % scale)
    ax.errorbar(x, p0*scale, yerr=dp0*scale, **pltParams)
    ax.set_xlim(-1, x[-1]+1)
    ax.set_xticks(x)
    ax.set_xticklabels(configs)

    ax = fig.add_subplot(122)
    ax.set_title('Dipole Phase vs Energy')
    ax.set_xlabel(r'log$_{10}$(E/GeV)')
    ax.set_ylabel(r'Dipole Phase ($\phi^{\circ}$)')

    phi, dPhi = phi*deg, dPhi*deg
    phi[phi>180] -= 360
    ax.errorbar(x, phi, yerr=dPhi, **pltParams)

    ax.set_xlim(-1, x[-1]+1)
    ax.set_xticks(x)
    ax.set_xticklabels(configs)

    plt.show()


#def energyPlot(scale=1e4, verbose=False):
def energyPlot(scale=1e4,offset=180, out=False, batch=False, verbose=False):

    prefix = '/data/user/fmcnally/anisotropy/maps/merged'
    alpha = 1/20.
    deg = 180./np.pi
    colorDict = {0:'b'}

    # IceCube information
    files, meds, sigL, sigR = [],[],[],[]
    ebins = getEbins()
    files = getEnergyMaps('IC')
    for i in range(len(ebins)-1):
        distInfo = readDist_IC('IC', ebins[i], ebins[i+1])
        meds += [distInfo[0]]
        sigL += [distInfo[1]]
        sigR += [distInfo[2]]

    # IceTop information
    idx0 = len(files)
    colorDict[idx0] = 'r'
    files += [['%s/IT_24H_sid_STA8.fits' % (prefix)]]
    distInfo = readDist_IT('IT', 8, 100)
    meds += [distInfo[0]]
    sigL += [distInfo[1]]
    sigR += [distInfo[2]]

    dipoleParams = np.zeros((len(files), 6))
    for i, files_i in enumerate(files):
        data = mapFunctions.getMap(*files_i, mapName='data')
        bg = mapFunctions.getMap(*files_i, mapName='bg')
        p = mapFunctions.multifit(1, data, bg, alpha, params=True)
        if verbose:
            print meds[i]
            print('key = {}'.format(sorted(p.keys())))
            #ptable = [[key, p[key], p['d'+key]] for key in sorted(p.keys())]
            ptable = [[key, p[key]] for key in sorted(p.keys())]
            #print tabulate(ptable, headers=['Ylm','Value','Error'])
            #print tabulate(ptable, headers=['Ylm','Value'])
            print('files_i = {}'.format(files_i))
            for key in sorted(p.keys()):
			    print('p[{}] = {}'.format(key,p[key]))
            print('p[Y(1,1)]+p[Y(1,-1)] = {}'.format(np.sqrt(3.0/(4.0*np.pi))*np.power(p['Y(1,1)'],2.0)+np.power(p['Y(1,-1)'],2.0)))			
        dipoleParams[i] = getDipole(p)

        p0, dp0, thetas, dThetas, phis, dPhis = dipoleParams.T

        if verbose:
            for i, median in enumerate(meds):
                print 'Median energy : %.3f' % median
                print 'Theta: %.3f +/- %.3f' % (thetas[i]*deg, dThetas[i]*deg)
                print 'Phi: %.3f +/- %.3f' % (phis[i]*deg, dPhis[i]*deg)

    #fig = plt.figure(figsize=(17,6))
    fig = plt.figure()
    #pltParams = {'fmt':'.','lw':2,'ms':14}
    pltParams = {'lw':2,'ms':14}

    #ax = fig.add_subplot(121)
    #ax.set_title('Dipole Strength vs Energy')
    #ax.set_xlabel(r'log$_{10}$(E/GeV)')
    #ax.set_ylabel('Dipole Strength (x %s)' % scale)
    #ax.errorbar(x, p0*scale, yerr=dp0*scale, **pltParams)

    #ax = fig.add_subplot(122)
    ax = fig.add_subplot(111)
    ax.set_title('Dipole Phase vs Energy')
    ax.set_xlabel(r'Reconstructed Energy (log$_{10}$(E/GeV))')
    ax.set_ylabel(r'Dipole Phase ($\phi^{\circ}$)')

    for i in range(0, len(meds), idx0):
        x = meds[i:i+idx0]
        xerr = [sigL[i:i+idx0], sigR[i:i+idx0]]
        phi, dPhi = phis[i:i+idx0]*deg, dPhis[i:i+idx0]*deg
        # Phi offset (requires renaming ylabels)
        pcut = phi > offset
        phi[pcut] -= offset
        phi[np.logical_not(pcut)] += (360-offset)
        ax.errorbar(x, phi, xerr=xerr, yerr=dPhi,
                fmt=colorDict[i]+'.', **pltParams)

    ax.set_ylim(-10,370)
    ax.set_yticks(range(0, 361, 45))
    ylabels  = [str(i) for i in range(offset, 360-1, 45)]
    ylabels += [str(i) for i in range(0, offset+1, 45)]
    ax.set_yticklabels(ylabels)

    plt.axhline(360-offset, linewidth=1.5, linestyle='--', color='black')

    if out:
        plt.savefig(out, dpi=300, bbox_inches='tight')
    if not batch:
        plt.show()


def minimumPlot(offset=180, out=False, batch=False, verbose=False):

    prefix = '/data/user/fmcnally/anisotropy/maps/merged'
    alpha = 1/20.
    deg = 180./np.pi
    colorDict = {'IC':'b','IT':'r'}

    # Create dictionaries to store IceCube and IceTop data separately
    files, meds, sigL, sigR = {},{},{},{}
    thetas, phis, dThetas, dPhis = {},{},{},{}

    # IceCube information
    meds['IC'], sigL['IC'], sigR['IC'] = [],[],[]
    ebins = getEbins()
    files['IC'] = getEnergyMaps('IC')
    for i in range(len(ebins)-1):
        distInfo = readDist_IC('IC', ebins[i], ebins[i+1])
        meds['IC'] += [distInfo[0]]
        sigL['IC'] += [distInfo[1]]
        sigR['IC'] += [distInfo[2]]

    # IceTop information
    files['IT'], meds['IT'], sigL['IT'], sigR['IT'] = [],[],[],[]
    files['IT'] += [['%s/IT_24H_sid_STA8.fits' % (prefix)]]
    distInfo = readDist_IT('IT', 8, 100)
    meds['IT'] += [distInfo[0]]
    sigL['IT'] += [distInfo[1]]
    sigR['IT'] += [distInfo[2]]

    for detector in files.keys():
        fileList = files[detector]
        theta, phi = np.zeros((2, len(fileList)))
        mapParams = {'mask':True,'smooth':20}
        for i, f in enumerate(fileList):
            sig   = mapFunctions.getMap(*f, mapName='signal', **mapParams)
            ri    = mapFunctions.getMap(*f, mapName='relint', **mapParams)
            rierr = mapFunctions.getMap(*f, mapName='relerr', **mapParams)
            c0 = (sig==hp.UNSEEN)
            sig[c0] = 0
            nside = hp.npix2nside(len(sig))
            minPix = sig.argmin()
            theta[i], phi[i] = hp.pix2ang(nside, minPix)

        if verbose:
            for i, median in enumerate(medians):
                print 'Median energy : %.3f' % median
                print 'Theta: %.3f' % (theta[i]*deg)
                print 'Phi: %.3f' % (phi[i]*deg)
                #print 'Theta: %.3f +/- %.3f' % (theta[i]*deg, dTheta[i]*deg)
                #print 'Phi: %.3f +/- %.3f' % (phi[i]*deg, dPhi[i]*deg)

        thetas[detector]  = theta
        #dThetas[detector] = dTheta
        phis[detector]  = phi
        #dPhis[detector] = dPhi

    #fig = plt.figure(figsize=(17,6))
    fig = plt.figure()
    #pltParams = {'fmt':'.','lw':2,'ms':14}
    pltParams = {'lw':2,'ms':14}

    #ax = fig.add_subplot(121)
    #ax.set_title('Dipole Strength vs Energy')
    #ax.set_xlabel(r'log$_{10}$(E/GeV)')
    #ax.set_ylabel('Dipole Strength (x %s)' % scale)
    #ax.errorbar(x, p0*scale, yerr=dp0*scale, **pltParams)

    #ax = fig.add_subplot(122)
    ax = fig.add_subplot(111)
    #ax.set_title('Dipole Phase vs Energy')
    ax.set_xlabel(r'Reconstructed Energy (log$_{10}$(E/GeV))')
    ax.set_ylabel(r'Dipole Phase ($\phi^{\circ}$)')

    for detector in files.keys():
        x = meds[detector]
        xerr = [sigL[detector], sigR[detector]]
        phi = phis[detector]*deg 
        #phi, dPhi = phis[detector]*deg, dPhis[detector]*deg
        # Phi offset (requires renaming ylabels)
        pcut = phi > offset
        phi[pcut] -= offset
        phi[np.logical_not(pcut)] += (360-offset)
        ax.errorbar(x, phi, xerr=xerr,
                fmt=colorDict[detector]+'.', **pltParams)
        #ax.errorbar(x, phi, xerr=xerr, yerr=dPhi,
        #        fmt=colorDict[detector]+'.', **pltParams)

    ax.set_ylim(-10,370)
    ax.set_yticks(range(0, 361, 45))
    ylabels  = [str(i) for i in range(offset, 360-1, 45)]
    ylabels += [str(i) for i in range(0, offset+1, 45)]
    ax.set_yticklabels(ylabels)

    plt.axhline(360-offset, linewidth=1.5, linestyle='--', color='black')

    if out:
        plt.savefig(out, dpi=300, bbox_inches='tight')
    if not batch:
        plt.show()


if __name__ == "__main__":

    #configPlot()
    energyPlot(verbose=True)

