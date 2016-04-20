#!/usr/bin/env python

#=========================================================================
# File Name     : SHdipole.py
# Description   :
# Creation Date : 04-08-2016
# Last Modified : Fri 08 Apr 2016 12:10:29 PM CDT
# Created By    : James Bourbeau
#=========================================================================

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import argparse

from mapfunctions.mapFunctions import getMap
from mapfunctions.projFunctions import *

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Calculate phase from 2D fit map')
    parser.add_argument('--init', dest='init', type=float, \
        help='Number of RA bins to be used in projection')
    parser.add_argument('--step', dest='step', type=float, \
        help='Number of RA bins to be used in projection')
    parser.add_argument('-f','--file', dest='file',\
        default='/data/user/jbourbeau/anisotropy/maps/merged/IC_24H_sid_4-4.25GeV.fits',
        help='.fits file to be analyzed')
    parser.add_argument('--projRI', dest='projRI', default=False, \
        action='store_true', help='Plot the projected RI from 2D fitmap')
    parser.add_argument('--dipole', dest='dipole', default=False, \
        action='store_true', help='Plot the dipole phase and amplitude from 2D fitmap')

    args = parser.parse_args()
    kwargs = vars(args).copy()
    defaults = {'smooth':5, 'stype':'tophat', 'swindow':3,\
                'verbose':False, 'decmin':-90., 'decmax':-25.}
    opts = {k:kwargs[k] for k in kwargs if k not in defaults}
    opts.update({k:defaults[k] for k in defaults})
    # print('opts = {}'.format(opts))

    init = 1e-15
    step = 1e-06

    chi2list = ['standard', 'tibetonly', 'RI']
    labels = {'standard':'IceCube only', 'RI':'Combined', 'tibetonly':'Tibet only'}
    colors = {'standard':'blue', 'RI':'green', 'tibetonly':'red'}

    if opts['dipole']:
        fig = plt.figure()
        ax = plt.gca()
        for chi2 in ['standard', 'RI']:
            amp = []
            amperr = []
            phase = []
            phaseerr = []
            for l in range(1,6):
                # Load appropreiate SH coeff dictionary
                p = np.load('SHCoeff/SHCoeff{}_{}/{}chi2_lmax_{}.npy'.format(init,step,chi2,l))
                p = p.item()

                # Useful breakdown of l, m values to be used
                nsph = sum([2*l_i+1 for l_i in range(l+1)])
                lvals = [[l_i]*(2*l_i+1) for l_i in range(l+1)]
                mvals = [[m for m in range(-l_i, l_i+1)] for l_i in range(l+1)]
                lvals = [item for sublist in lvals for item in sublist]
                mvals = [item for sublist in mvals for item in sublist]
                fitparams = ['Y(%i,%i)' % (lvals[i], mvals[i]) for i in range(nsph)]
                fiterrparams = ['dY(%i,%i)' % (lvals[i], mvals[i]) for i in range(nsph)]
                normedSH = np.load('/home/jbourbeau/anisotropy/dev/normedSH.npy')
                normedSH = normedSH.item()
                nside = 64
                npix = hp.nside2npix(nside)
                fitmap = np.zeros(npix)
                fiterrmap = np.zeros(npix)
                for i in range(nsph):
                    fitmap += p[fitparams[i]] * normedSH[fitparams[i]]
                    fiterrmap += p[fiterrparams[i]] * normedSH[fitparams[i]]

                # Dipole amplitude from 2D fit
                # amp.append((max(fitmap)-min(fitmap))/(2+max(fitmap)+min(fitmap)))
                amp.append(np.sqrt((3/(4*np.pi))*(p['Y(1,1)']**2+p['Y(1,0)']**2+p['Y(1,-1)']**2)))
                amperr.append(np.sqrt((3/(4*np.pi)))*2*np.sqrt((p['Y(1,1)']*p['dY(1,1)'])**2+(p['Y(1,0)']*p['dY(1,0)'])**2+(p['Y(1,-1)']*p['dY(1,-1)'])**2))
                # amp.append((max(relint)-min(relint))/(max(relint)+min(relint)))

                # Dipole phase from 2D fit
                phase.append(np.arctan(p['Y(1,1)']/p['Y(1,-1)'])*180./np.pi)

                prefactor = 1./(1+(p['Y(1,1)']/p['Y(1,-1)'])**2)
                first = (p['dY(1,1)']/p['Y(1,-1)'])**2
                second = ((p['dY(1,-1)']*p['Y(1,1)'])/(p['Y(1,-1)']**2))**2
                error = prefactor*np.sqrt(first+second)*180./np.pi
                phaseerr.append(error)

            x = np.arange(1,6)
            # Only include first three points for tibetonly
            if chi2=='tibetonly':
                x = x-0.05
                x = x[:3]
                phase = phase[:3]
                phaseerr = phaseerr[:3]
                amp = amp[:3]
            # Slight shift in x values so error bars can be seen
            if chi2=='RI':
                x = x+0.05
            # plt.errorbar(x,phase,yerr=phaseerr,marker='.',markersize=10,linestyle='-.',label=labels[chi2]+r' $\chi^2$')
            plt.errorbar(x,amp,yerr=amperr,marker='.',markersize=10,linestyle='-.',label=labels[chi2]+r' $\chi^2$')
            # plt.plot(x,amp,label=labels[chi2]+r' $\chi^2$',marker='.',linestyle='-.',markersize=10)
            # plt.title(r'Phase from 2D fit maps')
        tPars = {'fontsize':16}
        ax.set_xlabel(r'$l_{}$'.format('{max}'),**tPars)
        ax.set_ylabel(r'Dipole amplitude',**tPars)
        # ax.set_ylabel(r'Dipole phase $[^{}]$'.format('{\circ}'),**tPars)
        # plt.legend()
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
            ncol=3, mode="expand", borderaxespad=0.)
        # ax.set_ylim(-50,50)
        ax.set_ylim(0,0.01)
        ax.set_xlim(0.75,5.25)
        # ax.set_xlim(0.75,7.25)
        plt.savefig('/home/jbourbeau/public_html/figures/collaboration-meeting-2016/dipole-amp_SHCoeff{}_{}_lmax5.png'.format(init,step),dpi=300, bbox_inches='tight')
        # plt.savefig('/home/jbourbeau/public_html/figures/collaboration-meeting-2016/dipole-phase_SHCoeff{}_{}_lmax5.png'.format(init,step),dpi=300, bbox_inches='tight')

    #--------------------------------------------------------------

    if opts['projRI']:
        for l in range(1,8):
            fig = plt.figure()
            ax = plt.gca()
            # Get Tibet proj relint data points
            tibetdata = np.loadtxt('Tibet_data.dat')
            tibetRA = tibetdata[:,0]
            tibetRAerr = [10.]*18
            tibetRI = tibetdata[:,1]-1.
            tibetRIerr = np.array([2.8e-5]*18)
            tibetpopt, tibetperr, tibetchi2 = getHarmonicFitParams(tibetRA, tibetRI, 1, tibetRIerr)
            tibetprojRI = cosFit(tibetRA,*tibetpopt[:3])
            plt.errorbar(tibetRA,tibetRI,xerr=tibetRAerr,yerr=tibetRIerr,marker='.',
                linestyle='None',color='red')
            plt.plot(tibetRA,cosFit(tibetRA,*tibetpopt[:3]),color='red',label='Tibet proj. RI')
            for chi2 in ['standard', 'RI']:
                # Load appropreiate SH coeff dictionary
                p = np.load('SHCoeff/SHCoeff{}_{}/{}chi2_lmax_{}.npy'.format(init,step,chi2,l))
                p = p.item()

                # Useful breakdown of l, m values to be used
                nsph = sum([2*l_i+1 for l_i in range(l+1)])
                lvals = [[l_i]*(2*l_i+1) for l_i in range(l+1)]
                mvals = [[m for m in range(-l_i, l_i+1)] for l_i in range(l+1)]
                lvals = [item for sublist in lvals for item in sublist]
                mvals = [item for sublist in mvals for item in sublist]
                fitparams = ['Y(%i,%i)' % (lvals[i], mvals[i]) for i in range(nsph)]
                fiterrparams = ['dY(%i,%i)' % (lvals[i], mvals[i]) for i in range(nsph)]
                normedSH = np.load('/home/jbourbeau/anisotropy/dev/normedSH.npy')
                normedSH = normedSH.item()
                nside = 64
                npix = hp.nside2npix(nside)
                fitmap = np.zeros(npix)
                fiterrmap = np.zeros(npix)
                for i in range(nsph):
                    fitmap += p[fitparams[i]] * normedSH[fitparams[i]]
                    fiterrmap += p[fiterrparams[i]] * normedSH[fitparams[i]]

                opts.update({'lmax':1,'nbins':18,'ramin':0.,'ramax':360.,'decmin':-30.,'decmax':90.})
                ra, ri, ra_err, ri_err = getRIRAProj(fitmap, fiterrmap, **opts)
                # ra, ri, ra_err = getRIRAProj(fitmap, **opts)
                # Get fitmap dipole amp and phase that Tibet would see
                amp, amp_err, phase, phase_err = getProjDipole(fitmap,fiterrmap, **opts)
                print('\nlmax = {}'.format(l))
                print('chi2 = {}'.format(chi2))
                print('amp = {:e}'.format(amp))
                print('amp_err = {:e}'.format(amp_err))
                print('phase = {}'.format(phase))
                print('phase_err = {}\n'.format(phase_err))
                # amp, amp_err, phase, phase_err = getProjDipole(fitmap, **opts)
                popt, perr, chisquared = getHarmonicFitParams(ra, ri, 1)
                projRI = cosFit(ra,*popt[:3])

                plt.errorbar(ra,ri,xerr=ra_err,yerr=ri_err,marker='.',linestyle='None',color=colors[chi2])
                plt.plot(ra,cosFit(ra,*popt[:3]),color=colors[chi2],label='{} $\chi^2$'.format(labels[chi2]))
                # plt.plot(ra,cosFit(ra,*popt[:3]),color=colors[chi2],label='Dipole Fit')
                # ax.annotate(r'{} $\chi^2$'.format(labels[chi2]), xy=(125, -0.00085))
                # ax.annotate(r'Dipole = {:1.1e}$\pm${:1.1e}'.format(amp,amp_err), xy=(125, -0.0010))
                # ax.annotate(r'Phase = {:2.1f}$\pm${:1.1f}'.format(phase,phase_err), xy=(125, -0.00115))

            tPars = {'fontsize':16}
            ax.set_xlim(0.0,360.)
            ax.set_ylim(-0.0015,0.0015)
            #ax.set_ylim(-0.0025,0.0025)
            ax.set_ylabel(r'Relative Intensity', **tPars)
            ax.set_xlabel(r'Right Ascension $[^{}]$'.format('{\circ}'), **tPars)
            ax.invert_xaxis()
            # plt.legend()
            plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                   ncol=3, mode="expand", borderaxespad=0.)
            plt.savefig('/home/jbourbeau/public_html/figures/collaboration-meeting-2016/projRI/projRI_lmax{}.png'.format(l),dpi=300, bbox_inches='tight')
