#!/usr/bin/env python

from numpy import *
import matplotlib.pyplot as plt
from useful import getEbins


def getLLHhists(s, cut, binWidth):

    # Setup params
    r = log10(s['ML_energy'])
    Ebins = getEbins()
    pf_ratio = s['pLLH'] - s['fLLH']
    llh_bins = linspace(pf_ratio.min(), pf_ratio.max(), 200)
    llh_mids = (llh_bins[1:]+llh_bins[:-1]) / 2.
    e_start  = Ebins[0]
    newEbins = array([e_start])

    # Begin writing hists
    h = {}
    while e_start < Ebins[-1]:
        e_end = e_start + binWidth
        newEbins = append(newEbins, e_end)
        ecut = (r>=e_start) * (r<e_end)
        e_range = '%s_%s' % (e_start, e_end)
        e_start += binWidth

        # Build histograms for individual energy bins
        h[e_range] = {}
        for comp in ['P','Fe']:
            c0 = s[comp] * cut * ecut
            h[e_range][comp] = histogram(pf_ratio[c0], bins=llh_bins)[0]

    return h


def makeLLHcut(s, cut, binWidth, p, plot=False):

    # Setup params
    r = log10(s['ML_energy'])
    Ebins = getEbins()
    #Ebins = Ebins[40:]
    pf_ratio = s['pLLH'] - s['fLLH']
    llh_bins = linspace(pf_ratio.min(), pf_ratio.max(), 200)
    llh_mids = (llh_bins[1:]+llh_bins[:-1]) / 2.
    e_start = Ebins[0]
    newbins = array([e_start])
    pllh, fllh = array([]), array([])

    while e_start < Ebins[-1]:
        e_end = e_start + binWidth
        newbins = append(newbins, e_end)
        ecut = (r>=e_start) * (r<e_end)
        e_start += binWidth

        # Build histograms
        h = {}
        for true_comp in ['P','Fe']:
            c0 = s[true_comp] * cut * ecut
            h[true_comp] = histogram(pf_ratio[c0], bins=llh_bins)[0]

        #Calculate what fraction of total events are each composition
        p_fraction, f_fraction = zeros((2, len(llh_mids)))
        Ntot_p, Ntot_f = zeros((2, len(llh_mids)))
        for j in range(len(llh_mids)):
            Ntot_p[j] = float(h['P'][j:].sum() + h['Fe'][j:].sum())
            Ntot_f[j] = float(h['P'][:j].sum() + h['Fe'][:j].sum())
            p_fraction[j] = h['P'][j:].sum() / Ntot_p[j]
            f_fraction[j] = h['Fe'][:j].sum() / Ntot_f[j]
            if p_fraction[j]!=p_fraction[j]:
                p_fraction[j] = 0
            if f_fraction[j]!=f_fraction[j]:
                f_fraction[j] = 0

        # Limit to bins where over p% of the events are the desired comp
        p_cut = p_fraction > p
        f_cut = f_fraction > p
        # Select the bin among those with the highest counts
        pllh = append(pllh, llh_mids[(Ntot_p * p_cut).argmax()])
        fllh = append(fllh, llh_mids[(Ntot_f * f_cut).argmax()])

    if plot:
        f = plt.figure()
        ax1 = f.add_subplot(1,1,1)
        ax1.set_title('Likelihood ')
        ax1.set_xlabel('Energy (log10(GeV))')
        ax1.set_ylabel('pLLH - fLLH')
        x = (newbins[1:] + newbins[:-1]) / 2.
        plt.plot(x, pllh, 'b.', label='P')
        plt.plot(x, fllh, 'r.', label='Fe')
        plt.ylim(-10,10)
        plt.legend(loc='upper right')
        plt.show()

    # Take average if they pass each other?
    for i in range(len(pllh)):
        if pllh[i] < fllh[i]:
            pllh[i] = (pllh[i]+fllh[i])/2.
            fllh[i] = pllh[i]

    llhDict = {}
    llhDict['Ebins'] = newbins
    llhDict['pLLH']  = pllh
    llhDict['fLLH']  = fllh

    return llhDict


def getLLHcut(llhMax=10):

    # Load LLH cut information
    resourcedir = '/net/user/fmcnally/ShowerLLH/resources/'
    cutFile = resourcedir + 'llh_cuts.npy'
    llhDict = load(cutFile)
    llhDict = llhDict.item()

    # Restrict to energy range where likelihoods are stable
    ## Note : this will fail if there are bad likelihoods at highest energy
    pCut = (abs(llhDict['pLLH']) > llhMax)
    fCut = (abs(llhDict['fLLH']) > llhMax)
    lastGood = 0
    if pCut.sum()!=0 or fCut.sum()!=0:
        p_lastGood = where(pCut)[0][-1]
        f_lastGood = where(fCut)[0][-1]
        lastGood = max([p_lastGood, f_lastGood])

    # Apply cut
    for key in llhDict.keys():
        llhDict[key] = llhDict[key][lastGood:]

    return llhDict


## Notes : a cut at 0.75 seems to be smoother than at 0.66
## 17 points seems about right (roughly 4 energy bins per point)

if __name__ == "__main__":

    # Setup
    config = 'IT73'
    resourcedir = '/net/user/fmcnally/ShowerLLH/resources/'
    outFile = resourcedir + 'llh_cuts.npy'

    # Load simulation
    print 'Loading simulation...'
    s = load_sim(config)
    cut = s['cuts']['llh']
    binWidth = 0.2
    p = 0.75

    llhDict = makeLLHcut(s, cut, binWidth, p, plot=False)
    save(outFile, llhDict)


