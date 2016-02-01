#!/usr/bin/env python

from numpy import *
import pickle
import matplotlib.pyplot as plt
from data import bakhPlot

def plotter(month):

    fl = open('collab/'+month+'_hists.pkl', 'rb')
    q = pickle.load(fl)
    Emids = q['Emids']
    fl.close()

    out = month+'_myflux'

    # Starting table setup
    f = plt.figure()
    ax1 = f.add_subplot(1,1,1)
    ax1.set_title(month+' Energy Spectrum after Unfolding')
    #ax1.set_title(month+' Energy Spectrum after Unfolding (E^2.7)')
    ax1.set_xlabel('Log10(Energy/GeV)')
    ax1.set_ylabel('Flux (counts/(m^2 s ster GeV))')
    #ax1.set_ylabel('Flux (counts/(m^2 s ster GeV)) * E^2.7')

    flux, err = {},{}
    iter = '5'
    spec = 0
    scale = (10**Emids)**spec

    # Simplify plotting
    def myebar(key, m='.', c='Default'):
        colorDict = {'F':'r', 'P':'b', 'A':'k'}
        labelDict = {'F':'Iron', 'P':'Proton', 'A':'Total'}
        myFlux, myErr = flux[key], err[key]
        e = key[0]
        if c=='Default':
            c = colorDict[e]
        myMark  = c + m
        myLabel = labelDict[e]
        if len(key) > 1:
            myLabel += ' '+key[1:]
        ax1.errorbar(Emids, myFlux, yerr=myErr, fmt=myMark, label=myLabel)

    # Unfolded flux plot
    for key in ['P', 'F', 'A']:
        flux[key] = q['flux']['none_'+key+'_'+iter] * scale
        err[key] = q['err']['none_'+key+'_'+iter] * scale

    # Iteration plot
    #for key in ['P', 'F', 'A']:
    for key in ['A']:
        for i in range(1,len(q['chi2'])+1):
            i = str(i)
            flux[key+'iter'+i] = q['flux']['none_'+key+'_'+i] * scale
            err[key+'iter'+i] = q['err']['none_'+key+'_'+i] * scale

    # Hemisphere plot
    for key in ['P', 'F', 'A']:
        for dir in ['left', 'right']:
            flux[key+dir] = q['flux'][dir+'_'+key+'_'+iter] * scale
            err[key+dir] = q['err'][dir+'_'+key+'_'+iter] * scale

    # Odd/Even plot
    for key in ['P', 'F', 'A']:
        for dir in ['odds', 'evens']:
            flux[key+dir] = q['flux'][dir+'_'+key+'_'+iter] * scale
            err[key+dir] = q['err'][dir+'_'+key+'_'+iter] * scale

    # Chi2 plot
    #chi2List = []
    #old = flux['A1']
    #for i in range(2, len(q['chi2'])+1):
    #    new = flux['A'+str(i)]

    # Bakhtiyar's plot
    B_mids, B_flux, B_relup, B_reldn = bakhPlot()
    B_flux *= ((10**B_mids)**spec)
    B_errup = B_flux * B_relup
    B_errdn = B_flux * B_reldn
    #ax1.errorbar(B_mids, B_flux, yerr=(B_errup, B_errdn), fmt='gx', label='Bakhtiyar')

    # Sample plotting lists
    keyList = [k for k in flux.keys() if len(k)==1]
    iterList = [k for k in flux.keys() if 'iter' in k]
    leftList  = [k for k in flux.keys() if 'left' in k]
    rightList = [k for k in flux.keys() if 'right' in k]
    oddList  = [k for k in flux.keys() if 'odds' in k]
    evenList = [k for k in flux.keys() if 'evens' in k]

    # Plot
    myebar('A')
    #for key in leftList:
    #    myebar(key, m='-')
    #for key in rightList:
    #    myebar(key)

    # Refine plotting parameters
    ax1.set_yscale('log')
    ax1.legend(loc='lower left')
    ax1.set_xlim((6.2, 9.5))
    ax1.set_ylim((10**(-23), 10**(-11)))
    #ax1.set_ylim((10**(2), 10**(5)))
    plt.savefig('collab/pics/'+out+'.png')
    plt.show()


if __name__ == "__main__":

    #for month in ['All']:
    for month in ['201006', '201007', '201008', '201012', '201101', 'All']:

        plotter(month)
