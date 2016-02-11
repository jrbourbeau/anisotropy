#!/usr/bin/env python

import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
import os, argparse
from tabulate import tabulate

import simFunctions_IC as simFunctions
from dataFunctions import getConfigs
from anisotropy.icesim.analysis import load_sim as load_sim_IC
from anisotropy.icesim.analysis import readDist as readDist_IC
from anisotropy.icesim.analysis import getSimFiles
from anisotropy.topsim.analysis import load_sim as load_sim_IT
from anisotropy.topsim.analysis import readDist as readDist_IT
from anisotropy.topsim.analysis import getWeights as getWeights_IT
from anisotropy.mapFunctions.energyCuts import getEbins
from useful import getMids

def makeTable(outFile):

    d = {}
    d['IC'] = {}
    d['IT'] = {}

    configList = getConfigs('IC')
    ebins = getEbins()

    for cfg in configList:
 
        d['IC'][cfg] = {}

        # Load simulation -- part by part if separate files
        simFiles = getSimFiles(cfg)
        for i in range(len(simFiles)):
            ic = load_sim_IC(cfg, part=i)
            typeList = sorted(list(set(ic['type'])))
            weights = getWeights(ic, cfg)

            # Bin in energy
            for emin, emax in zip(ebins[:-1], ebins[1:]):
                bin = '%s-%s' % (emin, emax)
                if bin not in d['IC'][cfg].keys():
                    d['IC'][cfg][bin] = {}
                emin, emax = float(emin), float(emax)
                cut = (ic['energy_sp'] >= emin) * (ic['energy_sp'] < emax)

                # Increment counts for each composition
                for ptype in typeList:
                    if ptype not in d['IC'][cfg][bin].keys():
                        d['IC'][cfg][bin][ptype] = 0.
                    c0 = ic['type'][cut] == ptype
                    d['IC'][cfg][bin][ptype] += sum(weights[c0])

    it = load_sim_IT('IT73')
    c1 = (it['reco_zenith'] == it['reco_zenith'])
    weights = getWeights_IT(it)[c1]
    tot = float(weights.sum())
    for key in sorted(list(set(it['MC_type']))):
        c2 = (it['MC_type'][c1] == key)
        n = (weights[c2]).sum()
        d['IT'][key] = n/tot

    np.save(outFile, d)


def getWeights(s, config):

    if config in ['IC86','IC86-II','IC86-III','IC86-IV']:
        return np.ones(len(s['energy']))

    from icecube.weighting.fluxes import Hoerandel
    from icecube.weighting.weighting import from_simprod

    sim = int(simFunctions.cfg2sim(config))
    nfiles, generator = from_simprod(sim)
    nfiles = int(s['nFiles'].sum())
    generator *= nfiles
    flux = Hoerandel()

    weights = flux(s['energy'], s['type']) / \
              generator(s['energy'], s['type'])

    weights[weights!=weights] = 0
    weights[weights==np.inf]  = 0

    return weights


def fullTable(d):

    ic = d['IC']
    it = d['IT']
    mytable = []

    cutoffs = [1000020040, 1000060120, 1000200400]
    groups = ['P', 'He', 'CNO', 'Fe']
    for cfg in sorted(ic.keys()):
        for ebin in sorted(ic[cfg].keys()):
            row = [cfg, ebin]
            compList = np.array(ic[cfg][ebin].keys())
            counts = np.array([ic[cfg][ebin][comp] for comp in compList])
            ntot = float(counts.sum())
            bins = np.digitize(compList, cutoffs)
            for i, ptype in enumerate(groups):
                n_ptype = sum(counts[bins==i])
                row += [n_ptype/ntot]
            mytable += [row]

    # Repeat for IceTop
    row = ['IT73', '8-100']
    compList = np.array(it.keys())
    counts = np.array([it[comp] for comp in compList])
    ntot = float(counts.sum())
    bins = np.digitize(compList, cutoffs)
    for i, ptype in enumerate(groups):
        n_ptype = sum(counts[bins==i])
        row += [n_ptype/ntot]
    mytable += [row]

    return mytable


def paperTable(mytable):

    # Reduce to IC86
    data = [r[1:] for r in mytable if 'IC86' in r[0]]
    ebins = sorted(list(set([r[0] for r in data])))

    newtable = []
    for ebin in ebins:
        row = []
        emin, emax = ebin.split('-')
        median, sigL, sigR = readDist_IC('IC', emin, emax)
        row += [median]
        values = [r[1:] for r in data if r[0]==ebin]
        values = np.mean(values, axis=0)
        row += list(values)

        newtable += [row]

    # IceTop bin
    it = mytable[-1][1:]
    nmin, nmax = it[0].split('-')
    median, sigL, sigR = readDist_IT('IT', nmin, nmax)
    row = [median] + it[1:]
    newtable += [row]

    # Sort by energy bin
    newtable = sorted(newtable, key=lambda r: r[0])

    return newtable


def plotter(mytable):

    comps = ['P','He','CNO','Fe']
    for i, comp in enumerate(comps):
        p = {}
        configs = sorted(list(set([r[0] for r in mytable])))
        for config in configs:
            subtable = [row[1:] for row in mytable if row[0]==config]
            ebins = [r[0] for r in subtable]
            data = np.transpose([r[1:] for r in subtable])
            p[config] = [data[i][j] for j in range(len(ebins))]

        # Get x-values from bins
        ebins = np.array([ebin.split('-') for ebin in ebins])
        ebins = ebins.flatten().astype('float')
        ebins = sorted(list(set(ebins)))
        x = getMids(ebins, infvalue=100)

        fig, ax = plt.subplots()
        ax.set_title(comp)
        ax.set_xlabel('Reconstructed Energy')
        ax.set_ylabel('%')
        for config in configs:
            ax.plot(x, p[config], '.', label=config)
        ax.legend()
        plt.show()



if __name__ == "__main__":

    p = argparse.ArgumentParser()
    p.add_argument('-n', '--new', dest='new',
            default=False, action='store_true',
            help='Recalculate compTable.npy')
    p.add_argument('-o', '--out', dest='out',
            default='/home/fmcnally/anisotropy/icesim/compTable.npy',
            help='File name where composition information is stored')
    p.add_argument('--full', dest='full',
            default=False, action='store_true',
            help='Show full table including composition by config and ebin')
    p.add_argument('--paper', dest='paper',
            default=False, action='store_true',
            help='Output table for paper')
    args = p.parse_args()

    if args.new:
        yn = raw_input('This operation will overwrite the existing composition tables. Are you sure you want to continue? [y|n]: ')
        if yn == 'y':
            os.remove(args.out)
            makeTable(args.out)

    d = np.load(args.out)
    d = d.item()

    mytable = fullTable(d)

    if args.full:
        header = ['Config','Ebin'] + groups
        print tabulate(mytable, headers=header, floatfmt='.2f')

    if args.paper:
        newtable = paperTable(mytable)
        header = ['$E_{median}$', 'P','He','CNO','Fe']
        print
        print tabulate(newtable, headers=header, 
                floatfmt='.2f', tablefmt='latex')
    #plotter(mytable)
