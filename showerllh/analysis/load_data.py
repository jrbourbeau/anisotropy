#!/usr/bin/env python

import numpy as np
import tables, time, os, datetime
from dateutil.relativedelta import relativedelta

import myGlobals as my
import dataFunctions as df
import llh_cut
from llhtools import inPoly
from skypos import getDecRA


def load_data(config, yyyymm, bintype='logdist'):

    my.setupShowerLLH(verbose=False)

    # Load data from npy file 
    t0 = time.time()
    inFile = '%s/%s_data/DataPlot_%s_%s.npy' % \
            (my.llh_data, config, yyyymm, bintype)
    #if not os.path.isfile(inFile):
    #    continue

    print 'Loading', inFile
    d = np.load(inFile)
    d = d.item()
    print 'Time taken:', time.time()-t0

    # Calculate sky positions
    d['dec'], d['ra'] = getDecRA(d)

    ## Quality Cuts ##
    d['cuts'] = {}

    # ShowerLLH cuts
    it_geo = df.it_geo(config)
    d['cuts']['llh'] = inPoly(d['ML_x'], d['ML_y'], 0, config=it_geo)
    d['cuts']['llha'] = inPoly(d['ML_x'], d['ML_y'], -50, config=it_geo)
    d['cuts']['llht'] = inPoly(d['ML_x'], d['ML_y'], -90, config=it_geo)

    # Likelihood cuts
    #llhDict = llh_cut.getLLHcut()
    #r = np.log10(d['ML_energy'])
    #pf_ratio = d['pLLH'] - d['fLLH']
    # Make so all events outside energy range fail
    #llhDict['pLLH'] = np.append(llhDict['pLLH'], np.inf)
    #llhDict['fLLH'] = np.append(llhDict['fLLH'], -np.inf)
    #e_bins = np.digitize(r, llhDict['Ebins']) - 1
    #pllh_min = llhDict['pLLH'][e_bins]
    #fllh_min = llhDict['fLLH'][e_bins]
    #d['cuts']['pLLH'] = (pf_ratio > pllh_min)
    #d['cuts']['fLLH'] = (pf_ratio < fllh_min)

    # Laputop cuts
    #d['cuts']['lap1'] = (np.cos(d['lap_zenith']) >= 0.8)
    #d['cuts']['lap2'] = inPoly(d['lap_x'], d['lap_y'], -90)
    #d['cuts']['lapE'] = (d['lap_e_proton']==d['lap_e_proton'])
    #d['cuts']['lapE']*= (d['lap_e_proton']!=0)
    #d['cuts']['lap']  = d['cuts']['lap1']*d['cuts']['lap2']*d['cuts']['lapE']

    return d


if __name__ == "__main__":

    p = argparse.ArgumentParser(description='Load ShowerLLH data from file')
    p.add_argument('-c', '--config', dest='config',
            help='Detector configuration [IT59 | IT73 | IT81]')
    p.add_argument('-d', '--date', dest='date',
            help='Date to look at (yyyymm)')
    args = p.parse_args()

    d = load_data(args.config, args.date)


