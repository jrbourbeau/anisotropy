#!/usr/bin/env python

from numpy import *
import tables, sys, time, os
from useful import inPoly
import llh_cut
import datetime
from dateutil.relativedelta import relativedelta

def load_data(config, yyyymm):

    prefix = '/net/user/fmcnally/ShowerLLH/'+config+'_data'
    d = {}

    # Get previous and next month files
    yy = int(yyyymm[:4])
    mm = int(yyyymm[4:])
    last_month = datetime.date(yy, mm, 01) - relativedelta(months=1)
    ym_prev = last_month.strftime('%Y%m') + '_post'
    next_month = datetime.date(yy, mm, 01) + relativedelta(months=1)
    ym_post = next_month.strftime('%Y%m') + '_prev'

    ymList = [ym_prev, yyyymm, ym_post]
    filesFound = False
    for yyyymm in ymList:

        # Load data from npy file 
        t0 = time.time()
        inFile = '%s/DataPlot_%s.npy' % (prefix, yyyymm)
        if not os.path.isfile(inFile):
            continue

        print 'Loading', inFile
        q = load(inFile)
        q = q.item()
        for key in q.keys():
            if key not in d.keys():
                d[key] = q[key]
            elif isinstance(q[key], dict):
                for subkey in q[key].keys():
                    d[key][subkey] = append(d[key][subkey], q[key][subkey])
            else:
                d[key] = append(d[key], q[key])

        filesFound = True
        print 'Time taken:', time.time()-t0

    # Make sure some files were found
    if not filesFound:
        print 'No files found for', config, yyyymm
        sys.exit(1)

    # Set array types
    bools = ['p', 'f']
    for key in bools:
        d[key] = d[key].astype('bool')
    for weight in d['weights'].keys():
        d['weights'][weight] = d['weights'][weight].astype('bool')

    ## Quality Cuts ##
    d['cuts'] = {}

    # ShowerLLH cuts
    d['cuts']['llh1'] = (d['SubEventStream'] == 0)
    d['cuts']['llh2'] = inPoly(d['ML_x'], d['ML_y'], -50, config=config)
    d['cuts']['llh']  = d['cuts']['llh1'] * d['cuts']['llh2']

    # Likelihood cuts
    llhDict = llh_cut.getLLHcut()
    r = log10(d['ML_energy'])
    pf_ratio = d['pLLH'] - d['fLLH']
    # Make so all events outside energy range fail
    llhDict['pLLH'] = append(llhDict['pLLH'], inf)
    llhDict['fLLH'] = append(llhDict['fLLH'], -inf)
    e_bins = digitize(r, llhDict['Ebins']) - 1
    pllh_min = llhDict['pLLH'][e_bins]
    fllh_min = llhDict['fLLH'][e_bins]
    d['cuts']['pLLH'] = (pf_ratio > pllh_min)
    d['cuts']['fLLH'] = (pf_ratio < fllh_min)

    # Laputop cuts
    d['cuts']['lap1'] = (cos(d['lap_zenith']) >= 0.8)
    d['cuts']['lap2'] = inPoly(d['lap_x'], d['lap_y'], -90)
    d['cuts']['lapE'] = (d['lap_e_proton']==d['lap_e_proton'])
    d['cuts']['lapE']*= (d['lap_e_proton']!=0)
    d['cuts']['lap']  = d['cuts']['lap1']*d['cuts']['lap2']*d['cuts']['lapE']

    return d


if __name__ == "__main__":

    config = sys.argv[1]
    yyyymm = sys.argv[2]
    d = load_data(config, yyyymm)


