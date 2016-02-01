#!/usr/bin/env python

###########################################################################
# Takes a list of processed hdf5 files and returns desired information    #
# as a binary dictionary for rapid reading for use with plotting          #
# Usage: Save_Sim.py [IT59|IT73]                                          #
###########################################################################

import sys, time, tables, pickle, glob
from numpy import *

if __name__ == "__main__":

    d = {}

    if len(sys.argv) != 2:
        print 'Incorrect usage: Save_Sim.py [config]'
        raise
    config = sys.argv[1]    # 'IT59' or 'IT73'

    prefix = '/net/user/fmcnally/ShowerLLH/%s_sim/' % config
    simList = {}
    simList['P'] = ['7006','7351','7579']
    simList['H'] = ['7483','7241','7263','7791']
    simList['O'] = ['7486','7242','7262','7851']
    simList['F'] = ['7007','7394','7784']
    # Only have a low energy sample for iron and proton, better to remove it?
    eList = {}
    eList['lowE']  = ['7351','7394','7483','7486']
    eList['medE']  = ['7006','7007','7241','7242','7262','7263']
    eList['highE'] = ['7579','7784','7791','7851']

    # Lists of keys to import
    values = ['RecoRanCut', 'zenith']
    tList = ['P', 'H', 'O', 'F']
    rList = ['p', 'h', 'o', 'f']
    values += tList
    for r in rList:
        values += ['maxLLH_'+r]
    for key in ['x', 'y', 'energy']:
        values += ['pML_'+key, 'hML_'+key, 'oML_'+key, 'fML_'+key]
    for key in ['x', 'y', 'energy', 'zenith']:
        values += ['MC_'+key]
    for key in ['x', 'y', 'fit_status', 's125', 'log10_s125_err', 'beta']:
        #values += ['comb_'+key, 'plane_'+key]
        values += ['lap_'+key]
    values += ['lap_pEnergy', 'lap_fEnergy']
    values += ['LoudestOnEdge', 'Q1', 'Q2', 'Q3', 'Q4']
    values += ['lowE', 'medE', 'highE']

    # Option to output likelihood space - makes files very large
    #llhs = []
    #for grid in ['Coarse', 'Med', 'Fine']:
    #    llhs += ['p'+grid+'_llhs', 'f'+grid+'_llhs']

    fileList = glob.glob(prefix + 'files/SimLLH_????.hdf5')
    fileList.sort()

    t0 = time.time()

    # Initialize arrays
    for value in values:
        d[value] = array([])

    for fl in fileList:

        print "Working on", fl
        sim = fl[-9:-5]
        t = tables.openFile(fl)
        temp = {}

        # Get ShowerLLH cuts and info
        temp['RecoRanCut'] = t.root.ShowerLLHParams.col('RecoRanCut')
        temp['zenith'] = t.root.ShowerPlane.col('zenith')
        temp['maxLLH_p'] = t.root.ShowerLLHParams.col('maxLLH_proton')
        temp['maxLLH_h'] = t.root.ShowerLLHParams.col('maxLLH_helium')
        temp['maxLLH_o'] = t.root.ShowerLLHParams.col('maxLLH_oxygen')
        temp['maxLLH_f'] = t.root.ShowerLLHParams.col('maxLLH_iron')

        # Get position information
        for value in ['x', 'y', 'energy']:
            temp['pML_'+value] = t.root.ShowerLLH_proton.col(value)
            temp['hML_'+value] = t.root.ShowerLLH_helium.col(value)
            temp['oML_'+value] = t.root.ShowerLLH_oxygen.col(value)
            temp['fML_'+value] = t.root.ShowerLLH_iron.col(value)

        # MCPrimary information
        for value in ['x', 'y', 'energy', 'zenith']:
            temp['MC_'+value] = t.root.MCPrimary.col(value)

        # Composition information
        for e in tList:
            test = sim in simList[e]
            temp[e] = array([test for i in range(len(temp['MC_x']))])

        # Low/Med/High energy simulation info
        for erange in ['lowE', 'medE', 'highE']:
            test = sim in eList[erange]
            temp[erange] = array([test for i in range(len(temp['MC_x']))])

        # Get other reconstruction information
        for value in ['x', 'y', 'fit_status']:
            temp['lap_'+value] = t.root.Laputop.col(value)
        for value in ['beta', 's125', 'log10_s125_err']:
            temp['lap_'+value] = t.root.LaputopParams.col(value)
        temp['lap_pEnergy'] = t.root.LaputopParams.col('e_proton')
        temp['lap_fEnergy'] = t.root.LaputopParams.col('e_iron')

        # Get Bakhtiyar's cuts
        temp['LoudestOnEdge'] = t.root.LoudestOnEdge.col('value')
        temp['Q1'] = t.root.Q1.col('value')
        temp['Q2'] = t.root.Q2.col('value')
        temp['Q3'] = t.root.Q3.col('value')
        temp['Q4'] = t.root.Q4.col('value')

        # append to existing arrays
        new = len(temp['RecoRanCut'])
        old = len(d['RecoRanCut'])
        for value in values:
            d[value].resize(old+new)
            d[value][old:] = temp[value]

        t.close()

    # Create final elements
    full_llhs = array([d['maxLLH_'+r] for r in rList])
    max_llh = amax(full_llhs, axis=0)
    for r in rList:
        d[r] = (d['maxLLH_'+r] == max_llh)
    for value in ['x', 'y', 'energy']:
        d['ML_'+value] = sum([d[r+'ML_'+value]*d[r] for r in rList], axis=0)

    print "Time taken:", time.time()-t0
    print "Average time per run:", (time.time()-t0)/len(fileList)

    fl = open(prefix+'SimPlot_test.pkl', 'wb')
    pickle.dump(d, fl)
    fl.close()
