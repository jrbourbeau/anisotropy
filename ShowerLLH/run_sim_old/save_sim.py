#!/usr/bin/env python

###########################################################################
# Takes a list of processed hdf5 files and returns desired information    #
# as a binary dictionary for rapid reading for use with plotting          #
# Usage: save_sim.py [config]                                             #
###########################################################################

import sys, time, tables, glob
from numpy import *

if __name__ == "__main__":

    d = {}
    config = 'IT73'
    prefix = '/net/user/fmcnally/ShowerLLH/%s_sim/' % config
    fileList = glob.glob(prefix + 'files/SimLLH_????.hdf5')
    fileList.sort()

    # Basic composition and simulation information
    comp = {'p':'proton', 'h':'helium', 'n':'nitrogen', 'o':'oxygen'}
    comp.update({'a':'aluminum','f':'iron'})
    typeDict = {'none':0,'P':14,'He':402,'N':1407,'O':1608,'Al':2713,'Fe':5626}
    eList = {}
    eList['lowE']  = ['7351','7394','7483','7486']
    eList['medE']  = ['7006','7007','7241','7242','7262','7263']
    eList['highE'] = ['7579','7784','7791','7851']
    #tList = ['P', 'He', 'O', 'Fe']
    tList = ['P', 'Fe']
    rList = [t[0].lower() for t in tList]

    ## Keys to import ##
    values  = []
    values += tList

    # ShowerLLH values
    values += ['rrc']
    for r in rList:
        values += [r+'LLH', r+'ML_energy', r+'ML_x', r+'ML_y']

    # Laputop values
    for key in ['x','y','zenith','azimuth','s125','e_proton','e_iron','beta']:
        values += ['lap_'+key]

    # MC true values
    for key in ['x', 'y', 'energy', 'zenith']:
        values += ['MC_'+key]

    # Other reconstructions and cuts
    for key in ['zenith', 'azimuth']:
        values += ['ShowerPlane_'+key]
    values += ['LoudestOnEdge', 'Q1', 'Q2', 'Q3', 'Q4']
    values += ['lowE', 'medE', 'highE']

    # Option to output likelihood space - makes files very large
    #llhs = []
    #for grid in ['Coarse', 'Med', 'Fine']:
    #    llhs += ['p'+grid+'_llhs', 'f'+grid+'_llhs']

    # Initialize arrays
    for value in values:
        d[value] = array([])

    t0 = time.time()
    for file in fileList:

        print "Working on", file
        st = file.find('SimLLH_') + 7
        sim = file[st:st+4]
        t = tables.openFile(file)
        q = {}

        # Get ShowerLLH cuts and info
        q['rrc'] = t.root.ShowerLLHParams.col('RecoRanCut')
        for value in ['x', 'y', 'energy']:
            q['pML_'+value] = t.root.ShowerLLH_proton.col(value)
            #q['hML_'+value] = t.root.ShowerLLH_helium.col(value)
            #q['oML_'+value] = t.root.ShowerLLH_oxygen.col(value)
            q['fML_'+value] = t.root.ShowerLLH_iron.col(value)
        for r in rList:
            q[r+'LLH'] = t.root.ShowerLLHParams.col('maxLLH_'+comp[r])

        # Get Laputop info
        for value in ['x', 'y', 'zenith', 'azimuth']:
            q['lap_'+value] = t.root.Laputop.col(value)
        for value in ['s125', 'e_proton', 'e_iron', 'beta']:
            q['lap_'+value] = t.root.LaputopParams.col(value)

        # Other reconstruction info
        for value in ['zenith', 'azimuth']:
            q['ShowerPlane_'+value] = t.root.ShowerPlane.col(value)

        # MCPrimary information
        for value in ['x', 'y', 'energy', 'zenith', 'azimuth', 'type']:
            q['MC_'+value] = t.root.MCPrimary.col(value)
        for true_comp in typeDict.keys():
            q[true_comp] = (q['MC_type'] == typeDict[true_comp])

        # Low/Med/High energy simulation info
        for erange in ['lowE', 'medE', 'highE']:
            test = sim in eList[erange]
            q[erange] = array([test for i in range(len(q['MC_x']))])

        # Get cuts
        q['LoudestOnEdge'] = t.root.LoudestOnEdge.col('value')
        q['Q1'] = t.root.Q1.col('value')
        q['Q2'] = t.root.Q2.col('value')
        q['Q3'] = t.root.Q3.col('value')
        q['Q4'] = t.root.Q4.col('value')

        # append to existing arrays
        new = len(q['rrc'])
        old = len(d['rrc'])
        for value in values:
            d[value].resize(old+new)
            d[value][old:] = q[value]

        t.close()

    # Create final elements
    full_llhs = array([d[r+'LLH'] for r in rList])
    max_llh = amax(full_llhs, axis=0)
    for r in rList:
        d[r] = (d[r+'LLH'] == max_llh)
    for value in ['x', 'y', 'energy']:
        d['ML_'+value] = sum([d[r+'ML_'+value]*d[r] for r in rList], axis=0)

    print "Time taken:", time.time()-t0
    print "Average time per run:", (time.time()-t0)/len(fileList)

    save(prefix + 'SimPlot.npy', d)
