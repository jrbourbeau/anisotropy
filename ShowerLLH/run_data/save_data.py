#!/usr/bin/env python

###########################################################################
# Takes a list of processed hdf5 files and returns desired information    #
# as a binary dictionary for rapid reading for use with plotting          #
# Usage: save_data.py                                                     #
###########################################################################

from icecube import astro
import sys, time, tables
from numpy import *
from copy import deepcopy

def saver(config, yyyymm):

    ##=======================================================================
    ## Starting parameters

    prefix = '/net/user/zgriffith/ShowerLLH/%s_data/' % config
    inFile   = prefix + 'files/DataLLH_'+yyyymm+'.hdf5'
    fileList = [inFile]
    outFile  = prefix + 'DataPlot_'+yyyymm+'.npy'

    # Find start and end time for the month
    t1 = astro.Time()
    yy, mm = int(yyyymm[:4]), int(yyyymm[4:])
    t1.SetTime(yy, mm, 1, 0, 0, 0)
    mjd0 = t1.GetMJD()
    yy = yy+1 if mm==12 else yy
    mm = mm+1 if mm!=12 else 1
    t1.SetTime(yy, mm, 1, 0, 0, 0)
    mjd1 = t1.GetMJD()

    # Basic composition information
    comp = {'p':'proton', 'h':'helium', 'n':'nitrogen', 'o':'oxygen'}
    comp.update({'a':'aluminum','f':'iron'})
    # Dictionary of filters for a config and their respective weights
    filterDict = {}
    for c in ['IT73','IT81','IT81-II']:
        filterDict[c] = {}
    filterDict['IT73']['IceTopSTA3_10'] = 'w3'
    filterDict['IT73']['IceTopSTA3_InIceSMT_10'] = 'w1.5'
    filterDict['IT73']['IceTopSTA8_10'] = 'w1'
    filterDict['IT81']['IceTopSTA3_11'] = 'w3'
    filterDict['IT81']['IceTopSTA3_InIceSMT_11'] = 'w1.5'
    filterDict['IT81']['IceTopSTA8_11'] = 'w1'
    filterDict['IT81-II']['IceTopSTA3_12'] = 'w10'
    filterDict['IT81-II']['IceTopSTA5_12'] = 'w1'
    f = filterDict[config]
    rList = ['p','f']

    ##=======================================================================
    ## Keys to import

    d = {}
    arrays = []
    dicts  = {}

    # General event information
    arrays += ['SubEventStream', 'mjd', 'NStations']
    dicts['weights'] = f.values()

    # ShowerLLH values
    arrays += ['rrc']
    for r in rList:
        arrays += [r+'LLH', r+'ML_energy', r+'ML_x', r+'ML_y', r]
    arrays += ['ML_x', 'ML_y', 'ML_energy']

    # Laputop values
    for key in ['x','y','zenith','azimuth','s125','e_proton','e_iron','beta']:
        arrays += ['lap_'+key]

    # Other reconstructions and cuts
    for key in ['zenith', 'azimuth']:
        arrays += ['ShowerPlane_'+key]

    # Initialize arrays
    for value in arrays:
        d[value] = array([])
    for key in dicts.keys():
        d[key] = {}
        for value in dicts[key]:
            d[key][value] = array([])

    # Create separate dicts for events that occur in past or next months
    d_prev = deepcopy(d)
    d_post = deepcopy(d)

    t0 = time.time()
    for file in fileList:

        print "Working on", file
        t = tables.openFile(file)
        q = {}

        # General event information
        q['SubEventStream'] = t.root.I3EventHeader.col('SubEventStream')
        q['NStations'] = t.root.NStations.col('value')
        mjd_day = t.root.I3EventHeader.col('time_start_mjd_day')
        mjd_sec = t.root.I3EventHeader.col('time_start_mjd_sec')
        mjd_ns  = t.root.I3EventHeader.col('time_start_mjd_ns')
        q['mjd'] = zeros(len(mjd_day), dtype=float64)
        for i in range(len(mjd_day)):
            day = int(mjd_day[i])
            sec = int(mjd_sec[i])
            ns  = int(mjd_ns[i])
            t1.SetTime(day, sec, ns)
            q['mjd'][i] = t1.GetMJD()

        # Condition and prescale passed (filter[condition, prescale])
        # For notes on weights see bottom of file
        for filter in f.keys():
            q[f[filter]] = t.root.FilterMask.col(filter)
            q[f[filter]] = q[f[filter]][:,0]
        if config in ['IT73','IT81']:
            q['w1.5'] *= logical_not(q['w1'])
            q['w3'] *= logical_not(q['w1']) * logical_not(q['w1.5'])
        elif config == 'IT81-II':
            q['w10'] *= logical_not(q['w1'])

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

        # Get other reconstruction info
        for value in ['zenith', 'azimuth']:
            q['ShowerPlane_'+value] = t.root.ShowerPlane.col(value)

        # Create final elements
        full_llhs = array([q[r+'LLH'] for r in rList])
        max_llh = amax(full_llhs, axis=0)
        for r in rList:
            q[r] = (q[r+'LLH'] == max_llh)
        for value in ['x', 'y', 'energy']:
            q['ML_'+value] = sum([q[r+'ML_'+value]*q[r] for r in rList], axis=0)

        # Separate into actual months
        mjd_cut  = (q['mjd']>=mjd0) * (q['mjd']<mjd1)
        prev_cut = (q['mjd'] <  mjd0)
        post_cut = (q['mjd'] >= mjd1)

        # Append to existing arrays
        def dictUpdate(q, d, cut):
            new = len(q['rrc'][cut])
            old = len(d['rrc'])
            for value in arrays:
                d[value].resize(old+new)
                d[value][old:] = q[value][cut]
            for key in dicts.keys():
                for value in dicts[key]:
                    d[key][value].resize(old+new)
                    d[key][value][old:] = q[value][cut]

        if mjd_cut.sum() != 0:
            dictUpdate(q, d, mjd_cut)
        if prev_cut.sum() != 0:
            dictUpdate(q, d_prev, prev_cut)
        if post_cut.sum() != 0:
            dictUpdate(q, d_post, post_cut)

        t.close()


    print "Time taken:", time.time()-t0
    print "Average time per run:", (time.time()-t0)/len(fileList)

    save(outFile, d)
    if len(d_prev['rrc'])!=0:
        prevFile = outFile.replace('.npy', '_prev.npy')
        save(prevFile, d_prev)
    if len(d_post['rrc'])!=0:
        postFile = outFile.replace('.npy', '_post.npy')
        save(postFile, d_post)



if __name__ == "__main__":

    config = sys.argv[1]
    yyyymm = sys.argv[2]
    #config = 'IT73'
    #for yyyymm in ['201006','201007','201008','201009','201010','201011','201012','201101','201102','201103','201104','201105']:

    saver(config, yyyymm)


###############################################################################
## Notes on weights

"""
 - Events that pass STA8 condition have a prescale and weight of 1.
 - Events that pass STA3ii condition have a 1/2 chance to pass the STA3ii
prescale. Those that fail have a 1/3 chance to pass the STA3 prescale. So, the 
total chance of an event passing is 1/2+(1/3*1/2) = 2/3. Weight = 1.5
 - Events that pass STA3 condition but not STA3ii condition have a prescale and
weight of 3
"""
