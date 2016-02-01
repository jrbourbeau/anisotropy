#!/usr/bin/env python
import sys, glob

if __name__ == "__main__":

    yyyymm = sys.argv[1]
    n = int(sys.argv[2])

    # Create combined list of IT and GCD files
    yy = yyyymm[:4]
    mm = yyyymm[4:]
    prefix = '/data/exp/IceCube/'+yy+'/filtered/level2a/'
    itList  = glob.glob(prefix + mm+'*/*_IT*')
    gcdList = glob.glob(prefix + mm+'*/*GCD*')
    masterList = itList + gcdList
    masterList.sort()

    # Get good run list
    goodList = []
    f = open('/net/user/fmcnally/ShowerLLH/resources/IT73_GoodRuns.txt', 'r')
    for line in f:
        goodList.append(line[8:16])
    f.close()

    # Clean the master list, through the use of a bad list
    badList = []
    for item in masterList:
        run = item[66:74]
        if run not in goodList:
            badList.append(item)
    for item in badList:
        masterList.remove(item)

    # Create individual job lists
    ntot = len(masterList)
    if ntot%n != 0:
        print ntot/n + 1
    else:
        print ntot/n
