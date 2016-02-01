#!/usr/bin/env python

from numpy import *
import glob, tables, pickle
#from matplotlib import nxutils

##=========================================================================##
## Setup useful global variables

n_prefix = '/net/user/zgriffith/ShowerLLH/'
t = tables.openFile('/net/user/santander/ShowerLLH/resources/ShowerLLH_bins.hdf5')
Ebins = t.root.bins.col('Ebins')[0]
nE = len(Ebins) - 1
t.close()

# For use with geomeric containment cut
fl = open('/home/santander/ShowerLLH/useful/tankpos.pkl', 'rb')
tank_x, tank_y = pickle.load(fl)
fl.close()


##=========================================================================##
## Functions for getting file lists for jobs

""" Filelist and name for simulation """
def simList(config, sim, group=None, startnum=None, testnum=None, table=False):

    gcdPre  = '/data/sim/sim-new/downloads/GCD_'

    if config == 'IT59':
        prefix = '/data/sim/IceTop/2009/filtered/level2/CORSIKA-ice-top/'
        gcd = gcdPre+'20_04_10/GeoCalibDetectorStatus_IC59.55040_official.i3.gz'
    elif config == 'IT73':
        prefix = '/data/sim/IceTop/2010/filtered/level2a/CORSIKA-ice-top/'
        gcd = gcdPre+'31_08_11/GeoCalibDetectorStatus_IC79.55380_L2a.i3.gz'
    elif config == 'IC79':
        prefix = '/data/sim/IceCube/2010/filtered/level2a/CORSIKA-in-ice/'
        gcd = gcdPre+'31_08_11/GeoCalibDetectorStatus_IC79.55380_L2a.i3.gz'

    name = n_prefix + '%s_sim/files/SimLLH_%s.hdf5' % (config, sim)
    i3name = n_prefix + '%s_sim/SimLLH_%s.i3' % (config, sim)
    if table:
        name = n_prefix + 'resources/%s/LLHTable_%s.pkl' % (config, sim)
    simList = glob.glob(prefix + sim + '/*/*')
    if group:
        simList = glob.glob(prefix + sim + '/'+group+'/*')
    simList.sort()
    fileList = [gcd] + simList
    if startnum and testnum:
        fileList = [gcd] + fileList[startnum: startnum+testnum]
        name = name.replace('.','_%05d_%05d.' % (startnum,startnum+testnum-1))
    if group:
        st  = int(group[:5])  + 1
        end = int(group[-5:]) + 1
        name = name.replace('.','_%05d_%05d.' % (st, end))
        i3name = i3name.replace('.','_%05d_%05d.' % (st, end))

    return fileList, name, i3name


""" Filelist and name for data """
def dataList(yyyymm, index, n=30):

    # Create combined list of IT and GCD files
    yyyy = yyyymm[:4]
    mm   = yyyymm[4:]
    dataPrefix = '/data/exp/IceCube:'+yyyy+'/filtered/level2a/'
    itList  = glob.glob(dataPrefix + mm+'*/*_IT*')
    gcdList = glob.glob(dataPrefix + mm+'*/*GCD*')
    masterList = itList + gcdList
    masterList.sort()

    # Get good run list
    goodList = []
    prefix = '/net/user/santander/ShowerLLH/'
    f = open(prefix + 'resources/IT73_GoodRuns.txt', 'r')
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
    jobList = []
    for i in range(ntot/n):
        jobList.append(masterList[n*i:n*(i+1)])
    if ntot%n != 0:
        jobList.append(masterList[n*(i+1):])

    # Clean our specific job list
    files = jobList[index]
    run = files[0][66:74]
    if 'GCD' not in files[0]:
        gcd = [x for x in gcdList if run in x]
        files.insert(0, gcd[0])
    if 'GCD' in files[len(files)-1]:
        files.pop()

    # Get part numbers for file name
    pstart = files[1][79:87]
    name = prefix + 'IT73_data/%s_%s_%s.hdf5' % (yyyymm, run, pstart)

    return files, name


##=========================================================================##
## Cuts for use on data                          

""" Vertex finder for use with geometric cut """
def getVertex(d, config='IT73'):

    # Define the sides by the tanks they contain, array built bottom up, L->R
    # Note: there appears to only be one tank showing at station 39
    side = {}
    if config == 'IT73':
        side[0] = arange(0,10)
        side[1] = array([8,9,20,21,34,35,50,51,67,68,87,88])
        side[2] = array([87,88,105,106,121,122,135,136])
        side[3] = arange(131,137)
        side[4] = array([131,132,143,144])
        side[5] = arange(137,145)
        side[6] = array([69,70,89,90,107,108,123,124,137,138])
        side[7] = array([10,11,22,23,36,37,52,53,69,70])

    # Find slopes and y-intercepts
    nside = len(side)
    m, b = zeros((2,nside))
    for i in range(nside):
        m[i], b[i] = polyfit(tank_x[side[i]], tank_y[side[i]], 1)

    # Intersection function
    def x_intersect(m1, b1, m2, b2):
        x = (b2-b1)/(m1-m2)
        return x

    # d = distance moved perpendicular from original line position
    theta = arctan(m)
    b1 = b + sign(b)*d/cos(theta)
    x = array([x_intersect(m[i-1], b1[i-1], m[i], b1[i]) for i in range(nside)])
    y = m*x + b1

    return x, y


""" Geometric containment cut """
def inPoly(xx, yy, d, config='IT73'):
    # values used for the original attempt at recreating Bakhtiyar's cuts
    #vertx = array([-310, 0, 60, 330, 530, 300, -100, -420])
    #verty = array([365, 405, 320, 360, 140, -350, -415, 10])
    vertx, verty = getVertex(d, config)
    xypoints = array([xx,yy]).transpose()
    xyverts  = array([vertx, verty]).transpose()
    tf = nxutils.points_inside_poly(xypoints, xyverts)
    return tf



##=========================================================================##
## Other useful functions

""" Convert bin to energy in log10(GeV) """
def b2e(bin):
    step = (Ebins[1]-Ebins[0]) / 2.0
    center = Ebins[bin] + step
    return 10**center

Emids = log10(b2e(range(nE)))


""" Return an array of the number of counts per energy bin """
def Nfinder(a, cut, w=False):
    if not w:
        Npassed, bins = histogram(a[cut], Ebins)
        return Npassed.astype('float')
    else:
        w1, w24, w8 = w
        Npassed_1, bins = histogram(a[cut*w1], Ebins)
        Npassed_24, bins = histogram(a[cut*w24], Ebins)
        Npassed_8, bins = histogram(a[cut*w8], Ebins)
        return (Npassed_1 + 2.4*Npassed_24 + 8*Npassed_8).astype('float')


""" Return an array of the error per energy bin """
def Efinder(a, cut, w=False):
    if not w:
        x = a[cut]
        n, bins = histogram(x, Ebins)
        return sqrt(n)
    else:
        w1, w24, w8 = w
        n1, bins  = histogram(a[cut*w1], Ebins)
        n24, bins = histogram(a[cut*w24], Ebins)
        n8, bins  = histogram(a[cut*w8], Ebins)
        return sqrt(n1 + 2.4**2 * n24 + 8**2 * n8)



#############################################################################
# Less useful and/or just used for plotting tests
#############################################################################


""" Simple smoothing function for unfolded spectrum """
def smoother(a):

    temp = zeros(len(a))
    for i in range(len(a)):
        if i == 0:
            temp[i] = (2*a[i] + a[i+1]) / 3
        elif i == len(a)-1:
            temp[i] = (2*a[i] + a[i-1]) / 3
        else:
            temp[i] = (a[i-1] + 2*a[i] + a[i+1]) / 4
    return temp


""" VEM to GeV conversion taken from Bakhtiyar's work with S125 """
def vem2gev(s125, comp, zenith):
    a = cos(zenith)
    if comp == 'P':
        if a >= 0.95:
            p = [5.9977, 0.9589, -0.001]
        elif a >= 0.9:
            p = [6.0340, 0.9455, -0.0007]
        elif a >= 0.85:
            p = [6.0823, 0.9288, -0.0013]
        elif a >= 0.8:
            p = [6.1397, 0.9186, -0.0013]
        else:
            p = [0, 0, 0]
    elif comp == 'F':
        if a >= 0.95:
            p = [6.0893, 0.8695, 0.0143]
        elif a >= 0.9:
            p = [6.1525, 0.8516, 0.0164]
        elif a >= 0.85:
            p = [6.2241, 0.8427, 0.0147]
        elif a >= 0.8:
            p = [6.3066, 0.8372, 0.0134]
        else:
            p = [0, 0, 0]

    logE = p[2]*log10(s125)**2 + p[1]*log10(s125) + p[0]
    return logE


""" Get the information for IT26 spectrum """
def getBakh():
    file = genfromtxt('it73.dat')
    logE = file[:,0]
    Nevents = file[:,1]
    flux = file[:,2]
    sigma = file[:,3]
    syserr_upper = file[:,4]
    syserr_lower = file[:,5]
    return [logE, Nevents, flux, sigma, syserr_upper, syserr_lower]

