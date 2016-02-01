#!/usr/bin/env python

############################################################################
# Makes the bins for the LLH tables
############################################################################

from numpy import *
import sys

def makeBins(resourcedir):

    d = {}

    # Charge bins in VEM
    # Original Cbins:
    #d['Cbins'] = append(array([0]), 10**linspace(-3,4.5,46))
    #d['Cbins'] = append(d['Cbins'], inf)
    # modify binning - multiple bins on either end were empty
    d['Cbins'] = append(array([0]), 10**linspace(-3,4.5,46)[10:41])
    d['Cbins'] = append(d['Cbins'], inf)

    # Energy bins in Log10(Energy/GeV)
    d['Ebins'] = arange(4, 9.501, 0.05)

    # Distance bins in meters
    d['Dbins'] = append(arange(0,600,10), arange(600,1051,50))
    d['Dbins'] = append(d['Dbins'], inf)

    #Snow bins in meters
    d['Sbins'] = array([-1, .001, .5, .85])
    d['Sbins'] = append(d['Sbins'], inf)

    # Zenith Bins in radians (made with equal solid angle bins)
    sigma  = 2*pi*(1 - cos(40*pi/180))/4
    theta1 = arccos(1 - sigma/(2*pi))
    theta2 = arccos(cos(theta1) - sigma/(2*pi))
    theta3 = arccos(cos(theta2) - sigma/(2*pi))
    d['Zbins']  = array([0, theta1, theta2, theta3, pi/2])

    # Histogram shape will be (Energy, Zenith, Snow, Distance, Charge)
    outFile = '%s/ShowerLLH_bins.npy' % resourcedir
    save(outFile, d)


if __name__ == "__main__":

    resourcedir = '/net/user/zgriffith/ShowerLLH/resources'
    makeBins(resourcedir)



