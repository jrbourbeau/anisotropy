#!/usr/bin/env python

#=============================================================================
# File Name     : maker.py
# Description   : Plot fitmaps in batch
# Creation Date : 03-15-2016
# Last Modified : Thu 24 Mar 2016 08:49:51 AM CDT
# Created By    : James Bourbeau
#=============================================================================

import os

if __name__ == "__main__":

    # initlist = [0.,1e-15,1e-12,1e-9,1e-6]
    # steplist = [1e-6,1e-9,1e-12]
    # initlist = [0.,1e-15]
    # steplist = [1e-6,1e-9]
    # for init in initlist:
    #     for step in steplist:
            # cwd = os.getcwd()
            # if not os.path.isdir(cwd+'/SHCoeff{}_{}'.format(init,step)):
            #     os.mkdir(cwd+'/SHCoeff{}_{}'.format(init,step))
    chi2list = ['standard', 'RI', 'tibetonly']

    # os.system('/home/jbourbeau/anisotropy/dev/phaseplot.py --step {} --init {}'.format(step,init))
    scalelist= [1.26,1.23,1.85,1.66,7.12]
    for chi2 in chi2list:
        for l in range(1,6):
            # os.system('/home/jbourbeau/anisotropy/dev/almTibet.py -l {} --chi2 {}'.format(l,chi2))
            # os.system('submit_npx4 /home/jbourbeau/anisotropy/dev/almTibet.py --dipole -l {} --chi2 {} --step {} --init {}'.format(l,chi2,step,init))
            os.system('/home/jbourbeau/anisotropy/mapfunctions/plotFITS.py -n single -v --scale 3 -f fitmapsOptimal/fitmap_lmax{}_{}chi2.fits -d -90 -D -25 --outDir /home/jbourbeau/public_html/figures/collaboration-meeting-2016/fitmapsOptimalScaled/ --customOut fitmap_{}chi2_lmax{} -o -b -m {} -M {}'.format(l,chi2,chi2,l,-scalelist[l-1],scalelist[l-1]))
