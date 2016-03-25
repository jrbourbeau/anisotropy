#!/usr/bin/env python

#==============================================================================
# File Name     : maker.py
# Description   : Plot fitmaps in batch
# Creation Date : 03-15-2016
# Last Modified : Thu 24 Mar 2016 08:49:51 AM CDT
# Created By    : James Bourbeau
#==============================================================================

import os

if __name__ == "__main__":

    initlist = [0.,1e-15,1e-12,1e-9,1e-6]
    steplist = [1e-6,1e-9,1e-12]
    # initlist = [0.,1e-15]
    # steplist = [1e-6,1e-9]
    for init in initlist:
        for step in steplist:
            # cwd = os.getcwd()
            # if not os.path.isdir(cwd+'/SHCoeff{}_{}'.format(init,step)):
            #     os.mkdir(cwd+'/SHCoeff{}_{}'.format(init,step))
            # chi2list = ['standard', 'RI', 'tibetonly']

            os.system('/home/jbourbeau/anisotropy/dev/phaseplot.py --step {} --init {}'.format(step,init))
            # for chi2 in chi2list:
            #     for l in range(1,8):
                    # os.system('submit_npx4 /home/jbourbeau/anisotropy/dev/almTibet.py --dipole -l {} --chi2 {} --step {} --init {}'.format(l,chi2,step,init))
                    # os.system('/home/jbourbeau/anisotropy/mapfunctions/plotFITS.py -n single -v --scale 3 -f fitmaps/fitmap_lmax{}_{}chi2_dipole2.fits -d -90 -D -25 --outDir /home/jbourbeau/public_html/figures/2Dfit/dipolemaps/ --customOut fitmap_{}chi2_lmax{}_dipole -o -b'.format(l,chi2,chi2,l))
