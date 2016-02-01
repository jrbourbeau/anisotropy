#!/usr/bin/env python

import glob, os

def maker(njobs=4, config='IT73', yyyymm=False):

    cwd = os.getcwd()
    os.chdir('/home/zgriffith/ShowerLLH/npx4')

    # Split into quarters
    prefix = '/net/user/zgriffith/ShowerLLH/%s_data/' % config
    files = glob.glob(prefix + 'DataPlot_??????.npy')
    monthList = [f[-10:-4] for f in files]
    monthList.sort()
    #n = len(monthList)/njobs
    #n = n if len(monthList)%njobs==0 else n+1
    #monthList = [[monthList[i+n*j] for i in range(n)] for j in range(njobs)]
    #print monthList
    #if yyyymm:
    #    monthList = [yyyymm]

    for ymList in monthList:
        #ymList = ' '.join(ymList)
        ex  = '/net/user/fmcnally/offline/V04-05-00/build/env-shell.sh'
        ex += ' python %s/make_maps.py %s %s' % (cwd, config, ymList)
        os.system('./submit_npx4.sh ' + ex)

    os.chdir(cwd)


if __name__ == "__main__":

    maker(config='IT73')
    #maker(config='IT81')
