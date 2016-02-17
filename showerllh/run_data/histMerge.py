#!/usr/bin/env python

import numpy as np
import glob, os, argparse

import myGlobals as my

if __name__ == "__main__":

    # Setup global paths
    my.setupShowerLLH(verbose=False)

    p = argparse.ArgumentParser(description='Creates merged histogram files')
    p.add_argument('-c', '--configs', dest='configs', nargs='*',
            default=['IT73','IT81','IT81-II','IT81-III','IT81-IV'],
            help='Detector configuration(s) to run over')
    p.add_argument('-b', '--bintype', dest='bintype',
            default='logdist',
            help='ShowerLLH bin type')
    p.add_argument('--sky', dest='sky',
            default=False, action='store_true',
            help='Parameters based on sky position')
    args = p.parse_args()

    for config in args.configs:

        prefix = '%s/%s_data' % (my.llh_data, config)
        histFiles = glob.glob('%s/hists/*_%s.npy' % (prefix, args.bintype))
        if args.sky:
            histFiles = glob.glob('%s/hists/*_%s_sky.npy' % \
                    (prefix, args.bintype))
        histFiles.sort()

        # Build list of months
        months = [os.path.basename(f).split('_')[1][:6] for f in histFiles]
        months = sorted(list(set(months)))

        htotal = {}
        for month in months:

            if args.sky:
                continue

            print 'Working on %s %s...' % (config, month)
            mFiles = glob.glob('%s/hists/DataLLH_%s*_%s.npy' % \
                    (prefix, month, args.bintype))
            mFiles.sort()

            # Combine histograms
            hmonth = {}
            for i, histFile in enumerate(mFiles):
                h = np.load(histFile)
                h = h.item()
                if i == 0:
                    hmonth = h
                else:
                    for key in h.keys():
                        hmonth[key] += h[key]

            for key in hmonth.keys():
                tkey = '%s_%s_%s' % (config, month, key)
                htotal[tkey] = hmonth[key]

        if args.sky:
            hcfg = {}
            for i, histFile in enumerate(histFiles):
                print 'Working on %s...' % histFile
                h = np.load(histFile)
                h = h.item()
                if i == 0:
                    hcfg = h
                else:
                    for key in h.keys():
                        hcfg[key] += h[key]

            for key in hcfg.keys():
                tkey = '%s_%s' % (config, key)
                htotal[tkey] = hcfg[key]

        # Write to file
        outfile = '%s/%s_hists_%s.npy' % (prefix, config, args.bintype)
        if args.sky:
            outfile = '%s/%s_hists_%s_sky.npy' % (prefix, config, args.bintype)

        np.save(outfile, htotal)


