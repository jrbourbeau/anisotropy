#!/usr/bin/env python

import glob, re, os, argparse
import myGlobals as my


def checker(config, verbose=True):

    # Get all data files with run and part numbers
    pattern = re.compile('.*DataLLH_\d{8}_.*_\d{8}_\d{8}_\d{8}_\d{8}\.hdf5')
    prefix = '%s/%s_data/files' % (my.llh_data, config)
    fileList = glob.glob('%s/DataLLH_*.hdf5' % prefix)
    fileList = [f for f in fileList if pattern.match(f)]
    fileList.sort()

    basenames = [os.path.splitext(os.path.basename(f))[0] for f in fileList]
    dates = sorted(list(set([n.split('_')[1] for n in basenames])))

    if verbose:
        print 'May be a disjoint between the following files...'
    missing = []
    for date in dates:

        tempBases = [n for n in basenames if '_%s_' % date in n]
        starts = [n.split('_')[-3] for n in tempBases]
        ends   = [n.split('_')[-1] for n in tempBases]

        for i, (start, end) in enumerate(zip(starts[1:], ends[:-1])):
            if (int(start) != int(end)+1) and (int(start) != 0):
                if verbose:
                    print '--------------------------------------'
                    print tempBases[i]
                    print tempBases[i+1]
                # Create missing filename
                params = tempBases[i].split('_')
                params[-4:-2] = tempBases[i].split('_')[-2:]
                params[-2:] = tempBases[i+1].split('_')[-4:-2]
                params[-3] = '%08i' % (int(params[-3])+1)
                params[-1] = '%08i' % (int(params[-1])-1)
                if verbose:
                    print 'Suggested out:'
                    print '_'.join(params)
                missing += ['%s/%s.hdf5' % (prefix, '_'.join(params))]

    return missing


if __name__ == "__main__":

    my.setupShowerLLH(verbose=False)
    p = argparse.ArgumentParser(description='Searches for dropped npx4 jobs')
    p.add_argument('-c', '--config',
            default='IT81',
            help='Detector configuration')
    args = p.parse_args()

    checker(args.config)
