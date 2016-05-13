#!/usr/bin/env python

import os, glob, argparse

if __name__ == "__main__":

    p = argparse.ArgumentParser(description='Merges simple-dst files')
    p.add_argument('-c', '--config', dest='config',
            default='IC86-IV',
            help='Detector configuration to run over')
    p.add_argument('--overwrite', dest='overwrite',
            default=False, action='store_true', 
            help='Option to overwrite existing output files')
    args = p.parse_args()

    if args.config == 'IC86-IV':
        outdir = '/data/ana/CosmicRay/Anisotropy/IceCube/IC86/2014/simple-dst'

    # Get list of all part files
    fileList = glob.glob('%s/*_p??.root' % outdir)
    fileList.sort()
    # Reduce to list of filebases, ignoring part number
    baseList = sorted(list(set([f[:f.rfind('_')] for f in fileList])))

    cwd = os.getcwd()
    os.chdir(outdir)

    for filebase in baseList:

        # Check for (and deal with) existing outfile
        out = '%s.root' % filebase
        if os.path.isfile(out):
            if not args.overwrite:
                print 'File %s already exists!' % out
                continue
            os.remove(out)

        tempList = sorted(glob.glob('%s_p??.root' % filebase))
        nfiles = len(tempList)
        tempList = ' '.join(tempList)
        ex = ' '.join(['hadd', out, tempList])
        outBase = os.path.basename(filebase)
        os.system(ex)

    os.chdir(cwd)
