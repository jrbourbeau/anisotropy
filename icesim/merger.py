#!/usr/bin/env python

import glob, os, sys, re
import numpy as np
import myGlobals as my

def badFileMerger(overwrite=False):

    # Setup global path names
    my.setupAnisotropy(verbose=False)

    outFile = '%s/badFiles.txt' % my.ani_sim
    d = []
    if os.path.isfile(outFile):
        with open(outFile, 'r') as f:
            d = f.readlines()

    fileList = glob.glob('%s/*_badFiles.txt' % (my.ani_sim))
    fileList.sort()
    for file in fileList:
        with open(file, 'r') as f:
            tempLines = f.readlines()
        d += tempLines

    d = list(set(d))
    with open(outFile, 'w') as f:
        f.writelines(d)


def simMerger():

    # Setup global path names
    my.setupAnisotropy(verbose=False)

    fileList = glob.glob('%s/I*.npy' % my.ani_sim)
    fileList = [f for f in fileList if 'badFiles' not in f]
    fileList = [f for f in fileList if 'median' not in f]
    fileList.sort()

    simList = []
    for file in fileList:
        sim = os.path.basename(file)
        sim = '_'.join(re.split('_|\.', sim)[:2])
        if sim not in simList:
            simList += [sim]

    for sim in simList:

        d = {}

        fileList = glob.glob('%s/%s_*.npy' % (my.ani_sim, sim))
        fileList = [f for f in fileList if 'badFiles' not in f]
        fileList = [f for f in fileList if 'median' not in f]
        fileList.sort()

        if len(fileList) == 0:
            print 'No subfiles found for '+sim+'. Nothing to merge.'
            continue

        i = 0
        idx0, idx1 = 0, len(fileList)
        tempList = fileList[idx0:idx1]

        while len(tempList) > 1:

            fileSize = sum(np.array([os.stat(f).st_size for f in tempList]))
            if fileSize > 1e9:
                idx1 -= 1
                tempList = fileList[idx0:idx1]

            else:
                outFile = '%s/%s_Part%02d.npy' % (my.ani_sim, sim, i)
                if os.path.isfile(outFile):
                    print 'File %s already exists.' % outFile
                    print 'Skipping...'
                    continue
                print 'Outfile:', outFile

                for j, file in enumerate(tempList):

                    print 'Working on', file
                    # Load dictionary from npy file
                    temp = np.load(file)
                    temp = temp.item()
                    # Append to total dictionary
                    for key in temp.keys():
                        if j == 0:
                            d[key] = np.array([])
                        d[key] = np.append(d[key], temp[key])

                # Write to file
                np.save(outFile, d)

                # Move onto next set of files
                i += 1
                idx0 = idx1
                idx1 = len(fileList)
                tempList = fileList[idx0:idx1]


if __name__ == "__main__":

    simMerger()
    badFileMerger()

