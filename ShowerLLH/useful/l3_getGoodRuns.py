#!/usr/bin/env python

def fileCleaner(runFile, fileList):

    f = open(runFile, 'r')
    lines = f.readlines()
    goodRunList = []

    for line in lines:

        # Lines with run info need at least 6 parts
        line_info = line.split()
        #print line_info[0]
        #if len(line_info) < 6:
        #    continue
        # Make sure we're looking at a run line
        try:
            test = int(line_info[0])
        except ValueError:
            continue

        # Make sure run is good
        run  = '00'+line_info[0]
        #print run
        #IT   = line_info[3]
        #GOOD = line_info[5]

        #if (IT=='IT') and (GOOD=='GOOD'):
        goodRunList.append(run)
        #print goodRunList
    f.close()
    goodRunList.sort()
    cleanList = []
    for file in fileList:
        run_st = file.find('run') + 3
        run = file[run_st:run_st+8]
        if run in goodRunList:
            cleanList.append(file)

    return cleanList


