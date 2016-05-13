#!/usr/bin/env python
################################################################################
# Parse a good/bad run list for IC86, get a list of DST files from the data
# warehouse, and extract an analysis-grade ntuple from each run.
#
# Then set up simple DST processing for each day by splitting the day's subruns
# into groups of files, and then merging the result.  The script will output a
# DAGMan file that can be executed by running
#
# condor_submit_dag -config dagman.config ntmaker.dag
################################################################################
#
# Version: $Id: ntmaker.py 120535 2014-06-09 21:47:04Z sybenzvi $
#
################################################################################

import itertools
import math
import os
import re
import sys
import subprocess
import argparse
from collections import defaultdict

# IceCube environment
envs = "/data/user/sybenzvi/icecube/offline-software/V14-03-01/build/env-shell.sh"

# Local disk output directory
outp = "${_CONDOR_SCRATCH_DIR}/simple-dst"

# Network disk submission script directory
pdir = "/data/user/sybenzvi/npx/ic86-03/dstcut"

# Network disk final destination directory
prfx = "/data/ana/IC86/CRAnisotropy/2014/simple-dst"

# -----------------------------------------------------------------------------
def splitSequence(iterable, size):
    """Take an iterable sequence of arbitrary length and split it into chunks
    of a specified size
    """
    it = iter(iterable)
    item = list(itertools.islice(it, size))
    while item:
        yield item
        item = list(itertools.islice(it, size))

# -----------------------------------------------------------------------------
def getRunDates():
    """Get a list of dates with extracted DSTs
    """
    years = ["2013", "2014"]
    yymmdd = []
    for yy in years:
        dir = "/data/exp/IceCube/%s/filtered/level2" % yy
        #dir = "/".join([dataPrefix, yy])
        list = os.listdir(dir)
        list.sort()
        dirp = re.compile("\d{4}")
        for folder in list:
            if dirp.match(folder):
                mm = folder[:2]
                dd = folder[2:]
                yymmdd.append([yy, mm, dd])
    yymmdd.sort()
    return yymmdd

# -----------------------------------------------------------------------------
def getRunFiles(ymdInfo, runList):
    """Get a list of sub-run files in a folder that are on the good run list
    """
    runFiles = []
    yy, mm, dd = ymdInfo
    runKey = ",".join(ymdInfo)

    if runKey in runList:
        runFiles = []
        for folder in runList[runKey]:
            files = os.listdir(folder)
            files.sort()
            filep = re.compile(r"Level2_IC86\.\d{4}_data_Run\d{8}_Subrun\d{8}_DST\.root")
            for file in files:
                if filep.match(file):
                    tokens = file.split("_")
                    runId = int(tokens[3][3:])
                    runFiles.append("".join([folder, file]))
    return runFiles

# -----------------------------------------------------------------------------
def getGoodRunList(filename):
    """Parse a list of Good and Bad runs store in the input file.
       Return a list of good runs with the full IC86 detector.  Note that the
       format of the file tends to change from year to year; this function will
       parse the good run lists produced automatically during Level 2
       processing.  For more information see

       https://wiki.icecube.wisc.edu/index.php/Goodrunlist
    """
    runList = defaultdict(list)

    # Search the file for good records
    file = open(filename, "r")
    for line in file:
        if line.startswith("RunNum"):
            continue

        # Tokenize the records
        tokens = line.split()
        runId, goodI3, goodIT, liveTime = [int(t) for t in tokens[0:4]]

        # Get the data directory name and split it apart
        dataDir = tokens[4]
        tokens = dataDir.split("/")
        yyyy = tokens[4]
        mm, dd = [tokens[7][:2], tokens[7][2:]]
        runKey = ",".join([yyyy,mm,dd])

        if goodI3==1:
            runList[runKey].append(dataDir)

    return runList

# -----------------------------------------------------------------------------
def submitDay(ymdInfo, runFileList, splitSize=20):
    """Submit the run files for processing, splitting up the batches if needed
    """
    # Split processing of all data in one 24 hour period into small chunks,
    # then merge the chunks...
    print("Creating DAG info for %s" % "-".join(ymdInfo))
    yy, mm, dd = ymdInfo
    jchunks = []
    outputs = []
    scripts = []
    # Process chunks for this yyyy-mm-dd
    for i, batchList in enumerate(splitSequence(runFileList, splitSize)):
        jobId = "simpleDST_%04s-%02s-%02s_p%02d" % (yy, mm, dd, i+1)
        if len(batchList) > 0:
            outf, script = processLevel2Files(jobId, batchList)
            jchunks.append(jobId)
            outputs.append(outf)
            scripts.append(script)
    # Merge chunks
    jobId = "simpleDST_%04s-%02s-%02s" % (yy, mm, dd)
    dest, script = mergeChunks(jobId, outputs)

    # Create a job hierarchy and return it so we can make a DAG file
    dag = []
    for j, s in zip(jchunks, scripts):
        dag.append("JOB %s %s" % (j, s))
    dag.append("JOB %s %s" % (jobId, script))
    for j in jchunks:
        dag.append("PARENT %s CHILD %s" % (j, jobId))

    return dag

# -----------------------------------------------------------------------------
def processLevel2Files(jobId, fileList):
    """Set up a condor job to process the DST files in the fileList.
    
        Returns the name of the output file created and the condor script used
        to create it, so that we can put it into a DAG.
    """
    # Create output file
    outd = "%s" % (outp)
    outf = "%s/ic86_%s.root" % (outd, jobId)
    dest = "%s/ic86_%s.root" % (prfx, jobId)

    # Build Simple-DST generator shell script
    fnam = ["%s/%s" % (outd, os.path.basename(f)) for f in fileList]
    argf = " \\\n".join(fnam)
    ocmd = " \\\n".join([envs, "dstcut", outf, argf])
    exeScript  = [
        "#!/bin/bash\n",
        "date",
        "hostname\n",
        "# Copy ROOT files to local disk for I/O",
        "mkdir -p %s" % outd,
        "cp %s %s" % (" \\\n".join(fileList), outd),
        "\n# Process files",
        ocmd,
        "\n# Remove copied Level2 ROOT files from local disk",
        "rm -f %s" % " \\\n".join(fnam),
        "\n# Move output file to final destination",
        "mv %s %s" % (outf, dest),
        "\ndate"
    ]
    scriptName = "%s/npx/exe/%s.sh" % (pdir, jobId)
    script = open(scriptName, "w")
    for line in exeScript:
        script.write("%s\n" % line)
    script.close()
    subprocess.call(["chmod", "755", scriptName])

    # Build a condor submission script to run the shell script
    cdrScript = [
        "Universe = vanilla",
        "getenv = true",
        "Executable = %s" % scriptName,
        "Log = %s/npx/log/%s.log" % (pdir, jobId),
        "Output = %s/npx/out/%s.out" % (pdir, jobId),
        "Error = %s/npx/err/%s.err" % (pdir, jobId),
        "Notification = NEVER",
        "Queue"
    ]
    scriptName = "%s/npx/sub/%s.sub" % (pdir, jobId)
    script = open(scriptName, "w")
    for line in cdrScript:
        script.write("%s\n" % line)
    script.close()

    return dest, scriptName

# -----------------------------------------------------------------------------
def mergeChunks(jobId, fileList):
    """Merge the Simple-DST chunks generated from the Level-2 ROOT files.
    """
    # Create output file
    outd = "%s" % (outp)
    outf = "%s/ic86_%s.root" % (outd, jobId)
    dest = "%s/ic86_%s.root" % (prfx, jobId)

    # Build Simple-DST generator shell script
    #fnam = ["%s/%s" % (outd, os.path.basename(f)) for f in fileList]
    fnam = fileList
    argf = " \\\n".join(fnam)
    ocmd = " \\\n".join([envs, "dstmerge", outf, argf])
    exeScript  = [
        "#!/bin/bash\n",
        "date",
        "hostname\n",
        #"# Copy ROOT files to local disk for I/O",
        "mkdir -p %s" % outd,
        #"cp %s %s" % (" \\\n".join(fileList), outd),
        "\n# Process files",
        ocmd,
        "\n# Move output file to final destination",
        "mv %s %s" % (outf, dest),
        "\n# Remove chunks from remote disk",
        "rm -f %s" % " \\\n".join(fileList),
        "\ndate"
    ]
    scriptName = "%s/npx/exe/%s.sh" % (pdir, jobId)
    script = open(scriptName, "w")
    for line in exeScript:
        script.write("%s\n" % line)
    script.close()
    subprocess.call(["chmod", "755", scriptName])

    # Build a condor submission script to run the shell script
    cdrScript = [
        "Universe = vanilla",
        "getenv = true",
        "Executable = %s" % scriptName,
        "Log = %s/npx/log/%s.log" % (pdir, jobId),
        "Output = %s/npx/out/%s.out" % (pdir, jobId),
        "Error = %s/npx/err/%s.err" % (pdir, jobId),
        "Notification = NEVER",
        "Queue"
    ]
    scriptName = "%s/npx/sub/%s.sub" % (pdir, jobId)
    script = open(scriptName, "w")
    for line in cdrScript:
        script.write("%s\n" % line)
    script.close()

    return dest, scriptName

# -----------------------------------------------------------------------------
# Set up the command line options
p = argparse.ArgumentParser(description="Simple DST ntuple maker")
p.add_argument("goodRunList", nargs=1,
               help="IceCube/IceTOP good run list ASCII file")
p.add_argument("-s", "--split", dest="split", type=int,
               default=10,
               help="Split sub-run processing into chunks of n files/chunk")
p.add_argument("-u", "--use-run", dest="YMD",
               default=None,
               help="Force use of year,month,day (comma-separated)")
args = p.parse_args()

runList = getGoodRunList(args.goodRunList[0])

# DAGMan job hierarchy
dag = []

# Mode 1: user has specified the run parameters to process (by date)
if args.YMD:
    ymd = args.YMD.split(",")
    runFiles = getRunFiles(ymd, runList)
    dag += submitDay(ymd, runFiles, args.split)
# Mode 2: all available good runs will be processed
else:
    for runKey in runList:
        ymd = runKey.split(",")
        runFiles = getRunFiles(ymd, runList)
        dag += submitDay(ymd, runFiles, args.split)

# Output DAG file
outdag = open("ntmaker.dag", "w")
for line in dag:
    outdag.write("%s\n" % line)
outdag.close()

