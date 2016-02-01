#!/bin/bash
date
hostname
cd /home/zgriffith/ShowerLLH/run_data
eval export SROOTBASE="/cvmfs/icecube.opensciencegrid.org/py2-v1" ;
export SROOT="/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64" ;
export I3_PORTS="/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64/i3ports" ;
export I3_SITE_CMAKE_DIR="/cvmfs/icecube.opensciencegrid.org/py2-v1/site_cmake" ;
export I3_DATA="/cvmfs/icecube.opensciencegrid.org/py2-v1/../data" ;
export PATH="/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64/i3ports/root-v5.34.18/bin:/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64/i3ports/qt-4.6.4/bin:/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64/bin:/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64/i3ports/bin:/home/zgriffith/.local/lib/python2.7/site-packages:/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64/i3ports/root-v5.34.18/bin:/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64/i3ports/qt-4.6.4/bin:/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64/bin:/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64/i3ports/bin:/usr/local/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin" ;
export MANPATH="/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64/man:/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64/share/man:/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64/man:/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64/share/man:" ;
export PKG_CONFIG_PATH="/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64/lib/pkgconfig:/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64/lib/pkgconfig:" ;
export LD_LIBRARY_PATH="/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64/i3ports/root-v5.34.18/lib:/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64/i3ports/qt-4.6.4/lib:/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64/i3ports/lib/Minuit2-5.24.00:/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64/i3ports/lib/boost-1.38.0:/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64/i3ports/lib/log4cplus-1.0.4:/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64/lib:/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64/i3ports/lib:/cvmfs/icecube.opensciencegrid.org/py2-v1/../distrib/OpenCL_RHEL_6_x86_64/lib/RHEL_6_x86_64:/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64/i3ports/root-v5.34.18/lib:/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64/i3ports/qt-4.6.4/lib:/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64/i3ports/lib/Minuit2-5.24.00:/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64/i3ports/lib/boost-1.38.0:/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64/i3ports/lib/log4cplus-1.0.4:/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64/lib:/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64/i3ports/lib::/cvmfs/icecube.opensciencegrid.org/py2-v1/../distrib/jdk1.6.0_24_RHEL_6_x86_64/lib:/cvmfs/icecube.opensciencegrid.org/py2-v1/../distrib/jdk1.6.0_24_RHEL_6_x86_64/jre/lib:/cvmfs/icecube.opensciencegrid.org/py2-v1/../distrib/jdk1.6.0_24_RHEL_6_x86_64/jre/lib/amd64:/cvmfs/icecube.opensciencegrid.org/py2-v1/../distrib/jdk1.6.0_24_RHEL_6_x86_64/jre/lib/amd64/server:/cvmfs/icecube.opensciencegrid.org/py2-v1/../distrib/jdk1.6.0_24_RHEL_6_x86_64/jre/lib/i386:/cvmfs/icecube.opensciencegrid.org/py2-v1/../distrib/jdk1.6.0_24_RHEL_6_x86_64/jre/lib/i386/server:/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64/tools/gfortran:/cvmfs/icecube.opensciencegrid.org/py2-v1/../distrib/jdk1.6.0_24_RHEL_6_x86_64/lib:/cvmfs/icecube.opensciencegrid.org/py2-v1/../distrib/jdk1.6.0_24_RHEL_6_x86_64/jre/lib:/cvmfs/icecube.opensciencegrid.org/py2-v1/../distrib/jdk1.6.0_24_RHEL_6_x86_64/jre/lib/amd64:/cvmfs/icecube.opensciencegrid.org/py2-v1/../distrib/jdk1.6.0_24_RHEL_6_x86_64/jre/lib/amd64/server:/cvmfs/icecube.opensciencegrid.org/py2-v1/../distrib/jdk1.6.0_24_RHEL_6_x86_64/jre/lib/i386:/cvmfs/icecube.opensciencegrid.org/py2-v1/../distrib/jdk1.6.0_24_RHEL_6_x86_64/jre/lib/i386/server" ;
export PYTHONPATH="/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64/i3ports/root-v5.34.18/lib:/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64/lib/python2.7/site-packages:/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64/i3ports/lib/python2.7/site-packages:/home/zgriffith/.local/lib/python2.7/site-packages:/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64/i3ports/root-v5.34.18/lib:/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64/lib/python2.7/site-packages:/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64/i3ports/lib/python2.7/site-packages::/home/zgriffith/:/home/zgriffith/Modules:/home/zgriffith/ShowerLLH/segments:/home/zgriffith/showerllh" ;
export ROOTSYS="/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64/i3ports/root-v5.34.18" ;
export OS_ARCH="RHEL_6_x86_64" ;
export GCC_VERSION="4.4.7" ;
export JAVA_HOME="/cvmfs/icecube.opensciencegrid.org/py2-v1/../distrib/jdk1.6.0_24_RHEL_6_x86_64" ;
export GOTO_NUM_THREADS="1" ;
export PERL5LIB="/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64/lib/perl:/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64/lib/perl5:/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64/lib/perl5/site_perl:/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64/lib/perl:/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64/lib/perl5:/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64/lib/perl5/site_perl:" ;
export GLOBUS_LOCATION="/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64" ;
export PATH="/home/zgriffith/.local/lib/python2.7/site-packages:/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64/i3ports/root-v5.34.18/bin:/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64/i3ports/qt-4.6.4/bin:/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64/bin:/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64/i3ports/bin:/usr/local/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin" ;
export PYTHONPATH="/home/zgriffith/.local/lib/python2.7/site-packages:/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64/i3ports/root-v5.34.18/lib:/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64/lib/python2.7/site-packages:/cvmfs/icecube.opensciencegrid.org/py2-v1/RHEL_6_x86_64/i3ports/lib/python2.7/site-packages::/home/zgriffith/:/home/zgriffith/Modules:/home/zgriffith/ShowerLLH/segments:/home/zgriffith/showerllh" ;
/data/user/zgriffith/offline/build/env-shell.sh python merger.py IT81-II 201306
date
