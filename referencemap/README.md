
NOTE: IceTray environment must be loaded before compiling TimeScrable.cpp

trunk verison of healpix-cxx from:
http://code.icecube.wisc.edu/svn/projects/healpix-cxx/trunk/

astro from:
http://code.icecube.wisc.edu/projects/icecube/browser/IceCube/sandbox/cfinley/astro/trunk



For some reason libphotospline.so could not be found when linking TimeScrable.cpp.
To fix this, the following line needed to be added to my .zshrc (or .bashrc if using a bash shell):

`export LD_LIBRARY_PATH="/data/user/jbourbeau/offline-software-2/V04-08-00/build/lib:$LD_LIBRARY_PATH"`
