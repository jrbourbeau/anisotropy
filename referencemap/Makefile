SHELL         = /bin/sh
CXX           = g++
CXXFLAGS      = -g -Wall -O2 -DSDST_V3 -I. $(BOOST_FLAGS) $(ROOT_FLAGS) $(ICETRAY_FLAGS) $(SPLINE_FLAGS) $(HEALPIX_FLAGS)
LD            = g++
LDFLAGS       = -O2 $(BOOST_LIBS) $(ROOT_LIBS) $(ICETRAY_LIBS) $(SPLINE_LIBS) $(HEALPIX_LIBS)

ROOT_FLAGS    = $(shell root-config --cflags)
ROOT_LIBS     = $(shell root-config --libs)

ICETRAY_FLAGS = $(shell icetray-config astro healpix-cxx --cflags)
ICETRAY_LIBS  = $(shell icetray-config astro healpix-cxx --ldflags)
#ICETRAY_FLAGS = $(shell icetray-config astro --cflags)
#ICETRAY_LIBS  = $(shell icetray-config astro --ldflags)

#ICEREC = /data/user/fmcnally/offline/V04-08-00
#SPLINE_LIBS    = $(ICEREC)/build/lib/libphotospline.so
#SPLINE_FLAGS   = -I$(ICEREC)/src/photospline/public

ICEREC = /data/user/jbourbeau/offline/V04-08-00
# ICEREC = /data/user/jbourbeau/offline-software-2/V04-08-00
SPLINE_LIBS    = -L$(ICEREC)/build/lib -lphotospline
SPLINE_FLAGS   = -I$(ICEREC)/src/photospline/public

BOOST_FLAGS = -I$(I3_PORTS)/include/boost-1.38.0
BOOST_LIBS = -L$(I3_PORTS)/lib/boost-1.38.0 -lboost_program_options-mt-d -lboost_filesystem-mt-d

# HEALPIX_FLAGS = -I$(I3_PORTS)/include/healpix_cxx
# HEALPIX_LIBS = -L$(I3_PORTS)/lib -lhealpix_cxx
# Attempt to use updated boost version conflicts with cvmfs/icetray version
#BOOST_FLAGS = -I/data/user/fmcnally/boost_1_57_0/build/include
#BOOST_LIBS  = -L/data/user/fmcnally/boost_1_57_0/build/lib -lboost_program_options

.PHONY: all clean

all: timeScramble

timeScramble: TimeScramble.o SimpleDST.o
	$(LD) $^ $(LDFLAGS) -o $@

%.o:%.cpp
	$(CXX) $(CXXFLAGS) -c $^

clean:
	rm -f *.o a.out timeScramble
