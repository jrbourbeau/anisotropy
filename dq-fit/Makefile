SHELL    = /bin/sh
CXX      = g++
#CXXFLAGS = $(HEALPIXCPPFLAGS) $(HEALPIXCXXFLAGS) $(ROOT_FLAGS) -g -O2 -I. 
CXXFLAGS = -g -O2 -I. $(HPFLAGS) $(BOOST_FLAGS) $(ROOT_FLAGS)
LD       = g++
#LDFLAGS  = $(HEALPIX_LIBS) $(ROOT_LIBS) -lMinuit 
LDFLAGS  = $(HPLIBS) $(BOOST_LIBS) $(ROOT_LIBS) -lMinuit 

#HEALPIXCPPFLAGS = $(shell healpix-cxx-config --cppflags)
#HEALPIXCXXFLAGS = $(shell healpix-cxx-config --cxxflags)
#HEALPIX_LIBS    = $(shell healpix-cxx-config --ldflags)

HPFLAGS = -I$(I3_SRC)/healpix-cxx/public -DHEALPIX_STANDALONE
HPLIBS  = -L$(I3_BUILD)/lib -lfftpack -lcxxsupport -lhealpix-cxx

BOOST_FLAGS = -I$(I3_PORTS)/include/boost-1.38.0
BOOST_LIBS  = -L$(I3_PORTS)/lib/boost-1.38.0

ROOT_FLAGS      = $(shell root-config --cflags)
ROOT_LIBS       = $(shell root-config --libs)

.PHONY: all clean

all: dq-fit

dq-fit: dq-fit.cc
	$(CXX) $(LDFLAGS) $(CXXFLAGS) $^ -o $@

clean:
	rm -f dq-fit
