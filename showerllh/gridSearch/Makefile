SHELL    = /bin/sh
CXX		   = g++
CXXFLAGS = -g -Wall -O2 -DSDST_V3 -I. $(BOOST_FLAGS) $(PYTHON_FLAGS)
LD       = g++
LDFLAGS  = -O2 $(BOOST_LIBS) $(PYTHON_LIBS)

PYTHON_FLAGS = -I/usr/include/python2.6
PYTHON_LIBS  = -L/usr/lib/python2.6 -lpython

BOOST_FLAGS = -I$(I3_PORTS)/include/boost-1.38.0
BOOST_LIBS  = -L$(I3_PORTS)/lib/boost-1.38.0 -lboost_python

.PHONY: all clean

all: ctest

ctest: ctest.o
	$(LD) $^ $(LDFLAGS) -o $@

%.o:%.cc
	$(CXX) $(CXXFLAGS) -c $^

clean:
	rm -f *.o a.out ctest
