#!/usr/bin/env python
from timeScramble import timeScramble
import sys
data_map = sys.argv[1]
scramble_time = 24 #in hours
timeScramble(data_map, scramble_time, out=True)
