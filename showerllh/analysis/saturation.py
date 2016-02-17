#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

import myGlobals as my
from sim import eres


def load_sim():

    my.setupShowerLLH(verbose=False)
    s, files = {},{}
    files['old'] = '%s/original/IT73_sim/SimPlot_standard.npy' % my.llh_data
    files['new'] = '%s/IT73_sim/SimPlot_standard.npy' % my.llh_data
    for key in files:
        q = np.load(files[key])
        s[key] = q.item()

    return s


def my_eres(s):

    eres(s)


if __name__ == "__main__":

    s = load_sim()
    my_eres(s)
