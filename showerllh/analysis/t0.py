#!/usr/bin/env python

import matplotlib.pyplot as plt

def figtest(a, b, ax=None):
    if ax != None:
        ax.plot(a, b, '.')

if __name__ == "__main__":

    fig, ax = plt.subplots()
    a = range(5)
    b = range(5)
    figtest(a, b, ax=ax)
    a = range(5)
    b = range(1,6)
    figtest(a, b, ax=ax)
    plt.show()
