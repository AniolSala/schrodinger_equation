#!/usr/bin/env python3

# Script that plots the real part of the the solution in one
# subplot and the error in a second subplot.
#
# The error is estimated computing (for each time step)
# the area of the squared module of the solution (using the simps
# method of scipy.integrate).


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
from getData import DataParser


def main():
    filename = "schrEq1D_gaussian_wpacketpot.txt"
    dataP = DataParser(filename)
    nt, nx, dt, dx = dataP.getParams()
    plots = dataP.getPlotsList()

    indList = np.arange(0, nt, nt // 9, dtype=int)

    fig, [axsol, axerr] = plt.subplots(figsize=(15, 5), ncols=2)

    # Solution plot
    for ind in indList:
        plot = plots[ind]
        x, phiReal2 = plot[1, :], plot[4, :]
        axsol.plot(x, phiReal2, label="t={}".format(plot[0, 0]))
    axsol.legend()
    axsol.set_xlabel("x")
    axsol.set_ylabel("real(ψ)")

    # Error plot
    freq = 1
    iters = range(len(plots))
    area = np.array(
        [simps(pl[4, ::freq], pl[1, ::freq]) for pl in plots])
    error = 1. - area

    axerr.plot(iters, error)
    axerr.set_xlabel("Iteration")
    axerr.set_ylabel("error")

    plt.show()


if __name__ == '__main__':
    main()
