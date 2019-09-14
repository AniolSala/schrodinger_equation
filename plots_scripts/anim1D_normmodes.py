#!/usr/bin/env python3

import os
import seaborn as sns
import matplotlib.pyplot as plt
from getData import DataParser
from matplotlib.animation import FuncAnimation
plt.rcParams.update({'font.size': 18})


# Plot the first four modeReals
nPlots = 4
filenames = ["schrEq1D_normmode_{}_box.txt".format(i + 1)
             for i in range(nPlots)]
outputName = "schrEq1D_normmodes_box.gif"
dataPars = [DataParser(filename) for filename in filenames]
nFrames = [dataP.getPlotsList() for dataP in dataPars]
# nFrames[modeReal][frame] = (nrows, ncols) ndarray

# Define variables for the animation
sns.set()
sns.set_context("paper")
sns.set(font_scale=2)
fig, axs = plt.subplots(2, 2, figsize=(15, 15), sharex=True, sharey=True)

# The data to plot
# Real part
nmodeReal1, = axs[0, 0].plot(
    nFrames[0][0][1, :], nFrames[0][0][2, :], "r-", label=r'real($\psi$)')
nmodeReal2, = axs[0, 1].plot(
    nFrames[1][0][1, :], nFrames[1][0][2, :], "r-", label=r'real($\psi$)')
nmodeReal3, = axs[1, 0].plot(
    nFrames[2][0][1, :], nFrames[2][0][2, :], "r-", label=r'real($\psi$)')
nmodeReal4, = axs[1, 1].plot(
    nFrames[3][0][1, :], nFrames[3][0][2, :], "r-", label=r'real($\psi$)')

# Imaginary part
nmodeImag1, = axs[0, 0].plot(
    nFrames[0][0][1, :], nFrames[0][0][3, :], "r-", color='g',
    label=r'imag($\psi$)')
nmodeImag2, = axs[0, 1].plot(
    nFrames[1][0][1, :], nFrames[1][0][3, :], "r-", color='g',
    label=r'imag($\psi$)')
nmodeImag3, = axs[1, 0].plot(
    nFrames[2][0][1, :], nFrames[2][0][3, :], "r-", color='g',
    label=r'imag($\psi$)')
nmodeImag4, = axs[1, 1].plot(
    nFrames[3][0][1, :], nFrames[3][0][3, :], "r-", color='g',
    label=r'imag($\psi$)')

axs[0, 0].legend(loc='upper right')
axs[0, 1].legend(loc='upper right')
axs[1, 0].legend(loc='upper right')
axs[1, 1].legend(loc='upper right')

# x-y labels
for i, ax in enumerate(axs.flat):
    ax.set(xlabel=r'$x (eV^{-1})$',
           title=r'modeReal $n_x$ = {}'.format(i + 1),
           xlim=(-.6, .6))
    ax.label_outer()

print("All plot settings done")


def update(frame):
    # Real part
    nmodeReal1.set_ydata(nFrames[0][frame][2, :])
    nmodeReal2.set_ydata(nFrames[1][frame][2, :])
    nmodeReal3.set_ydata(nFrames[2][frame][2, :])
    nmodeReal4.set_ydata(nFrames[3][frame][2, :])

    # Imaginary part
    nmodeImag1.set_ydata(nFrames[0][frame][3, :])
    nmodeImag2.set_ydata(nFrames[1][frame][3, :])
    nmodeImag3.set_ydata(nFrames[2][frame][3, :])
    nmodeImag4.set_ydata(nFrames[3][frame][3, :])
    return nmodeReal1, nmodeReal2, nmodeReal3, nmodeReal4, nmodeImag1, nmodeImag2, nmodeImag3, nmodeImag4,


output = os.path.join(dataPars[0].getDataDir(), "plots", "normmodes")
outputName = os.path.join(output, outputName)
animation = FuncAnimation(
    fig, update, frames=range(1, len(nFrames[0]) // 4, 1), interval=400, blit=True)
animation.save(outputName, writer="imagemagick", fps=50, dpi=40)
