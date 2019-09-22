#!/usr/bin/env python3

#
# Script to make an animated subplot of the solution and its area
#
# It is a mess, but main idea is to define two subplots. In the
# first subplot we will plot the real and complex part of the
# solution.
# In the second subplot we will plot the squared module of the
# solution. If the solution is a mixed state (built from various
# states), then in this subplot show also the squared module of
# these states. These states are got in the `modes` list.
#


import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from getData import DataParser
from matplotlib.animation import FuncAnimation
from scipy.integrate import simps


def getError(plotsList):
    dx = 10 * plotsList[0][0, 1] - plotsList[0][0, 0]
    areas = [simps(pl[4, :], pl[1, :], dx) for pl in plotsList]
    return np.array(areas)


# Files to plot
filename = "schrEq1D_gaussian_wpacketpot_high.txt"
outputName = os.path.splitext(filename)[0] + "_area.gif"
modes = [char for char in filename.split(
    "_") if char.isdigit()]  # Take the modes we want to plot from the filename

# Extract the data from the txt file
dataPars = DataParser(filename)
data = dataPars.getData()
frames = dataPars.getPlotsList()

# Define variables for the animation
sns.set()
sns.set_context("paper")
sns.set(font_scale=2)
fig, plots = plt.subplots(
    2,
    figsize=(20, 15),
    sharex=True,
    gridspec_kw={'hspace': 0}
)

axs, area = plots[0], plots[1]

# Solution plot (real part)
phiReal, = axs.plot(
    frames[0][1, :],        # x data
    frames[0][2, :],        # y data
    "r-",                   # style
    label=r'real($\psi$)'   # legend
)

# Solution plot (imaginary part)
phiImag, = axs.plot(
    frames[0][1, :],        # x data
    frames[0][3, :],        # y data
    "r-",                   # style
    color='g',              # color
    label=r'imag($\psi$)'   # legend
)

# Set axis properties
x0 = round(max(data[1, :]) * 1.1, 2)
axs.set(
    xlabel=r'$x (eV^{-1})$',
    ylabel=r'$\psi$',
    xlim=(-x0, x0),
    ylim=(min(data[2, :]) - .1, max(data[2, :]) + .1)
)


# Error / area plot
yErrLabel = r'$\int\vert\psi\vert^2\,dx$'
xErrLabel = r'time $eV^{-1}$'
yAreaLabel = r'$\vert\psi\vert^2$'
xAreaLabel = r'x $(eV^{-1})$'


pl2, = area.plot(frames[0][1, :],
                 frames[0][4, :],
                 "r-",
                 rasterized=True)

# Plot also psi² for the desired modes
for i in modes:
    modePhi2 = DataParser(
        "schrEq1D_normmode_{}_box.txt".format(i)).getPlotsList()[0]
    area.plot(modePhi2[1, :],
              modePhi2[4, :],
              '--',
              label=r'$\vert \psi \vert^2$ for mode n = ' + '{}'.format(i))

# Plot also the potential barrier
xb, d, heigh = 0., .05, max(data[4, :] * 1.05)  # 1.05 / 0.3
xBarrier = [xb, xb, xb + d, xb + d]
yBarrier = [0, heigh, heigh, 0]
area.plot(xBarrier, yBarrier, "--", color="black")

for plot in plots:
    plot.label_outer()

x0 = .5  # round(max(data[1, :]) * 1.1, 1)
y0 = round(max(data[4, :]) * 1.1, 1)
area.set(
    xlabel=xAreaLabel,
    ylabel=yAreaLabel,
    xlim=(-x0, x0),
    ylim=(-.1, y0)
)

axs.legend(loc='upper right')
area.legend(loc='upper right')

print("Plot settings done")


def update(frame):
    phiReal.set_ydata(frames[frame][2, :])
    phiImag.set_ydata(frames[frame][3, :])

    pl2.set_ydata(frames[frame][4, :])

    return phiReal, phiImag, pl2,


nFrames = int(len(frames) * .6)
output = os.path.join(dataPars.getDataDir(), "plots", "areas")
outputName = os.path.join(output, outputName)
animation = FuncAnimation(fig, update, frames=range(1, nFrames, 5), blit=True)
animation.save(outputName, writer="imagemagick", fps=40, dpi=25, metadata=None)
