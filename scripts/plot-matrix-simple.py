#!/usr/bin/python

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import numpy as np
import general_scripts as gs
import argparse

parser = argparse.ArgumentParser(description="Simple tool to plot matrix data.")
parser.add_argument('-f', dest='inFile', type=str, default="./fitted-CorMapMatrix-1.dat",
                    help="Name of the matrix file to be plotted.")
parser.add_argument('--period', type=float, default=None,
                    help="Period between major ticks")
#parser.add_argument('--cmin', type=float, default=None,
#                    help="Specified minimum of the colorbar.")
parser.add_argument('--cmax', type=float, default=None,
                    help="Specified maximum of the colorbar.")
parser.add_argument('--cmap', type=str, default="Greys",
                    help="Color-map name. Other examples like viridis, plasma, and hot will also work."
                         "See, e.g.: https://matplotlib.org/examples/color/colormaps_reference.html")
parser.add_argument('-o', dest='outFile', type=str, default="",
                    help="Optional output file to print. Otherwise defaults to stdout.")

args = parser.parse_args()

inFile=args.inFile
X = gs.load_matrix(inFile)

fig, ax = plt.subplots()
ax.xaxis.tick_top()

minorLocator = MultipleLocator(1)
ax.xaxis.set_minor_locator(minorLocator)
ax.yaxis.set_minor_locator(minorLocator)
if not args.period is None:
    plt.xticks(np.arange(0,X.shape[0],args.period))
    ax.xaxis.set_major_formatter(FormatStrFormatter('Titr. %d'))
    plt.yticks(np.arange(0,X.shape[1],args.period))
    ax.yaxis.set_major_formatter(FormatStrFormatter('Titr. %d'))

if not args.cmax is None:
    p = plt.imshow(X, cmap=args.cmap, vmax=args.cmax, interpolation='nearest')
else:
    p = plt.imshow(X, cmap=args.cmap, interpolation='nearest')

fig.colorbar(p)

if args.outFile != '':
    plt.savefig(args.outFile)
else:
    plt.show()
