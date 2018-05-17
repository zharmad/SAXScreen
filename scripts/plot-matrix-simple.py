#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
import general_scripts as gs
import argparse

parser = argparse.ArgumentParser(description="Simple tool to plot matrix data.")
parser.add_argument('-f', dest='inFile', type=str, default="./fitted-CorMapMatrix-1.dat",
                    help="Name of the matrix file to be plotted.")
parser.add_argument('-o', dest='outFile', type=str, default="",
                    help="Optional output file to print. Otherwise defaults to stdout.")

args = parser.parse_args()

inFile=args.inFile
X = gs.load_matrix(inFile)

fig, ax = plt.subplots()
ax.xaxis.tick_top()
plt.imshow(X, cmap='Greys',  interpolation='nearest')

if args.outFile != '':
    plt.savefig(args.outFile)
else:
    plt.show()
