#!/usr/bin/python

import sys, time
import numpy as np
from math import *
import argparse
import general_scripts as gs

from distutils.version import LooseVersion

def read_file(fn, field, key="none"):
    legs=[]
    nplots=0
    nlines=[]
    xlist=[]
    ylist=[]
    tmp=0
    tmpx=[]
    tmpy=[]
    if key=="none":
        dumLeg=1
        with open(fn) as fp:
            for line in fp:
                l = line.split()
                if line[0] == '#' or line[0] == '@' or line == '\n':
                    continue
                if l[0].isalnum():
                    continue
                if line[0].isalpha():
                    # Catch Diamond q(Angs) header.
                    continue
                if l[0].endswith(';') or l[0].endswith(':'):
                    # Hamburg and Australian Synchrotron specific, to catch
                    # opening arguements that are not technically alphanumeric.
                    continue
                elif (line[0] == '&'):
                    nlines.append(tmp)
                    xlist.append(tmpx)
                    ylist.append(tmpy)
                    legs.append(str(dumLeg))
                    dumLeg+=1
                    tmpx=[]
                    tmpy=[]
                    tmp=0
                    nplots+=1
                else:
                    tmpx.append(l[0])
                    tmpy.append(l[field])
                    tmp += 1
        if len(tmpx) > 0:
            nplots+=1
            legs.append(str(dumLeg))
            nlines.append(tmp)
            xlist.append(tmpx)
            ylist.append(tmpy)
    else:
        with open(fn) as fp:
            for line in fp:
                l = line.split()
                if (line[0] == '#' or line[0] == '@' or line == '\n'):
                    if (line[0] == '@' and key in line):
                        legs.append( l[-1].strip('"') )
                        nplots += 1
                    continue
                elif (line[0] == '&'):
                    nlines.append(tmp)
                    xlist.append(tmpx)
                    ylist.append(tmpy)
                    tmpx=[]
                    tmpy=[]
                    tmp=0
                else:
                    tmpx.append(l[0])
                    tmpy.append(l[field])
                    tmp += 1
        if len(tmpx) > 0:
            nlines.append(tmp)
            xlist.append(tmpx)
            ylist.append(tmpy)

    #Sanity check
    if len(legs) != len(nlines):
        print( "= = ERROR: number of legends(%d) is not equal to the number of plots (%d)!" % (len(legs), len(nlines)), file=sys.stderr )
        sys.exit(1)
    #for i in range(len(nlines)):
    #    for j in range(nlines[i]):

    return legs, xlist, ylist


if __name__ == '__main__':
#Parsing.
    parser = argparse.ArgumentParser(description='Takes a number of xmgrace-like files containing equivalent data'
                                                 'each from a different replicate, and perform averaging across'
                                                 'these files.'
                                                 'More info: For each file containing sets of difference curves,'
                                                 'e.g. s0, s1, s2... perform average across different files while'
                                                 'preserving the set layout.',
                                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('filelist', metavar='file', type=str, nargs='+',
                        help='List of file names of the data files to be averaged.')
    parser.add_argument('-o', '--outfile', type=str, dest='out_file', default='out',
                        help='Output file for the averaged data.')
    parser.add_argument('-s', '--search_key', type=str, dest='key', default='legend',
                        help='String to search to identify keys that delineate different graphs.'
                        'For example, "legend" will isolate xmgrace command: @s# legend "txt".'
                        'Put in "none" to have this not use keys.' )
    parser.add_argument('-f', '--field', type=str, dest='field', default='end',
                        help='Which field to use as a key. Defaults to the last field. This is 0-indexed like python standard.')

    time_start = time.time()
    args = parser.parse_args()

    out_filename=args.out_file
    nfiles = len(args.filelist)
    if nfiles < 2:
        print( "= = ERROR: this script averages data from multiple curves!", file=sys.stderr )
        sys.exit(-1)
    if args.field == "end":
        field=-1
    else:
        field=int(args.field)

    bFirst = True
    for i in range(nfiles):
        fileName = args.filelist[i]
        if fileName[0] == '#':
            print( "= = NOTE: skipping file argument %s" % fileName, file=sys.stderr )
            continue
        l, x, y = read_file(fileName, field, args.key)
        x=np.array(x,dtype=np.float64)
        y=np.array(y,dtype=np.float64)
        print( " ...plot %s read." % fileName, file=sys.stderr )
        nplot = len(l)
        ndat  = len(x[0])
        if bFirst:
            bFirst = False
            leglist=[]
            xlist=np.zeros((nfiles,nplot,ndat),dtype=np.float64)
            ylist=np.zeros((nfiles,nplot,ndat),dtype=np.float64)
            check=(nplot,ndat)
            leglist.append(l)
            xlist[i]=x
            ylist[i]=y
            continue
        #Sanity check
        if check[0] != nplot:
            print( "= = ERROR: Input data files do not contain the same number of plots!", file=sys.stderr )
            sys.exit(2)
        # Check if X-values are identical!
        # print( xlist[0].shape, x.shape )
        # print( np.array_equal( xlist[0], x) )
        if not np.array_equal( xlist[0], x):
            # Try to interpolate instead.
            print( "= = WARNING: The latest input data file %s does not contain identical X-values!" % fileName, file=sys.stderr )
            print( "= = ...will use interpolation.", file=sys.stderr )
            y2=np.zeros(check)
            for j in range(nplot):
                #print( "Debug:", type(xlist[i,j]), type(x[j]) )
                yTemp = np.interp(xlist[0,j,:],x[j,:],y[j,:])
            for k in range(ndat):
                y2[0,k] = yTemp[k] 
            leglist.append(l)
            xlist[i]=xlist[0]
            ylist[i]=y2
        else:
            leglist.append(l)
            xlist[i]=x
            ylist[i]=y
    print( " ...all plots read. Conducting averaging.", file=sys.stderr )
    yavg=ylist.mean(axis=0)
    ystd=ylist.std(axis=0)/sqrt(nfiles-1)
    print( " ...average finished.", file=sys.stderr )

    if LooseVersion(np.version.version) >= LooseVersion('1.10'):
        gs.print_sxylist(out_filename, leglist[0], x[0], np.stack((yavg,ystd), axis=-1) )
    else:
        shape=list(yavg.shape)
        shape.append(2)
        tmp = np.zeros( shape, dtype=yavg.dtype )
        tmp[...,0] = yavg
        tmp[...,1] = ystd
        gs.print_sxylist(out_filename, leglist[0], x[0], tmp )
