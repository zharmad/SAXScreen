#!/usr/bin/python
import sys, math, itertools
import numpy as np
from scipy.optimize import fmin_powell
import random
import argparse
from ast import literal_eval as make_tuple

def read_xvg(fn):
    fp = open(fn, 'r')
    x  = []
    y  = []
    dy = []
    for line in fp.readlines():
        if (not (line[0] == '#' or line[0] == '@' or line[0] == '&' or line == '\n')):
            l = line.split()
            x.append(float(l[0]))
            y.append(float(l[1]))
            if (len(l) > 2):
                dy.append(float(l[2]))
    fp.close()
    return np.array( [x, y, dy] )
#    return np.stack( (x, y, dy), axis=0 )

def read_xvgs_all(filelist):
    out_list=[]
    for i in range(len(filelist)):
        block = read_xvg( filelist[i] )
        print >> sys.stderr, "  ...file read, shape: %s" % str(block.shape)
        out_list.append(block)
    return out_list

def write_to_xvg(fp, x, y, dy):
    if (len(dy) > 0):
        print >> fp, '@type xydy'
        for i in range(len(x)):
            print >> fp, '%g %g %g' % (x[i], y[i], dy[i])
    else:
        print >> fp, '@type xy'
        for i in range(len(x)):
            print >> fp, '%g %g' % (x[i], y[i])
    print >> fp, '&'

def write_xvg_raw(fn, db, header='', bReverse=False, bAmpersand=True):
    fp = open(fn, 'w')
    if header != '':
        print >> fp, header
    if isinstance( db, list):
        for i in len(db):
            sh = db[i].shape
            if len(sh)==2:
                if bReverse:
                    for j in range(sh[-1]):
                        for k in range(sh[0]):
                            print >> fp, '%f ' % db[i][k,j],
                        print >> fp, ''
                else:
                    for j in range(sh[0]):
                        for k in range(sh[-1]):
                            print >> fp, '%f ' % db[i][j,k],
                        print >> fp, ''
            if bAmpersand:
                print >> fp, '&'
    elif isinstance( db, np.ndarray):
        sh = db.shape
        if len(sh)==2:
            if bReverse:
                for j in range(sh[-1]):
                    for k in range(sh[0]):
                        print >> fp, '%f ' % db[k,j],
                    print >> fp, ''
            else:
                for j in range(sh[0]):
                    for k in range(sh[-1]):
                        print >> fp, '%f ' % db[j,k],
                    print >> fp, ''
        if bAmpersand:
            print >> fp, '&'

    fp.close()

def write_xvg(fn, x, dataBlock, header='', bType=False, bTarget=False, graph=0):
    fp = open(fn, 'w')
    sh = dataBlock.shape
    nVals = len(x)
    if header != '':
        print >> fp, header
    # 3D - multiple plots sharing the same X
    if len(sh) == 3:
        if bType:
            if sh[1] == 2:
                print >> fp, '@type xydy'
            elif sh[1] == 1:
                print >> fp, '@type xy'
        series=0
        # Each axis-0 in dataBlock represents a graph that must be plotted.
        for i in range(sh[0]):
            if bTarget:
                print >> fp, '@target g%i.s%i' % (graph, series)
            for h in range(nVals):
                print >> fp, '%g ' % x[h]
                for j in range(sh[1]):
                    print >> fp, '%g ' % dataBlock[i,j,h],
                print >> fp, ''
            print >> fp, '&'
            series+=1
    # 2D - a single plot with a separate X.
    elif len(sh) == 2:
        if bType:
            if sh[0] == 2:
                print >> fp, '@type xydy'
            elif sh[0] == 1:
                print >> fp, '@type xy'
        if bTarget:
            print >> fp, '@target g%i.s%i' % (graph, series)
        for h in range(nVals):
            print >> fp, '%g ' % x[h],
            for j in range(sh[0]):
                print >> fp, '%g ' % dataBlock[j,h],
            print >> fp, ''
        print >> fp, '&'
    # 1D - single plot again with separate X.
    elif len(sh) == 1:
        if bType:
            print >> fp, '@type xy'
        for h in range(nVals):
            print >> fp, '%g %g' % (x[h], dataBlock[h])
        print >> fp, '&'

    fp.close()

def find_subrange(x, xmin, xmax):
    return x[ (x>=xmin) & (x<=xmax) ]

def find_subrange_all( dataRaw, initBounds):
    nPlots = len(dataRaw)
    minAll = np.max( (initBounds[0], np.max([ d[0,0] for d in dataRaw ])) )
    maxAll = np.min( (initBounds[1], np.min([ d[0,-1] for d in dataRaw ])) )

    return find_subrange( dataRaw[0][0], minAll, maxAll )

def build_interpolated_block( xBasis, dataRaw ):
    nVals  = len(xBasis)
    nPlots = len(dataRaw)
    # Assume Each file has the same numver of y-values.
    nCols  = dataRaw[0].shape[0]
    outBlock = np.zeros( (nPlots, nCols-1, nVals) )

    for i in range(nPlots):
        for j in range(1, nCols):
            outBlock[i,j-1] = np.interp( xBasis, dataRaw[i][0], dataRaw[i][j])

    return outBlock

def subtract_spectra(yA, yB):
    # A-B
    sh  = yA.shape
    out = np.zeros( sh )
    out[0] = yA[0] - yB[0]
    out[1] = np.sqrt( np.square(yA[1]) + np.square(yB[1]) )
    return out

def intensityDiff(pos, *args):
    y1    = args[0][0]
    y1sig = args[0][1]
    y2    = args[1][0]
    y2sig = args[1][1]
    bLog  = args[2]
    bUseWeights = args[3]
    bNoConst    = args[4]

    #print y1[0], y1sig[0], y2[0], y2sig[0]
    # Final sanity check..
    if (len(y1) != len(y2)):
        printf >> sys.stderr, 'Different lengths of y1 and y2'
        sys.exit(1)
    if bNoConst:
        f = pos[0] ; c=0.0
    else:
        f, c = pos

#    print 'Inside IntensityDiff: f / c = %g / %g - bLog %r - bUseWeights %r ' % (f, c, bLog, bUseWeights)
    sumChi2 = 0.
    sumSig2 = 0.
    for i in range(len(y1)):
        thisy2    = f*y2[i] + c
        if len(y2sig) > 0:
            thisy2sig = f*y2sig[i]
        # doing the difference between the logarithms of the intensities
        if (thisy2 < 0):
            # don't allow negative I(q)
            return 1e20
        if bLog:
            tmp  = (math.log(y1[i])-math.log(thisy2))**2
            if bUseWeights:
                sig2 = (thisy2sig/thisy2)**2 + (y1sig[i]/y1[i])**2
            else:
                sig2 = 1.
        else:
            tmp  = (y1[i]-thisy2)**2
            if bUseWeights:
                sig2 = thisy2sig**2 + y1sig[i]**2
            else:
                sig2 = 1.
        sumChi2 += tmp/sig2
        sumSig2 += sig2
        #print '%d  / %12g %12g - chi2 %12g' % (i, tmp, sig2, sumChi2)
#    print 'In function IntensityDiff, returning %g' % (sumChi2/len(y1))
#    meanSig2 = sumSig2/len(y1)
#    return sumChi2/len(y1)*meanSig2
    return sumChi2/len(y1)

def intChiFree(x1, y1, y2, y2sig, f, c, bLog, dmax, nk, nparams):
    if (len(y1) != len(y2)):
        printf >> sys.stderr, 'Different lengths of y1 and y2'
        sys.exit(1)
    #Define channon channel width
    swidth = math.pi/dmax
    #Define number of q-points
    nq   = len(y1)
    qmin = x1[0]
    qmax = x1[nq-1]
    nbins = math.floor((qmax-qmin)/swidth)
    print 'Shannon-Nyqist data: %g channels, %g width' % (nbins, swidth)
    nu = nbins - nparams
    print 'Estimated degrees of freedom nu: %g' % nu
    #qlist contains the lists of q's for each single shannel channel
    qlist = []
    #q contains the q's for a single channel.
    q = []
    for i in range(nq):
        if( x1[i] - qmin <= swidth):
            q.append(i)
        else:
            print ' = = points in bin %g : %g' % (len(qlist), len(q))
            qlist.append(q)
            q = []
            qmin=x1[i]

    # Now obtain k estimates of chi2_reduced
    chilist = []
    for k in range(nk):
        chi2 = 0
        for qs in qlist:
            ind=random.randint(0,len(qs)-1)
            i=qs[ind]
            thisy2=f*y2[i]+c
            if (bLog):
                tmp  = (math.log(y1[i])-math.log(thisy2))**2
                sig2 = (f*y2sig[i]/thisy2)**2
            else:
                tmp  = (y1[i]-thisy2)**2
                sig2 = (f*y2sig[i])**2
            #Debug
            if ( len(chilist) == 0 ):
                print 'Debug: %g %g' % (x1[i],  tmp/sig2)
            chi2 += tmp/sig2

        chi2red = chi2 / nu
        chilist.append(chi2red)

    list.sort(chilist)

    return chilist

#####################################
# MAIN PROGRAM ######################

random.seed();

#scriptname=os.path.basename(__file__)
parser = argparse.ArgumentParser(description="Fits two curves (the second to the first) by Powell minimisation of the mean square difference.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('files', metavar='N', type=str, nargs='*',
                    help="The files to be manipulated and fitted. The first file will supply the basis-X on which other files will be interpolated.")
parser.add_argument('-mode', type=int, default=0, help="Mode 0: Fit all other files to the first"
                    "Mode 1: Fit first file to all others. Mode 2: Fit all files to each other (!).")
parser.add_argument('-o', type=str, dest='outpref', default='fitted', help='Output prefix. Files will be written as <out>-fittedExp.xvg, etc. ')
parser.add_argument('-debug', dest='bDebug', action='store_true', help='Debug mode.')
parser.add_argument('-v', dest='bVerbose', action='store_true', help='Verbose mode.')
parser.add_argument('-qmin', type=float, default=0.0, help='Minimum x to include in fitting.')
parser.add_argument('-qmax', type=float, default=10.0, help='Maximum x to include in fitting.')
parser.add_argument('-qrange', type=str, default='', help='X-range to include in fitting, overrides above. Give as a 2-tuple, e.g. (0,3).')
parser.add_argument('-lin', dest='bLin', action='store_true', help='Switch to linear-scale fit, instead of log-scale fit.')
parser.add_argument('-now', dest='bNoWeights', action='store_true', help='Do not weight each point by the errors in column 3.')
parser.add_argument('-noc', dest='bNoConst', action='store_true', help='Do not use constant background subtraction.')
parser.add_argument('-diffspec', dest='bDiffSpec', action='store_true', help='Write the difference spectra.')
parser.add_argument('-c0', type=float, default=np.nan, help='Start fitting with given background constant C0 instead of by estimation.')
parser.add_argument('-f0', type=float, default=np.nan, help='Start fitting with given scaling constant F0 instead of by estimation.')
parser.add_argument('-dmax', type=float, default=0.0, help='Additionally calculate Tainer\'s chi_free,'
                             'with this given Dmax for use in Shannon channel determination.')
parser.add_argument('-nk', type=int, default=500, help='In Tainer\'s chi_free modelling, the number of replicates to determine median chi.')
#parser.add_argument('-xmult', type=str, help='X-value pre-multiplier for each imput file.')

args = parser.parse_args()

fitMode = args.mode
outpref=args.outpref

dmax=args.dmax
nk=args.nk
nparams=2


fInit   = args.f0
bEstimateF0 = np.isnan(fInit)
f0EstInterval = 5

cInit   = args.c0
bEstimateC0 = np.isnan(cInit)

bVerbose = args.bVerbose
bUseWeights = not args.bNoWeights
bLog     = not args.bLin
bDebug   = args.bDebug
bNoConst = args.bNoConst
bDiffSpec = args.bDiffSpec

if (args.qrange != ''):
    tmp = make_tuple(args.qrange)
    qmin=tmp[0] ; qmax=tmp[1]
else:
    qmin=args.qmin ; qmax=args.qmax

# Read all files now.
fileList=args.files
nFiles = len(fileList)
if nFiles < 2:
    print >> sys.stderr, '= = ERROR: Script requires at least two files for fitting.'
    sys.exit(1)

# Read and assemble interpolated graphs
# dataRaw is list of numpy arrays, not guaranteed to be of the same sizes.
dataRaw = read_xvgs_all(fileList)
qBasis = find_subrange_all( dataRaw, (qmin, qmax) )
if bVerbose:
    print >> sys.stderr, '= = Input x1 trimmed from ( %g , %g ) to ( %g , %g ) to preserve coverage.' \
                          % ( dataRaw[0][0,0], dataRaw[0][0,-1], qBasis[0], qBasis[-1])

# dataBlock is the trimmed single numpy 3D array.
# It is of arrangement( file, y&dy&etc., vals )
dataBlock = build_interpolated_block( qBasis, dataRaw )


# run through fitting modes.
dataModel = []
dataModelInterp = np.zeros( dataBlock.shape )
# Put original default data into this second block.
dataModelInterp[0] = dataBlock[0]
if fitMode != 2:
    # Fit data between 1st and all others.
    for i in range(1,len(fileList)):
        if bVerbose:
            print >> sys.stderr, "= = = Beginning analysis of file %i (%s)..." % (i, fileList[i] )

        if bEstimateF0:
            fInit=np.mean( dataBlock[0,0,0:f0EstInterval] ) / np.mean( dataBlock[i,0,0:f0EstInterval] )
        if bVerbose:
            print >> sys.stderr, '      ...estimated F0 to be %g' % fInit
        if bEstimateC0:
            cInit=0
            y2min = np.min( dataBlock[i,0] )
            if (fInit*y2min + cInit < 0 ):
                print '= = WARNING: Adjusted initial c to prevent legative logarithms.'
                cInit = 1 - fInit*y2min
            if bVerbose:
                print >> sys.stderr, '      ...estimated C0 to be %g' % cInit

        # Initial values done
        if not bNoConst:
            pos  =  [ fInit, cInit ]
        else:
            pos = [ fInit ]

        if fitMode == 0:
            fminOut = fmin_powell(intensityDiff, pos, args=( dataBlock[0], dataBlock[i], bLog, bUseWeights, bNoConst), full_output=True)
        elif fitMode == 1:
            fminOut = fmin_powell(intensityDiff, pos, args=( dataBlock[i], dataBlock[0], bLog, bUseWeights, bNoConst), full_output=True)
        else:
            print >> sys.stderr, "= = = ERROR: Impossible fit mode found, cannot proceed! (%i) " % fitMode
            sys.exit(20)

        if not bNoConst:
            xopt    = fminOut[0]
            funcopt = fminOut[1]
            f, c = fminOut[0]
        else:
            xopt    = fminOut[0]
            funcopt = fminOut[1]
            f = fminOut[0]
            c = 0

        if fitMode == 0:
            dataModel.append( np.stack( (dataRaw[i][0], f*dataRaw[i][1]+c, f*dataRaw[i][2] ), axis=0 ) )
            dataModelInterp[i,0] = f*dataBlock[i,0]+c
            dataModelInterp[i,1] = f*dataBlock[i,1]
        elif fitMode == 1:
            dataModel.append( np.stack( (dataRaw[0][0], f*dataRaw[0][1]+c, f*dataRaw[0][2] ), axis=0 ) )
            dataModelInterp[i,0] = f*dataBlock[0,0]+c
            dataModelInterp[i,1] = f*dataBlock[0,1]

        if bDebug:
            outFile='%s-Debug-%i.xvg' % (outpref, i)
            header='# Debug file. Fit mode %i\n# Original file A: %s\n# Original file B: %s\n# Optimised Params: %s' \
                    % (fitMode, fileList[0], fileList[i], xopt)
            write_xvg( outFile, qBasis, dataModelInterp[0:i+1:i], header)

        if bDiffSpec:
            outFile='%s-DiffSpec-%i.xvg' % (outpref, i)
            header='# Difference spectra. Fit mode %i\n# Original file A: %s\n# Original file B: %s\n# Optimised Params: %s' \
                    % (fitMode, fileList[0], fileList[i], xopt)
            diffSpec = subtract_spectra( dataModelInterp[0], dataModelInterp[i])
            write_xvg( outFile, qBasis, diffSpec, header)


        #Now calculate Tainer's chi_free measure for the fitted
        if ( dmax != 0 ):
            #fminOut = fmin_powell(intensityDiff, pos, args=( dataBlock[0], dataBlock[i], bLog, bUseWeights, bNoConst), full_output=bVerbose)
            #chilist = intChiFree(x, y1int, y2int, y2errint, f, c, bLog, dmax, nk, nparams)
            if fitMode == 0:
                chiList = intChiFree(qBasis, dataBlock[0,0], dataModelInterp[i,0], dataModelInterp[i,1], 1, 0, bLog, dmax, nk, nparams)
            elif fitMode == 1:
                chiList = intChiFree(qBasis, dataModelInterp[i,0], dataBlock[0,0], dataModelInterp[0,1], 1, 0, bLog, dmax, nk, nparams)

            outFile='%s-ChiFree-%i.xvg' % (outpref, i)
            header='# Chi-free. Fit mode %i\n# Original file A: %s\n# Original file B: %s\n# Chi_free (median): %g' \
                    % (fitMode, fileList[0], fileList[i], math.sqrt(chiList[nk/2]) )
            write_xvg( outFile, np.arange(nk), np.sqrt(chiList), header, bType=True)

            #Print output when done with curve.
            outFile='%s-Overlay-%i.xvg' % (outpref, i)
            header='# Spectral overlay. Fit mode %i\n# Original file A: %s\n# Original file B: %s\n# Optimised Params: %s' \
                    % (fitMode, fileList[0], fileList[i], xopt)
            if fitMode == 0:
                write_xvg_raw( outFile, list( [dataRaw[0], dataModel[-1]] ), header, bReverse=True)
            elif fitMode == 1:
                write_xvg_raw( outFile, list( [dataRaw[i], dataModel[-1]] ), header, bReverse=True)

        fp = open('%s-Exp-%i.xvg' % (outpref, i), 'w')
        print >> fp, '#npts_fit = %8i'  % len(qBasis)
        print >> fp, '#qmin_fit = %8g' % qBasis[0]
        print >> fp, '#qmax_fit = %8g' % qBasis[-1]
        print >> fp, '#f    = %12g' % f
        print >> fp, '#c    = %12g' % c
        print >> fp, '#chi2 = %12g' % funcopt
        print >> fp, '#chi  = %12g' % math.sqrt(funcopt)
        if ( dmax != 0 ):
            print >> fp, '#chiFree = %12g' % math.sqrt( chilist[nk/2] )
        write_to_xvg(fp, dataModel[-1][0], dataModel[-1][1], dataModel[-1][2])
        fp.close()

# Write chi matrix.
elif fitMode == 2:
    chiMatrix = np.zeros( (nFiles, nFiles) )
    fMatrix = np.zeros( (nFiles, nFiles) )
    for pair in itertools.permutations( np.arange(nFiles), r=2 ):
        i, j = pair

        if bEstimateF0:
            fInit=np.mean( dataBlock[i,0,0:f0EstInterval] ) / np.mean( dataBlock[j,0,0:f0EstInterval] )
        if bVerbose:
            print >> sys.stderr, '      ...estimated F0 to be %g' % fInit
        if bEstimateC0:
            cInit=0
            y2min = np.min( dataBlock[j,0] )
            if (fInit*y2min + cInit < 0 ):
                print '= = WARNING: Adjusted initial c to prevent legative logarithms.'
                cInit = 1 - fInit*y2min
            if bVerbose:
                print >> sys.stderr, '      ...estimated C0 to be %g' % cInit
        if not bNoConst:
            pos  =  [ fInit, cInit ]
        else:
            pos = [ fInit ]

        fminOut = fmin_powell(intensityDiff, pos, args=( dataBlock[i], dataBlock[j], bLog, bUseWeights, bNoConst), full_output=True)
        if bVerbose:
            print >> sys.stderr, '    ...minimisation results:'
            print >> sys.stderr, fminOut
        else:
            print >> sys.stderr, '    ....chi(%i,%i)^2: %f' % (i, j, fminOut[1])
        chiMatrix[i,j] = math.sqrt(fminOut[1])
        if not bNoConst:
            fMatrix[i,j] = fminOut[0][0]
        else:
            fMatrix[i,j] = fminOut[0]

    print qBasis[0], dataBlock[:,0,0], dataBlock[:,1,0]
    print fMatrix
    print chiMatrix
    header='# Chi materix between files: %s\n# npts_fit = %8i\n# qmin_fit = %8g\n# qmax_fit = %8g' % \
            ( str(fileList), len(qBasis), qBasis[0], qBasis[-1] )
    outFile='%s-chiMatrix.dat' % outpref
    write_xvg_raw( outFile, chiMatrix, header, bReverse=False, bAmpersand=False)
