#!/usr/bin/python
import sys, math, itertools
import numpy as np
from scipy.optimize import fmin_powell
import random
import argparse
from ast import literal_eval as make_tuple
# Used to call datgnom and outsource Dmax determination.
import general_scripts as gs
import saxscalc as sc

def read_xvgs_all(filelist):
    out_list=[]
    for i in range(len(filelist)):
        block = np.array( gs.load_xydy(filelist[i]) )
        print( "  ...file read, shape: %s" % str(block.shape), file=sys.stderr )
        out_list.append(block)
    return out_list

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
    y1    = args[0][0] ; y1sig = args[0][1]
    y2    = args[1][0] ; y2sig = args[1][1]
    fitMetric   = args[2]
    bUseWeights = args[3]
    bNoConst    = args[4]
    stride      = args[5]
    numRounds   = args[6]
    if fitMetric == 'chi' or fitMetric == 'chi_red':
        if bNoConst:
            f = pos[0] ; c=0.0
        else:
            f, c = pos
        if bUseWeights:
            value=sc.chi_square(y1,f*y2+c, dx1=y1sig, dx2=f*y2sig)
        else:
            value=sc.chi_square(y1,f*y2+c)
    elif fitMetric == 'log_chi':
        if bNoConst:
            f = pos[0] ; c=0.0
        else:
            f, c = pos
        # Check for negative values first.
        if np.any(f*y2+c < 0.0):
            print( "= = WARNING: for values of f %g and c %g there exists invalid logs." % (f, c), file=sys.stderr )
            return 1e20
        if bUseWeights:
            value=sc.log_chi_square(y1,f*y2+c, dx1=y1sig, dx2=f*y2sig)
        else:
            value=sc.log_chi_square(y1,f*y2+c)
    elif fitMetric == 'chi_free':
        if bNoConst:
            f = pos[0] ; c=0.0
            nParams=1
        else:
            f, c = pos
            nParams=2
        if bUseWeights:
            value = sc.chi_square_free(y1, f*y2+c, dx1=y1sig, dx2=f*y2sig, stride=stride, nParams=nParams, nRounds=numRounds)
        else:
            value = sc.chi_square_free(y1, f*y2+c, stride=stride, nParams=nParams, nRounds = numRounds)
    elif fitMetric == 'vr':
        c = pos
        value = sc.volatility_ratio(y1-c, y2+c, stride=stride )
        # If the entire strip is masked return 1e99. This usually indicates some form of domain error in the underlying calculation.
        if np.ma.is_masked(value) or np.isnan(value):
            value=1e99
    elif fitMetric == 'cormap' or fitMetric == 'cormap_matrix':
        # Will need a two-stage eliminator, perhaps, but essentially we'll need a system that minimises the number
        # of sequential runs.
        f, c = pos
        runs = sc.run_distribution( y1 > f*y2+c )
        value = len(runs)
        # prob = sc.cormap_value( y1, f*y2+c )
        # value = 1.0 - sc.cormap_value( y1, f*y2+c )
    else:
        print( "= = = ERROR, metric not recognised! %s" % fitMetric, file=sys.stderr )
        sys.exit(1)
    return value

def populationIntensityDiff(pos, *args):
    yT = args[0][0,0]  ; yTsig = args[0][0,1]
    yP = args[0][1:,0] ; yPsig = args[0][1:,1]
    fitMetric   = args[1]
    bUseWeights = args[2]
    bNoConst    = args[3]
    stride      = args[4]
    numRounds   = args[5]

    if fitMetric == 'vr':
        # = = Try to maintain symmetry.
        c = pos[-1]
        yT -= c 
        y2 = np.mean(yP, axis=0)+c ; y2sig = np.mean(yPsig, axis=0)
    elif bNoConst:
        f = pos ; c=0.0
        y2 = np.mean(f[:,None]*yP, axis=0)+c ; y2sig = np.mean(f[:,None]*yPsig, axis=0)
    else:
        f = pos[:-1] ; c = pos[-1]
        y2 = np.mean(f[:,None]*yP, axis=0)+c ; y2sig = np.mean(f[:,None]*yPsig, axis=0)

    # Don't suport chi_red just yet.
    if fitMetric == 'chi':
        if bUseWeights:
            value=sc.chi_square(yT,y2, dx1=yTsig, dx2=y2sig)
        else:
            value=sc.chi_square(yT,y2)
    elif fitMetric == 'log_chi':
        # Check for negative values first.
        if np.any(y2 < 0.0):
            #print( "= = WARNING: for values of f %g and c %g there exists invalid logs." % (f, c), file=sys.stderr )
            return 1e20
        if bUseWeights:
            value=sc.log_chi_square(yT,y2, dx1=yTsig, dx2=y2sig)
        else:
            value=sc.log_chi_square(yT,y2)
    elif fitMetric == 'chi_free':
        nParams=len(pos)
        if bUseWeights:
            value = sc.chi_square_free(yT, y2, dx1=yTsig, dx2=y2sig, stride=stride, nParams=nParams, nRounds=numRounds)
        else:
            value = sc.chi_square_free(yT, y2, stride=stride, nParams=nParams, nRounds = numRounds)
    elif fitMetric == 'vr':
        value = sc.volatility_ratio(yT, y2, stride=stride )
    elif fitMetric == 'cormap' or fitMetric == 'cormap_matrix':
        runs = sc.run_distribution( yT > y2 )
        value = len(runs)
        # prob = sc.cormap_value( yT, y2 )
        # value = 1.0 - sc.cormap_value( yT, y2 )
    else:
        print( "= = = ERROR, metric not recognised! %s" % fitMetric, file=sys.stderr )
        sys.exit(1)
    print( value )
    return value

#####################################
# MAIN PROGRAM ######################

random.seed();

#scriptname=os.path.basename(__file__)
parser = argparse.ArgumentParser(description="Fits a set of curves by minimising differences according to some popular metrics,"
                        "utilising Powell-minimisation (sequential step-wise search over each variable) over one of two of the "
                        "variables overall scaling *f* and constant subtraction *c*. "
                        "These two variables approximate uncertainties in sample concentration/beam intensities and basic buffer subtraction."
                        "Note that there are some combinations of arguments that won't work well, like cormap with Powell minimization.",
                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('files', metavar='N', type=str, nargs='*',
                    help="The list of files to be manipulated and fitted. For modes 0 and 1, the first file given acts as the reference point"
                        "by which other files are comparied with. The first file will also serve to define the set of independent-variable"
                        "(scattering-angle q) to act as the basis for interpolation.")
parser.add_argument('-mode', type=int, default=0, help="Determine which curves to fit with which. Mode 0: Fit all other files to the first"
                        "Mode 1: Fit first file to all others. Mode 2: Fit all files to each other (!)."
                        "Mode 3: Use population mode. Fit other files as a population to the first."
                        "(ToDo) Mode 4: Use population mode. Use the first again as the target, fit combination of others." )
parser.add_argument('-metric', type=str, default='chi', help="Determine the metric to use as comparison."
                        "Options are: chi | log_chi | chi_free | volatility | cormap | cormap_matrix "
                        "chi_free is as defined by Rambo and Tainer, Nature, 2013."
                        "Volatility is the Volatility Ratio as defined by Hura et al. (2013). Can be given as 'vr' instead for short."
                        "cormap is the Correlation Map as defined by Franke et al. (2015)")
parser.add_argument('-o', type=str, dest='outpref', default='fitted', help='Output prefix. Files will be written as <out>-fittedExp.xvg, etc. ')
parser.add_argument('-debug', dest='bDebug', action='store_true', help='Debug mode.')
parser.add_argument('-v', dest='bVerbose', action='store_true', help='Verbose mode.')
parser.add_argument('-diffspec', dest='bDiffSpec', action='store_true', help='Write the difference spectra after fitting.')

# = = = Parameters for all fitting modes.
parser.add_argument('-qmin', type=float, default=0.0, help='Minimum q to include in fitting.')
parser.add_argument('-qmax', type=float, default=10.0, help='Maximum q to include in fitting.')
parser.add_argument('-qrange', type=str, default='', help='Q-range to include in fitting, overrides above. Give as a 2-tuple, e.g. (0,3).')

# = = = Parameters for fitting modes.
parser.add_argument('-now', dest='bNoWeights', action='store_true', help='Do not weight each point by the errors in column 3.')
parser.add_argument('-noc', dest='bNoConst', action='store_true', help='Do not use constant background subtraction.')
parser.add_argument('-c0', type=float, default=np.nan, help='Start fitting with given background constant C0 instead of by estimation.')
parser.add_argument('-f0', type=float, default=np.nan, help='Start fitting with given scaling constant F0 instead of by estimation.')

# = = = Parameters for chi_free and/or vr
parser.add_argument('-Dmax', type=float, default=np.nan, help='Give the maximum molecular extent D_max, required to compute number of Shannon channels '
                    'as used in chi_free and volatility ratio computations. NB: can compute using other utilities, e.g. ATSAS datgnom.')
parser.add_argument('-nRounds', type=int, default=500, help='In Tainer\'s chi_free modelling, the number of replicates to determine median chi.')
#parser.add_argument('-xmult', type=str, help='X-value pre-multiplier for each imput file.')

args = parser.parse_args()

fitMode   = args.mode
fitMetric = args.metric
outpref=args.outpref

Dmax=args.Dmax
numRounds=args.nRounds
nparams=2

fInit   = args.f0
bEstimateF0 = np.isnan(fInit)
f0EstInterval = 5

cInit   = args.c0
bEstimateC0 = np.isnan(cInit)

bVerbose = args.bVerbose
bUseWeights = not args.bNoWeights
bDebug   = args.bDebug
bNoConst = args.bNoConst
bDiffSpec = args.bDiffSpec

if (args.qrange != ''):
    tmp = make_tuple(args.qrange)
    qmin=tmp[0] ; qmax=tmp[1]
else:
    qmin=args.qmin ; qmax=args.qmax

# = = = Sanitise fitMetric and check for inconsistent given arguments.
fitMetric = fitMetric.lower()
if fitMetric == 'v_r' or fitMetric == 'volatrat' or fitMetric == 'volatility':
    fitMetric = 'vr'
elif fitMetric == 'chifree':
    fitMetric = 'chi_free'
elif fitMetric == 'log' or fitMetric == 'logchi':
    fitMetric = 'log_chi'
elif fitMetric == 'chired' or fitMetric == 'chi_reduced':
    fitMetric = 'chi_red'

if fitMetric == 'vr':
    if bNoConst:
        print( '= = = ERROR: The volatility ratio has only 1 free parameter for constant subtraction. Cannot be used with argument -noc.', file=sys.stderr )
        sys.exit(1)
    if np.isnan(Dmax):
        print( '= = = ERROR: Volatility ratio requires the definition of Dmax in order to operate!', file=sys.stderr )
        sys.exit(1)

if fitMetric == 'chi_free' or fitMetric == 'chi_red':
    if np.isnan(Dmax):
        print( '= = = ERROR: Chi_free and chi_red ratios requires the definition of Dmax in order to operate!', file=sys.stderr )
        sys.exit(1)

if fitMetric == 'cormap_matrix' and fitMode == 2:
    print( '= = = ERROR: Refuse to do CorMap matrices for a matrix of fits. This is not currently supported.', file=sys.stderr )
    sys.exit(1)

# Read all files now.
fileList=args.files
nFiles = len(fileList)
if nFiles < 2:
    print( '= = ERROR: Script requires at least two files for fitting.', file=sys.stderr )
    sys.exit(1)

# Read and assemble interpolated graphs
# dataRaw is list of numpy arrays, not guaranteed to be of the same sizes.
dataRaw = read_xvgs_all(fileList)
qBasis = find_subrange_all( dataRaw, (qmin, qmax) )
if bVerbose:
    print( '= = Input x1 trimmed from ( %g , %g ) to ( %g , %g ) to preserve coverage.' \
            % ( dataRaw[0][0,0], dataRaw[0][0,-1], qBasis[0], qBasis[-1]), file=sys.stderr )
    print( '= = ...number of q-points remaining: %i' % len(qBasis), file=sys.stderr )

# Assume again that the q-points are evenly spaced. Determine number of point per Shannon channel
if not np.isnan(Dmax):
    deltaQ = qBasis[1] - qBasis[0]
    numPointsPerChannel = int( np.pi / Dmax / deltaQ )
    numChannels = (qBasis[-1]-qBasis[0])*Dmax/np.pi
    #1.0*len(qBasis)/numPointsPerChannel
    if bVerbose:
        print( '= = Determined the number of Channels to be %g based on input data, resulting in %i points per channel' \
                % ( numChannels,  numPointsPerChannel ), file=sys.stderr )
    if bNoConst:
        redFactor=numChannels/(numChannels-1)
    else:
        redFactor=numChannels/(numChannels-2)
else:
    numPointsPerChannel = 1

# dataBlock is the trimmed single numpy 3D array.
# It is of arrangement( file, y&dy&etc., vals )
dataBlock = build_interpolated_block( qBasis, dataRaw )
if bVerbose:
    print( '= = Interpolated data block has been built with size:', dataBlock.shape, file=sys.stderr )

# run through fitting modes. dataModel is the converted raw data, while dataModelInterp is the converted interpolated data.
dataModel = []
dataModelInterp = np.zeros( dataBlock.shape )
# Put original default data into this second block.
dataModelInterp[0] = dataBlock[0]

if fitMode < 2:
    # = = = Fit data between 1st and all others. = = =
    for i in range(1,len(fileList)):
        if bVerbose:
            print( "= = = Beginning analysis of file %i (%s)..." % (i, fileList[i] ), file=sys.stderr )

        if bEstimateF0:
            fInit=np.mean( dataBlock[0,0,0:f0EstInterval] ) / np.mean( dataBlock[i,0,0:f0EstInterval] )
            if fitMode == 1:
                fInit=1/fInit
        if bVerbose:
            print( '      ...estimated F0 to be %g' % fInit, file=sys.stderr )
        if bEstimateC0:
            cInit=0
            y2min = np.min( dataBlock[i,0] )
            if fitMetric == 'log_chi' and (fInit*y2min + cInit < 0 ):
                print( '= = WARNING: Adjusted initial c to prevent legative logarithms.' )
                cInit = 1 - fInit*y2min
            if bVerbose:
                print( '      ...estimated C0 to be %g' % cInit, file=sys.stderr )

        #if fitMetric == 'chi' or fitMetric == 'chi_free' or fitMetric == 'cormap':
        if fitMetric == 'vr':
            pos = [ cInit ]
        elif not bNoConst:
            pos = [ fInit, cInit ]
        else:
            pos = [ fInit ]

        if fitMode == 0:
            fminOut = fmin_powell(intensityDiff, pos, args=( dataBlock[0], dataBlock[i], \
                    fitMetric, bUseWeights, bNoConst, numPointsPerChannel, numRounds), full_output=True)
        elif fitMode == 1:
            fminOut = fmin_powell(intensityDiff, pos, args=( dataBlock[i], dataBlock[0], \
                    fitMetric, bUseWeights, bNoConst, numPointsPerChannel, numRounds), full_output=True)
        else:
            print( "= = = ERROR: Impossible fit mode found, cannot proceed! (%i) " % fitMode, file=sys.stderr )
            sys.exit(20)

        xopt = fminOut[0] ; funcopt = fminOut[1]

        if fitMetric == 'vr':
            c = fminOut[0]
            # Here it's simply defined as the fill ratio.
            f = np.mean(dataBlock[0][0]/dataBlock[i][0])
            if fitMode==1:
                f = 1.0/f
        elif not bNoConst:
            xopt = fminOut[0] ; funcopt = fminOut[1]
            f, c = fminOut[0]
        else:
            xopt = fminOut[0] ; funcopt = fminOut[1]
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
            gs.print_xvg_complex( outFile, qBasis, dataModelInterp[0:i+1:i], header)

        if bDiffSpec:
            outFile='%s-DiffSpec-%i.xvg' % (outpref, i)
            header='# Difference spectra. Fit mode %i\n# Original file A: %s\n# Original file B: %s\n# Optimised Params: %s' \
                    % (fitMode, fileList[0], fileList[i], xopt)
            diffSpec = subtract_spectra( dataModelInterp[0], dataModelInterp[i])
            gs.print_xvg_complex( outFile, qBasis, diffSpec, header)

        outFile='%s-Exp-%i.xvg' % (outpref, i)
        fileHeader=[]
        fileHeader.append( '#npts_fit = %8i'  % len(qBasis) )
        fileHeader.append( '#qmin_fit = %8g' % qBasis[0]    )
        fileHeader.append( '#qmax_fit = %8g' % qBasis[-1]   )
        fileHeader.append( '#f    = %12g' % f               )
        fileHeader.append( '#c    = %12g' % c               )
        if fitMetric == 'chi' or fitMetric == 'log_chi':
            fileHeader.append( '#chi2 = %12g' % funcopt )
            fileHeader.append( '#chi  = %12g' % math.sqrt(funcopt) )
        elif fitMetric == 'chi_red':
            fileHeader.append( '#chi2 = %12g' % (funcopt*redFactor*redFactor) )
            fileHeader.append( '#chi  = %12g' % (math.sqrt(funcopt)*redFactor) )
        elif fitMetric == 'chi_free':
            fileHeader.append( '#chi2Free = %12g' % funcopt )
            fileHeader.append( '#chiFree  = %12g' % math.sqrt( funcopt ) )
            fileHeader.append( '#DMaxPar  = %12g' % Dmax )
            fileHeader.append( '#nPoints  = %i' % numPointsPerChannel )
            fileHeader.append( '#nRounds  = %i' % numRounds )
        elif fitMetric == 'vr':
            fileHeader.append( '#volatRat = %12g' % funcopt )
            fileHeader.append( '#DMaxPar  = %12g' % Dmax )
            fileHeader.append( '#nPoints  = %i' % numPointsPerChannel )
            # = = = Hack: report the individual ratio bins by recomputing the volatility ratio
            # value = sc.volatility_ratio(y1-c, y2+c, stride=stride )
            if   fitMode == 0:
                ratios = sc.volatility_ratio( dataBlock[0][0]-c, dataBlock[i][0]+c, stride=numPointsPerChannel, bElementwise=True)
            elif fitMode == 1:
                ratios = sc.volatility_ratio( dataBlock[i][0]-c, dataBlock[0][0]+c, stride=numPointsPerChannel, bElementwise=True)
            for nR, R in enumerate(ratios):
                fileHeader.append( '#VRBins %i   = %g' % (nR, R) )

        elif fitMetric == 'cormap':
            fileHeader.append( '#maxRun  = %i' % funcopt )
            prob = sc.probability_cormap_either( len(qBasis), funcopt )
            if prob > 0.0:
                log10p = -1*np.log10(prob)
            else:
                log10p = 99
            fileHeader.append( '#probExc = %g' % prob )
            fileHeader.append( '#-log10P = %g' % log10p )
            print( "= = Probability that difference is due completely to random noise = %g"  % prob )
        elif fitMetric == 'cormap_matrix':
            if fitMode == 0:
                mat = sc.cormap_matrix( dataBlock[0][0], dataModelInterp[i][0] )
            else:
                mat = sc.cormap_matrix( dataBlock[i][0], dataModelInterp[i][0] )
            outFile='%s-CorMapMatrix-%i.dat' % (outpref, i)
            np.savetxt( outFile, mat )
            # Break out into the next combo!
            continue
        else:
            print( "= = = ERROR, metric not recognised! %s" % fitMetric, file=sys.stderr )
            sys.exit(1)

        gs.print_xy(outFile, dataModel[-1][0], dataModel[-1][1], dy=dataModel[-1][2], header=fileHeader)

# Write chi matrix.
elif fitMode == 2:
    valueMatrix = np.zeros( (nFiles, nFiles) )
    fMatrix = np.zeros( (nFiles, nFiles) )

    for pair in itertools.permutations( np.arange(nFiles), r=2 ):
        i, j = pair
        if bEstimateF0:
            fInit=np.mean( dataBlock[i,0,0:f0EstInterval] ) / np.mean( dataBlock[j,0,0:f0EstInterval] )
        if bVerbose:
            print( '      ...estimated F0 to be %g' % fInit, file=sys.stderr )
        if bEstimateC0:
            cInit=0
            y2min = np.min( dataBlock[j,0] )
            if fitMetric == 'log_chi' and (fInit*y2min + cInit < 0 ):
                print( '= = WARNING: Adjusted initial c to prevent legative logarithms.' )
                cInit = 1 - fInit*y2min
            if bVerbose:
                print( '      ...estimated C0 to be %g' % cInit, file=sys.stderr )

        if fitMetric == 'vr':
            pos = [ cInit ]
        elif not bNoConst:
            pos  = [ fInit, cInit ]
        else:
            pos = [ fInit ]

        fminOut = fmin_powell(intensityDiff, pos, args=( dataBlock[i], dataBlock[j], \
            fitMetric, bUseWeights, bNoConst, numPointsPerChannel, numRounds), full_output=True)
        if bVerbose:
            print( '    ...minimisation results:', file=sys.stderr )
            print( fminOut, file=sys.stderr )
        else:
            print( '    ....value(%i,%i): %g' % (i, j, fminOut[1]), file=sys.stderr )
        # = = = Enter final value into matrix.
        if fitMetric == 'chi' or fitMetric == 'log_chi' or fitMetric == 'chi_free' :
            valueMatrix[i,j] = math.sqrt(fminOut[1])
        elif fitMetric == 'chi_red':
            valueMatrix[i,j] = math.sqrt(fminOut[1])*redFactor
        elif fitMetric == 'vr':
            valueMatrix[i,j] = fminOut[1]
        elif fitMetric == 'cormap':
            #valueMatrix[i,j] = fminOut[1]
            prob = sc.probability_cormap_either( len(qBasis), fminOut[1] )
            if prob > 0.0:
                valueMatrix[i,j] = -np.log10(prob)
            else:
                valueMatrix[i,j] = 20.0
        else:
            print( "= = ERROR: fitMetric not recognised in loop! This is a coding error.", file=sys.stderr )
            sys.exit(1)

        if fitMetric == 'vr':
            fMatrix[i,j] = np.mean(dataBlock[i][0]/dataBlock[j][0])
        elif not bNoConst:
            fMatrix[i,j] = fminOut[0][0]
        else:
            fMatrix[i,j] = fminOut[0]

    #print( qBasis[0], dataBlock[:,0,0], dataBlock[:,1,0] )
    #print( fMatrix )
    print( valueMatrix )
    header='# Chi matrix between files: %s\n# npts_fit = %8i\n# qmin_fit = %8g\n# qmax_fit = %8g' % \
            ( str(fileList), len(qBasis), qBasis[0], qBasis[-1] )
    outFile='%s_matrix.dat' % outpref
    gs.print_xvg_matrix( outFile, valueMatrix, header, bReverse=False, bAmpersand=False)

if fitMode == 3:
    # = = = Begin population fit model.
    # There will be N scaling factor, plus an optional 1 constant subtraction.
    # chi ~ ( I_exp - (f1*I1 + f2*I2 + ... + c) )
    if bVerbose:
        print( "= = = Starting population modelling...", file=sys.stderr )

    if len(fileList) < 3:
        print( "= = = ERROR: Population fitting requires at last two component curves and one target curve!", file=sys.stderr )
        sys.exit(1)

    if bEstimateF0:
        fInit = np.zeros(len(fileList)-1, dtype=dataBlock.dtype )
        for i in range(1,len(fileList)):
            if bEstimateF0:
                fInit[i-1] =  np.mean( dataBlock[0,0,0:f0EstInterval] ) / np.mean( dataBlock[i,0,0:f0EstInterval] )
            if bVerbose:
                print( '      ...estimated F0 to be %g' % fInit[-1], file=sys.stderr )
    else:
        fInit = np.ones( len(fileList)-1, dtype=dataBlock.dtype )
    pos = fInit/(len(fileList)-1.0)

    cInit = 0
    #if bEstimateC0:
    #    y2min = np.min( dataBlock[i,0] )
    #    if fitMetric == 'log_chi' and (fInit*y2min + cInit < 0 ):
    #        print( '= = WARNING: Adjusted initial c to prevent legative logarithms.' )
    #        cInit = 1 - fInit*y2min
    #    if bVerbose:
    #        print( '      ...estimated C0 to be %g' % cInit, file=sys.stderr )

    #if fitMetric == 'chi' or fitMetric == 'chi_free' or fitMetric == 'cormap':
    if fitMetric == 'vr':
        # Remove one scaling factor variable since vr autoscales, and set all ratios to be equal at first.
        pos = np.append( np.ones(len(fileList)-2, dtype=dataBlock.dtype ), cInit )
    elif not bNoConst :
        pos = np.append( pos, cInit )

    fminOut = fmin_powell(populationIntensityDiff, pos, args=( dataBlock, \
                    fitMetric, bUseWeights, bNoConst, numPointsPerChannel, numRounds), full_output=True)

    print( xopt )

    # = = = Parse results
    xopt = fminOut[0] ; funcopt = fminOut[1]
    yTarget = dataBlock[0,0]
    if fitMetric == 'vr':
        # The fit for volatility ratios is different from the others.
        f = np.insert( fminOut[0][:-1], 0, 1.0 )
        c = fminOut[0][-1]
        fact = sc.volatility_ratio_scaling( yTarget , np.mean( f[:,None]*dataBlock[1:,0], axis=0 ) + c )
        f *= fact ; c*= fact
        yModel  = np.mean( f[:,None]*dataBlock[1:,0], axis=0 ) + c
        dyModel = np.mean( f[:,None]*dataBlock[1:,1], axis=0 )
        #y2 = np.mean(f[:,None]*yP, axis=0)+c ; y2sig = np.mean(f[:,None]*yPsig, axis=0)
    else:
        if not bNoConst:
            f = fminOut[0][:-1] ; c = fminOut[0][-1]
        else:
            xopt = fminOut[0] ; funcopt = fminOut[1]
            f = fminOut[0] ; c = 0
        #dataModel.append( np.stack( (dataRaw[i][0], f*dataRaw[i][1]+c, f*dataRaw[i][2] ), axis=0 ) )
        yModel  = np.mean( f[:,None]*dataBlock[1:,0], axis=0 ) + c
        dyModel = np.mean( f[:,None]*dataBlock[1:,1], axis=0 )

    if bDiffSpec:
        outFile='%s-DiffSpec-pop.xvg' % (outpref)
        header='# Difference spectra. Fit mode %i\n# Original file A: %s\n# Original file B: %s\n# Optimised Params: %s' \
                % (fitMode, fileList[0], fileList[i], xopt)
        diffSpec = subtract_spectra( yTarget, yModel)
        gs.print_xvg_complex( outFile, qBasis, diffSpec, header)

    outFile='%s-Exp-pop.xvg' % (outpref)
    fileHeader=[]
    fileHeader.append( '#npts_fit = %8i'  % len(qBasis) )
    fileHeader.append( '#qmin_fit = %8g' % qBasis[0]    )
    fileHeader.append( '#qmax_fit = %8g' % qBasis[-1]   )
    for i in range(len(f)):
        fileHeader.append( '#f%i    = %12g' % (i+1, f[i] ) )

    fileHeader.append( '#c    = %12g' % c               )
    if fitMetric == 'chi' or fitMetric == 'log_chi':
        fileHeader.append( '#chi2 = %12g' % funcopt )
        fileHeader.append( '#chi  = %12g' % math.sqrt(funcopt) )
    elif fitMetric == 'chi_free':
        fileHeader.append( '#chi2Free = %12g' % funcopt )
        fileHeader.append( '#chiFree  = %12g' % math.sqrt( funcopt ) )
        fileHeader.append( '#DMaxPar  = %12g' % Dmax )
        fileHeader.append( '#nPoints  = %i' % numPointsPerChannel )
        fileHeader.append( '#nRounds  = %i' % numRounds )
    elif fitMetric == 'vr':
        fileHeader.append( '#volatRat = %12g' % funcopt )
        fileHeader.append( '#DMaxPar  = %12g' % Dmax )
        fileHeader.append( '#nPoints  = %i' % numPointsPerChannel )
    elif fitMetric == 'cormap':
        fileHeader.append( '#maxRun  = %i' % funcopt )
        prob = sc.probability_cormap_either( len(qBasis), funcopt )
        if prob > 0.0:
            log10p = -1*np.log10(prob)
        else:
            log10p = 99
        fileHeader.append( '#probExc = %g' % prob )
        fileHeader.append( '#-log10P = %g' % log10p )
        print( "= = Probability that difference is due completely to random noise = %g"  % prob )
    elif fitMetric == 'cormap_matrix':
        mat = sc.cormap_matrix( yTarget, yModel )
        outFile='%s-CorMapMatrix-pop.dat' % (outpref)
        np.savetxt(outFile, mat)
        # Break out into the next combo!
    else:
        print( "= = = ERROR, metric not recognised! %s" % fitMetric, file=sys.stderr )
        sys.exit(1)

    gs.print_xy(outFile, qBasis, yModel, dy=dyModel, header=fileHeader)
