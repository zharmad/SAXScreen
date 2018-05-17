#!/usr/bin/python
import sys
import numpy as np
import argparse
from ast import literal_eval as make_tuple
# Used to call datgnom and outsource Dmax determination.
import saxscalc as sc
import general_scripts as gs

def read_xvgs_all(filelist):
    out_list=[]
    for i in range(len(filelist)):
        block = np.array( gs.load_xydy(filelist[i]) )
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

#####################################
# MAIN PROGRAM ######################

#scriptname=os.path.basename(__file__)
parser = argparse.ArgumentParser(description="Fits two curves (the second to the first) by Powell minimisation of the mean square difference.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('files', metavar='N', type=str, nargs='*',
                    help="The files to be manipulated and fitted. The first file will supply the basis-X on which other files will be interpolated.")
parser.add_argument('-metric', type=str, default='chi', help="Determine the metric to use as comparison."
                        "Options are: chi | log_chi | chi_free | V_R | cormap | cormap_matrix ."
                        "log_chi is the logarithmic version of chi to remove scale dependence in fits."
                        "chi_free is as defined by Rambo and Tainer, Nature, 2013."
                        "V_R is the Volatility Ratio as defined by Hura et al. (2013)."
                        "cormap is the Correlation Map as defined by Franke et al. (2015). No Bonferronni correction has yet been applied.")
parser.add_argument('-qmin', type=float, default=0.0, help='Minimum x to include in fitting.')
parser.add_argument('-qmax', type=float, default=10.0, help='Maximum x to include in fitting.')
parser.add_argument('-qrange', type=str, default='', help='X-range to include in fitting, overrides above. Give as a 2-tuple, e.g. (0,3).')
parser.add_argument('-now', dest='bNoWeights', action='store_true', help='Ignore the uncertainties values (column 3) even if present in data.')
parser.add_argument('-Dmax', type=float, default=np.nan, help='Give the maximum molecular extent D_max, required to compute number of Shannon channels '
                             'as used in chi_free and volatility ratio computations.')
parser.add_argument('-nRounds', type=int, default=500, help='In Tainer\'s chi_free modelling, the number of replicates to determine median chi.')

args = parser.parse_args()

Dmax=args.Dmax
numRounds=args.nRounds
bNoWeights = args.bNoWeights
fitMetric=args.metric

if (args.qrange != ''):
    tmp = make_tuple(args.qrange)
    qmin=tmp[0] ; qmax=tmp[1]
else:
    qmin=args.qmin ; qmax=args.qmax

if fitMetric == 'V_R' or fitMetric=='chi_free':
    if np.isnan(Dmax):
        print >> sys.stderr, '= = = ERROR: A definition of maximum molecular extent Dmax is required use Volatility ratio or chi_free.'
        sys.exit(1)

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

# Assume again that the q-points are evenly spaced. Determine number of point per Shannon channel
if not np.isnan(Dmax):
    deltaQ = qBasis[1] - qBasis[0]
    numPointsPerChannel = int( np.pi / Dmax / deltaQ )
    numChannels = 1.0*len(qBasis)/numPointsPerChannel
else:
    numPointsPerChannel = 1

# dataBlock is the trimmed single numpy 3D array.
# It is of arrangement( file, y&dy&etc., vals )
dataBlock = build_interpolated_block( qBasis, dataRaw )

for i in range(1,nFiles):
    if fitMetric=='chi':
        if bNoWeights:
            value=sc.chi_square(dataBlock[0,0], dataBlock[i,0])
        else:
            value=sc.chi_square(dataBlock[0,0], dataBlock[i,0], dx1=dataBlock[i,1], dx2=dataBlock[i,1])
        value=np.sqrt(value)
    elif fitMetric=='log_chi':
        if bNoWeights:
            value=sc.log_chi_square(dataBlock[0,0], dataBlock[i,0])
        else:
            value=sc.log_chi_square(dataBlock[0,0], dataBlock[i,0], dx1=dataBlock[i,1], dx2=dataBlock[i,1])
        value=np.sqrt(value)
    elif fitMetric == 'chi_free':
        if bNoWeights:
            value = sc.chi_square_free(dataBlock[0,0], dataBlock[i,0], \
                    stride=numPointsPerChannel, nRounds = numRounds)
        else:
            value = sc.chi_square_free(dataBlock[0,0], dataBlock[i,0], dx1=dataBlock[i,1], dx2=dataBlock[i,1], \
                    stride=numPointsPerChannel, nRounds=numRounds)
        value=np.sqrt(value)
    elif fitMetric == 'V_R':
        value = sc.volatility_ratio(dataBlock[0,0], dataBlock[i,0], stride=numPointsPerChannel )
    elif fitMetric == 'cormap':
        #runs = sc.run_distribution( y1 > f*y2+c )
        value = sc.cormap_value( dataBlock[0,0], dataBlock[i,0] )
    elif fitMetric == 'cormap_matrix':
        mat = sc.cormap_matrix( dataBlock[0,0], dataBlock[i,0] )
        sc.print_cormap_matrix( sys.stdout, mat )
        value=""
    else:
        print >> sys.stderr, "= = = ERROR, metric not recognised! %s" % fitMetric
        sys.exit(1)

    print value

