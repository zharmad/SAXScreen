import sys, os, math
import argparse
import general_scripts as gs
import numpy as np
from scipy import optimize
#from scipy.optimize import fmin_powell
#from scipy.optimize import curve_fit

def obtain_linear_model(params, datablock):
    return np.sum( params[...,None]*datablock[1], axis=0)

# block is of shape( x/y/dy, data)
def gaussian_noise_curve(block, scale=1.0):
    sh     = block.shape
    nPts   = sh[-1]
    delta  = np.random.normal( 0.0, scale, nPts )
    out    = np.copy(block)
    out[1] = block[1]+block[2]*delta
    return out

# raw_block is of (x/y/dy, files, data)
def gaussian_noise_datablock(db, scale=1.0):
    sh      = db.shape
    out     = np.zeros( db.shape )
    nCurves = sh[1]
    for i in range(nCurves):
        out[:,i,:] = gaussian_noise_curve( db[:,i,:], scale)
    return out

def chi_linear_fit(params, *args):
    # = = Use the slower nanmean at the end to ignore any random nans that pop out of L-BFGS !!
    n_par = len(params)
    targ      = args[0]
    #datablock = args[1]
    scaled_y = params[...,None]*args[1][1]
    scaled_dy = params[...,None]*args[1][2]
    chi_raw = np.power( np.sum(scaled_y, axis=0) - targ[1],  2 )
    weights = np.power( np.sum( np.power(scaled_dy, 2) + np.power(targ[2], 2), axis=0 ), 0.5 )
    chisq = np.nanmean( chi_raw / weights )
    #print "= = ...optimization step over %i variables, chi%s: %g" % (n_par, params, chisq)
    return chisq

def chi_linear_fit_grad(params, *args):
    n_par = len(params)
    targ      = args[0]
    #datablock = args[1]
    scaled_y = params[...,None]*args[1][1]
    scaled_dy = params[...,None]*args[1][2]
    shared_fact = 2.0* ( np.sum(scaled_y, axis=0) - targ[1] )
    weights = np.power( np.sum( np.power(scaled_dy, 2) + np.power(targ[2], 2), axis=0 ), 0.5 )
    grads = [ np.average( ( args[1][1][i] * shared_fact )  / weights) for i in range( n_par) ]
    #print "= = ...optimization step over %i variables, grad%s: %s" % (n_par, params, grad)
    return grad

def load_input_file(infn, min=np.nan, max=np.nan):
    x, y, dy = gs.load_xydy(infn)
    if np.isnan(min):
        if np.isnan(max):
            return np.array(x), np.array(y), np.array(dy)
        else:
            ind = np.where( x<max )[0]
    elif np.isnan(max):
            ind = np.where( x>min )[0]
    else:
        ind = np.where( np.logical_and( x>min, x<max ) )[0]
    return np.array(x[ind]), np.array(y[ind]), np.array(dy[ind])

def load_component_files(fnlist, min, max):
    print "= = ...Loading list of component curves..."
    n_exp = len( fnlist )
    x, y, dy = load_input_file( fnlist[0], min, max )
    n_vals = len(x)
    ox  = np.zeros( (n_exp, n_vals) )
    oy  = np.zeros( (n_exp, n_vals) )
    ody = np.zeros( (n_exp, n_vals) )
    ox[0] = x ; oy[0] = y ; ody[0] = dy

    for i in np.arange(1,n_exp):
        x, y, dy = load_input_file( fnlist[i], min, max )
        if not np.array_equal( ox[0], x ):
            print >> sys.stderr, "= = ERROR encountered: the set of sources do not have the same x-values!"
            print >> sys.stderr, "    ...the failed comparison occurred between curves %i and %i,  with length %i and %i." % ( 0, i, n_vals, len(x) )
            sys.exit(1)
        ox[i] = x ; oy[i] = y ; ody[i] = dy

    for i in range(n_exp):
        print "= = Debug: x,y,dy (%s): %f,%f,%f" % (fnlist[i], ox[i,20], oy[i,20], ody[i,20])
    print "= = ...Loading complete for %i components." % n_exp
    return np.stack( (ox, oy, ody), axis=0 )

def print_model_fit(fn, x, ty, tdy, my, header=''):
    fp = open(fn, 'w')
    print >> fp, header
    print >> fp, '@type xydy'
    npts=len(x)
    for i in range(npts):
        print >> fp, x[i], ty[i], tdy[i]
    print >> fp, '&'
    print >> fp, '@type xy'
    for i in range(npts):
        print >> fp, x[i], my[i]
    print >> fp, '&'
    fp.close()

def print_chi_file(fn, chi, header=''):
    fp = open(fn, 'w')
    print >> fp, header
    n_pts = len(chi)
    for i in range(n_pts):
        print >> fp, "%i %g" % (i, chi[i])
    fp.close()

scriptname=os.path.basename(__file__)

parser = argparse.ArgumentParser(description='Calculate the extent to which a curve can be described by '
                                'a linear sum of another set of curves over the same x-coordinates, where the '
                                'fitted amplitude for each are forced to be non-negative.',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-f', '-t', type=str, dest='tfile', default='', required=True,
                    help='Target composite scattering curve to be linearly modelled by its constituents.')
parser.add_argument('-x', type=str, dest='bfile', default='', required=False,
                    help='Empty-cell curve to be included in modeling. This is permitted to be negative.')
parser.add_argument('cfiles', metavar='cfiles', type=str, nargs='+',
                    help='Files that comprise the expected linear components.')
parser.add_argument('-o', '--outpref', type=str, dest='opref', default='output', required=False,
                    help='Output prefix, for recording data such as linearity, residual difference curves, etc.')
parser.add_argument('-a', '--xmin', type=float, dest='xmin', default=np.nan,
                    help='Minimum x-value to consider.')
parser.add_argument('-b', '--xmax', type=float, dest='xmax', default=np.nan,
                    help='Maximum x-value to consider.')
parser.add_argument('--positive', dest='bPos', action='store_true',
                    help='Restrict all component factors to be positive.')
parser.add_argument('--doError', dest='bError', action='store_true',
                    help='Randomly alter the y-values by dy-values in all input curve(s), and repeat the chi fitting N-times.')
parser.add_argument('--ntrials', type=int, dest='nk', default=100,
                    help='Number of trials of random errors to input.')

args = parser.parse_args()
model_file=args.opref+'_model.dat'
chi_file  =args.opref+'_chis.dat'

bError   = args.bError
n_trials = args.nk

# Reading input files
tx, ty, tdy = load_input_file( args.tfile, min=args.xmin, max=args.xmax )
targ = np.stack( (tx,ty,tdy), axis=0 )
print "= = Debug: x,y,dy (%s): %f,%f,%f" % (args.tfile, tx[20], ty[20], tdy[20])

bPositive = args.bPos
# Check for blank.
bBlank = False
if args.bfile != '':
    bBlank = True
    args.cfiles.insert( 0, args.bfile)

n_comps = len( args.cfiles )
# Datablock has dimensions( x/y/dy, files, data)
datablock = load_component_files( args.cfiles, min=args.xmin, max=args.xmax )

# Leave option to interpolate data.

if not np.array_equal( datablock[0,0], tx ):
    print >> sys.stderr, "= = ERROR encountered: the sources and target do not have the same x-values!"
    print >> sys.stderr, "    ...the lengths are %i and %i." % ( len(datablock[0,0]), len(tx) )
    for i in range(len(tx)):
        if datablock[0,0,i] !=  tx[i]:
            print >> sys.stderr, datablock[0,0,i], "!=", tx[i]
    sys.exit(1)
else:
    print "= = ...X-value check complete, all points at same x."

# Conduct linear fit. Assume after this that the x values are now identical.
shape = datablock.shape
n_comps = shape[1]
n_vals  = shape[2]

if bBlank:
    guess  = np.repeat( 1.0/(n_comps-1), n_comps-1)
    guess  = np.insert(guess, 0, 0.0)
    if bPositive:
        bounds = np.tile( (0.0, None), (n_comps-1,1) )
    else:
        bounds = np.tile( (None, None), (n_comps-1,1) )
    bounds = np.insert(bounds, 0, (None, None), axis=0)
else:
    guess  = np.repeat(1.0/n_comps, n_comps)
    if bPositive:
        bounds = np.tile( (0.0, None), (n_comps,1) )
    else:
        bounds = np.tile( (None, None), (n_comps,1) )

print "= = ...fitting started, using L-BFGS-B..."
if not bError:
    func_opt = optimize.fmin_l_bfgs_b( chi_linear_fit, x0=guess, bounds=bounds, args=(targ, datablock), approx_grad=True )
    #func_opt = optimize.fmin_powell( chi_linear_fit, x0=guess, args=(targ, datablock) )
else:
    opt_pars  = np.zeros( (n_trials, n_comps) )
    chi_block = np.zeros( n_trials )
    for i in range(n_trials):
        targ_this = gaussian_noise_curve( targ, scale=1.0 )
        # print "= = Debug: y= %f+-%f y= %f+-%f " % ( targ[1,20], targ[2,20], targ_this[1,20], targ_this[2,20] )
        func_opt = optimize.fmin_l_bfgs_b( chi_linear_fit, x0=guess, bounds=bounds, args=(targ_this, datablock), approx_grad=True )
        if func_opt[2]['warnflag']>0:
            print >> sys.stderr, "= = WARNING: error encountered in L-BFGS-B minimization, code:", func_opt[2]['warnflag']
            print >> sys,stderr, "   Grad/#Calls/#Iterations:", func_opt[2]['grad'], func_opt[2]['funcalls'], func_opt[2]['nits']
            sys.exit(1)
        else:
            opt_pars[i]  = func_opt[0]
            chi_block[i] = np.sqrt( func_opt[1] )

#opt_pars=func_opt[0]
#opt_chi=np.sqrt(func_opt[1])
chi_mean = np.mean( chi_block )
chi_std  = np.std(  chi_block )
print "= = ...fitting complete. Chi: %g +- %g" % ( chi_mean, chi_std )
pars_mean = np.mean( opt_pars, axis=0 )
pars_std  = np.std(  opt_pars, axis=0 )

model_y = obtain_linear_model( pars_mean, datablock )

header_lines = '# Fitted chi: %g +- %g\n' % ( chi_mean, chi_std )
for i in range( n_comps ):
    header_lines = header_lines + '# Scaling factor for file %i: %g +- %g\n' % ( i+1, pars_mean[i], pars_std[i] )

print "= = ...finishing up and printing to %s and %s" % ( model_file, chi_file )
print_model_fit( model_file, tx, ty, tdy, model_y, header=header_lines )
print_chi_file( chi_file, chi_block, header=header_lines )
