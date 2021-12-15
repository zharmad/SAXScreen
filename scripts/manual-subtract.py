
import sys, os, math
import argparse
import general_scripts as gs
import numpy as np
from scipy.optimize import fmin_powell

def guinier_fit_with_buffer_subtraction(params, *args):
    """
    Conduct a Guinier fit of the SAXS curve. Optimise the following straight line function:
    log(I) = log(I0) - q^2.Rg^2/3 ,
    
    Three parameters are optimised:
    - f2, the scale of the buffer.
    - logI0, the extrapolated forward scattering
    - Rg, the fitted Radius of Gyration

    Using only points below q.Rg < cutoff, which
    - should be set to 1.3 for gobular proteims
    - or 0.8 for elongated proteins.
    """
    q  = args[0]
    I1 = args[1]
    I2 = args[2]
    #I1 = args[1][0] ; dI1 = args[1][1]
    #I2 = args[2][0] ; dI2 = args[2][1]
    f1 = args[3]
    cutoff = args[4]
    f2 = params[0]
    logI0 = params[1]
    Rg = params[2]
    nPoint = len(q)
    # = = = Simple version for now. Parallelise later.
    chiSum = 0.0
    for i in range(nPoint):
        if q[i]*Rg >= cutoff:
            break
        chi = np.power( np.log( f1*I1[i] - f2*I2[i] ) - logI0 + q[i]*q[i]*Rg*Rg /3.0 , 2 )
        chiSum += chi

    return np.sqrt(chiSum)/i
    #ody = f3*( np.power( f1*f1*np.power(sdy,2) + f2*f2*np.power(bdy_av,2) , 0.5 ) )

scriptname=os.path.basename(__file__)

parser = argparse.ArgumentParser(description='Does a manual subtraction using buffer/sample/buffer scheme.',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-s', '--sample', type=str, dest='sfile', default='', required=True,
                    help='Name of the sample intensity file.')
parser.add_argument('buffers', metavar='Buffer_file', type=str, nargs='+',
                    help='Files of all relevant buffers. Please use at least 2 when taking from the raw synchrotron measurements.')
parser.add_argument('-o', '--out', type=str, dest='ofile', default='', required=True,
                    help='Name of the output intensity file.')
parser.add_argument('-f1', '--scale_sample', type=float, dest='scale1', default=1.0,
                    help='multipication factor of the *sample* intensity when subtracting. This should rarely be used.')
parser.add_argument('-f2', '--scale_buffer', type=float, dest='scale2', default=1.0,
                    help='multipication factor of the averaged *buffer* intensity when subtracting. This should rarely be used.')
parser.add_argument('-f3', '--scale_after', type=float, dest='scale3', default='1.0',
                    help='multipication factor after subtraction has been carried out. This should rarely be used.')
parser.add_argument('--add_file',type=str, dest='addfile', default='',
                    help='Add this SAXS curve to the final result. Useful for manual correction of over/undersubtraction.')
parser.add_argument('--add_mult',type=float, dest='addfact', default=1.0,
                    help='Multiply Added curve by this factor before combining with subtracted data.')
parser.add_argument('--fit_Guinier',action='store_true',dest='bFitGuinier',
                    help='Search for improved buffer subtraction by demanding that the Guinier-fit be optimised. '
                         'The parameter -f2 will be rescaled, along with the Rg of the molecule'
                         'This ideally removes inflections at zero-angle.')

args = parser.parse_args()

f1 = args.scale1
f2 = args.scale2
f3 = args.scale3
sx, sy, sdy = gs.load_xydy(args.sfile)
sx = np.array(sx)
sy = np.array(sy)
sdy= np.array(sdy)
#bx1, by1, bdy1 = gs.load_xydy(args.bfile1)
#bx2, by2, bdy2 = gs.load_xydy(args.bfile2)
#load all buffer files.

if args.addfile != '' and args.addfact != 0 :
    bMod=True
    ax, ay, ady = gs.load_xydy(args.addfile)
    if ax[0] != sx[0]:
        print( "= = = ERROR: The first X-value between sample and additonal-correction curves are NOT the same! Aborting." )
        sys.exit()
    ay  *= args.addfact
    ady *= args.addfact
else:
    bMod=False


nq=len(sx)
bfiles=args.buffers
nbuf=len(bfiles)
bx=np.zeros((nbuf, nq))
by=np.zeros((nbuf, nq))
bdy=np.zeros((nbuf, nq))
for i in range(nbuf):
    bx[i], by[i], bdy[i] = gs.load_xydy(bfiles[i])
    if bx[i,0] != sx[0]:
        print( "= = = ERROR: The first X-value between buffer and sample curves are NOT the same! Aborting." )
        sys.exit()

#bx_av=np.average(bx,0)
if nbuf > 1:
    by_av=np.average(by,0)
    bdy_av=np.average(bdy,0)
    print( by.shape, by_av.shape )
else:
    by_av=by[0]
    bdy_av=bdy[0]

# = = Look at fitting options now.

if args.bFitGuinier:
    # = = = Functionality limitation?
    #if bMod:
    #    print( "= = ERROR: Guinier optimisation with buffer modulation not implemented.", file=sys.stderr )
    #    sys.exit(1)

    # Estimate initial parameters. Use the first N points, and then 2N-points after that.
    # Assume that this is an experimental curve with many q-points.
    # Set the nPoints to 10
    estPoints=10
    oy = f1*sy - f2*by_av
    logI0Init = np.log( np.average( oy[2*estPoints:3*estPoints] ) )
    estQ2 = np.average( sx[3*estPoints:4*estPoints] )
    logEstI2 = np.log( np.average( oy[3*estPoints:4*estPoints] ) )
    RgInit = 0.5* np.sqrt( 3.0*(logI0Init - logEstI2)/(estQ2*estQ2) )

    params = [ f2, logI0Init, RgInit ]
    print( "= = = Initial paramter estimates:", params )
    constArgs=( sx, sy, by_av, f1, 1.3)
    # = = = Conduct grid search.
    minChi = 1e99
    for i in np.logspace(-0.3,0.3,21)*f2:
        for j in np.linspace(-1,1,11)+logI0Init:
            for k in np.logspace(-0.3,0.3,11)*RgInit:
                val = guinier_fit_with_buffer_subtraction([i,j,k], *constArgs)
                if val < minChi:
                    params = [i,j,k]
                    minChi = val
                #print( i,j,k,val )
                #print( i, j, k, guinier_fit_with_buffer_subtraction([i,j,k], *constArgs) )
    print( "= = = Refined after grid-search:", params )
    #sys.exit(0)
    print( "= = = .. begin Powell minimization" )

    fminOut = fmin_powell(guinier_fit_with_buffer_subtraction,
                params, args=constArgs,
                full_output=True)

    optParams = fminOut[0] ; funcopt = fminOut[1]

    header="#Fitted chi: %g\n# Fitted I0: %g\n:# Fitted Rg: %g\n" \
            % ( funcopt, np.exp(optParams[1]), optParams[2] )
    f2 = optParams[0]

    oy  = f3*( f1*sy - f2*by_av )
    ody = f3*( np.power( f1*f1*np.power(sdy,2) + f2*f2*np.power(bdy_av,2) , 0.5 ) )

    if bMod:
        oy  = oy + ay
        ody = np.sqrt( np.square(ody) + np.square(ady) )

    gs.print_xydy(args.ofile, sx, oy, ody, header)

else:
    # = = = Run a simple buffer averging and subtraction scheme.
    oy  = f3*( f1*sy - f2*by_av )
    ody = f3*( np.power( f1*f1*np.power(sdy,2) + f2*f2*np.power(bdy_av,2) , 0.5 ) )
    if bMod:
        oy  = oy + ay
        ody = np.sqrt( np.square(ody) + np.square(ady) )

    gs.print_xydy(args.ofile, sx, oy, ody)

#oy=np.zeros(nq)
#ody=np.zeros(nq)
#for i in range(nq):
#    oy[i]  = sy[i] - f*0.5*( by1[i] + by2[i] )
#    ody[i] = math.sqrt( sdy[i]**2.0 + f*f*0.25*(bdy1[i]+bdy2[i])**2.0 )

