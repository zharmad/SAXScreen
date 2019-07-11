import general_scripts as gs
import general_maths as gm
import numpy as np
import argparse
#import scipy.integrate as S

def prune_data(min, max, ref, array=[]):
    ind=np.where(np.logical_and(min<ref, ref<max))
    out=[ ref[ind] ]
    for a in range(len(array)):
        out.append( array[a][ind] )
    return out

def run_integration(x, y, dy=[], int_type='1'):
    if int_type=='1':
        yl=y
    elif int_type=='x':
        yl=y*x
    elif int_type=='x^2':
        yl=y*x*x
    val=np.trapz(yl, x)
    if dy==[]:
        return val
    else:
        e1=y-dy ; e2=y+dy
        if int_type=='x':
            e1=e1*x ; e2=e2*x
        elif int_type=='x^2':
            e1=e1*x*x ; e2=e2*x*x
        b1=np.trapz(e1,x) ;  b2=np.trapz(e2,x)
        error=0.5*(b2-b1)
        return (val,error)

def run_Pr_integration(r, Pr, dPr=[], qmin=0.0, qmax=5.0, qpts=251):
    tpi=np.pi
    q=np.linspace(qmin/tpi, qmax/tpi, qpts)
    Iq=np.zeros(qpts)

    if dy==[]:
        for i in range(qpts):
            Iq[i] = 4.0*np.pi*np.trapz( Pr*np.sinc(q[i]*r), r)
        return q, Iq
    else:
        dIq=np.zeros(qpts)
        for i in range(qpts):
            Iq[i] = 4.0*np.pi*np.trapz( Pr*np.sinc(q[i]*r), r)
            d1 = 4.0*np.pi*np.trapz( (Pr-dPr)*np.sinc(q[i]*r), r)
            d2 = 4.0*np.pi*np.trapz( (Pr+dPr)*np.sinc(q[i]*r), r)
            dIq[i] = 0.5*(d2-d1)
        return q*tpi, Iq, dIq

def run_normalised_ratio(y1, y2):
    out = np.divide(y2, y1)
    return out/np.mean(out)

def run_Porod_integration(x, y):
    """
    Run trapozoidal integration over the function y(x)*x^2, returning the cumulative sum.
    The output array length is one less than the input since integration is performed.
    """
    x2y = np.multiply(y,np.power(x,2.0))
    out = np.multiply( 0.5*(x2y[1:]+x2y[:-1]), x[1:]-x[:-1] )
    return x[1:], np.cumsum(out)

def run_transformation(x, y, dy, int_type='1'):
    if int_type=='1':
        yl=y ; dyl = dy
    elif int_type=='x':
        yl=y*x ; dyl = dy*x
    elif int_type=='x^2':
        yl=y*x*x ; dyl= dy*x*x

    return x, yl, dyl

parser = argparse.ArgumentParser(description='Performs actions on an input curve/distribution. '
                                             'This can include regularisation, integration, and transformation. '
                                             'E.g. one can enforce mean of 0, integration volume of 1, sigma of 1, etc. - '
                                             'using the central moments as the definition:\n'
                                             '   mu_n = E[(X-E[X])^n] ',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-f', '--file', type=str, dest='ifile', default='', required=True,
                    help='Input data file.')
parser.add_argument('-g', '--file2', type=str, dest='ifile2', default=None, required=False,
                    help='Second data file for comparative analyses.')
parser.add_argument('--integrate', dest='bInt', action='store_true',
                    help='Report integrated values of the distribution.')
parser.add_argument('--transform', dest='bTrans', action='store_true',
                    help='Return the curve transformed by the function int_type.')
parser.add_argument('--porod', dest='bPorod', action='store_true',
                    help='Return the numerical integration of Int_0^x dx y(x)x^4 ')
parser.add_argument('--normRatio', dest='bNormRatio', action='store_true',
                    help='Return y(2)/y(1) normalised such at its mean is unity.')
parser.add_argument('--xscale', type=float, default=1.0,
                    help='Multiply the transformed curve by this value.')
parser.add_argument('--yscale', type=float, default=1.0,
                    help='Multiply the transformed curve by this value.')
parser.add_argument('--xmin',type=float,default=-np.inf)
parser.add_argument('--xmax',type=float,default=np.inf)
parser.add_argument('--int_type', type=str, dest='int_type', default='1',
                    help='Modify the integral by some prefactor, e.g. "x" for I=S_min^max dx x*y.'
                    'The Pr option returns the SAS intensity based on P(r) integral, rather than a single value.'
                    'Valid values at the moment are: 1, x, x^2, Pr.' )
parser.add_argument('--moments', dest='bMom', action='store_true',
                    help='Report central moments of the distribution.')
parser.add_argument('--symm', dest='bSymm', action='store_true',
                    help='For the moments, assume that the input data is a one-sided function of x>0, '
                         'with -x=x. Thus, always a symmetric distribution.')
parser.add_argument('--err', '--error', dest='bError', action='store_true',
                    help='Assume there is an error column in the input file, and then compute the result with errors.')
parser.add_argument('--qmin',type=float,default=0.0)
parser.add_argument('--qmax',type=float,default=3.0)
parser.add_argument('--qpts',type=float,default=151)

args = parser.parse_args()

if args.xmin == -np.inf and args.xmax == np.inf:
    bPrune=False
else:
    bPrune=True

x=None ; y=None ; dy=None
x2=None ; y2=None ; dy2=None
if args.bError:
    x,y,dy = gs.load_xydy(args.ifile)
    if bPrune:
        [x,y,dy] = prune_data(args.xmin, args.xmax, x, [y,dy] )
    if not args.ifile2 is None:
        x2,y2,dy2 = gs.load_xydy(args.ifile2)
        if bPrune:
            [x2,y2,dy2] = prune_data(args.xmin, args.xmax, x2, [y2,dy2] )
else:
    x, y = gs.load_xy(args.ifile)
    if bPrune:
        [x,y] = prune_data(args.xmin, args.xmax, x, [y] )
    if not args.ifile2 is None:
        x2, y2 = gs.load_xy(args.ifile2)
        if bPrune:
            [x2,y2] = prune_data(args.xmin, args.xmax, x2, [y2] )

# Integrals section
if args.bInt:
    if args.bError:
        if args.int_type != 'Pr':
            vol=run_integration(x, y, dy, int_type=args.int_type)
            print vol[0], vol[1]
        else:
            q, Iq, dIq = run_Pr_integration(x, y, dy, qmin=args.qmin, qmax=args.qmax, qpts=args.qpts)
            print '# I(q) derived from P(r)'
            for i in range(len(q)):
                print q[i], Iq[i], dIq[i]
    else:
        if args.int_type != 'Pr':
            vol=run_integration( x, y, int_type=args.int_type)
            print vol
        else:
            q, Iq = run_Pr_integration(x, y, qmin=args.qmin, qmax=args.qmax, qpts=args.qpts)
            print '# I(q) derived from P(r)'
            for i in range(len(q)):
                print q[i], Iq[i]

if args.bPorod:
    xo, yo = run_Porod_integration(x, y)   
    for i in range(len(xo)):
        print xo[i], yo[i]

if args.bNormRatio:
    yo = run_normalised_ratio(y, y2)
    for i in range(len(yo)):
        print x[i], yo[i]

# Moments of inertia section
if args.bMom:
    mom = gm.calc_central_moments(x, y, args.bSymm)
    print str(mom).strip('[]')

# Transformation section
if args.bTrans:
    if not args.bError:
        dy = np.zeros(y)

    q, Iq, dIq = run_transformation(x, y, dy, int_type=args.int_type)
    print '# Transformed curve'
    for i in range(len(q)):
        print args.xscale*q[i], args.yscale*Iq[i], args.yscale*dIq[i]
