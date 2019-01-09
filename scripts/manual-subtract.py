
import sys, os, math
import argparse
import general_scripts as gs
import numpy as np

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
                    help='multipication factor after subtraction has been caried out. This should rarely be used.')
parser.add_argument('--add_file',type=str, dest='addfile', default='',
                    help='Add this SAXS curve to the final result. Useful for manual correction of over/undersubtraction.')
parser.add_argument('--add_mult',type=float, dest='addfact', default=1.0,
                    help='Multiply Added curve by this factor before combining with subtracted data.')

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
        print "= = = ERROR: The first X-value between sample and additonal-correction curves are NOT the same! Aborting."
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
        print "= = = ERROR: The first X-value between buffer and sample curves are NOT the same! Aborting."
        sys.exit()

#bx_av=np.average(bx,0)
if nbuf > 1:
    by_av=np.average(by,0)
    bdy_av=np.average(bdy,0)
    print by.shape, by_av.shape
else:
    by_av=by[0]
    bdy_av=bdy[0]

#Average buffer and subtract
oy  = f3*( f1*sy - f2*by_av )
ody = f3*( np.power( f1*f1*np.power(sdy,2) + f2*f2*np.power(bdy_av,2) , 0.5 ) )
if bMod:
    oy  = oy + ay
    ody = np.sqrt( np.square(ody) + np.square(ady) )

#oy=np.zeros(nq)
#ody=np.zeros(nq)
#for i in range(nq):
#    oy[i]  = sy[i] - f*0.5*( by1[i] + by2[i] )
#    ody[i] = math.sqrt( sdy[i]**2.0 + f*f*0.25*(bdy1[i]+bdy2[i])**2.0 )

gs.print_xydy(args.ofile, sx, oy, ody)
