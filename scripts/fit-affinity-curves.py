
import sys, os, math
import argparse
import general_scripts as gs
import numpy as np
from scipy.optimize import fmin_powell
from scipy.optimize import minimize

def load_input_file(infn):
    x, y, dy = gs.load_xydy(infn)
    return np.array(x), np.array(y), np.array(dy)

def find_x_subrange(xinp, xtarg, startFit, endFit):
    i0 = 0
    while (xinp[i0] < xtarg[0] or xinp[i0] < startFit):
        i0 += 1
    i1 = len(xinp)
    while (xinp[i1-1] > xtarg[-1] or xinp[i1-1] > endFit):
        i1 -= 1
    return xinp[i0:i1], i0, i1

# = = = = Debugging
def debug_data_block ( db ):
    print "= = = Debug data_block."
    print "...number of conc series :", len(db)
    for i in range(len(db)):
        print "......series %i has %i values." % (i, len(db[i]))
        print "......first entry:", db[i][0]
    return

def write_xvg_header_simexp(fp, ncurves):
    cmax=31
    s=0
    c=0
    print >> fp, '@with g0'
    print >> fp, '@ view 0.1, 0.1, 0.7, 0.95'
    print >> fp, '@ legend 0.7, 0.95'
    print >> fp, '@ legend char size 0.75'
    print >> fp, '@ legend length 2'
    for i in range(ncurves):
        print >> fp, "@s%i line type 1" % s
        print >> fp, "@s%i line linewidth 1.0" % s
        print >> fp, "@s%i line color %i" % (s, 1+c%cmax )
        print >> fp, "@s%i symbol 0" % s
        print >> fp, "@s%i symbol size 0.25" % s
        print >> fp, "@s%i symbol color %i" % (s, 1+c%cmax  )
        s+=1
        print >> fp, "@s%i line type 0" % s
        print >> fp, "@s%i line linewidth 1.0" % s
        print >> fp, "@s%i line color %i" % (s, 1+c%cmax  )
        print >> fp, "@s%i errorbar color %i" % (s, 1+c%cmax )
        print >> fp, "@s%i errorbar size 0.20" % s
        print >> fp, "@s%i symbol 1" % s
        print >> fp, "@s%i symbol size 0.25" % s
        print >> fp, "@s%i symbol color %i" % (s, 1+c%cmax )
        s+=1 ; c+=1

def write_xvg_legends(fp, prefix, suffix, vals, skip=1):
    s=0
    for i in range(len(vals)):
        print >> fp, "@s%i legend \"%s %5.3g %s\"" % (s, prefix, vals[i], suffix)
        s+=skip

def write_score_file(fn, legs, logKd, header=''):
    fp=open(fn,'w')
    if header != '':
        print >> fp, header
    print >> fp, '@type bardy'
    s=0
    for i in range(len(legs)):
        print >> fp, '@s%i line type 0' % s
        print >> fp, '@s%i symbol fill pattern 6' % s
        print >> fp, '@s%i symbol size 1.0' % s
        print >> fp, '@s%i baseline type 3' % s
        print >> fp, "@s%i errorbar size 0.20" % s
        xval = ''.join([c for c in legs[i] if c in '1234567890.'])
        if xval == '':
            xval=str(s)
        print >> fp, "%s %g %g" % (xval, logKd[i,0], logKd[i,1])
        s+=1
    fp.close()

def write_distrib_file(fn, legs, vlist, header=''):
    fp=open(fn,'w')
    if header != '':
        print >> fp, header
    print >> fp, '@type xy'

    # Determine dimensions.
    shape=vlist.shape
    if len(shape)==1:
        print >> fp, "@s%i legend \"%s\"" % ( 0, legs[0] )
        print >> fp, "# MEAN +- SIG: %g %g" % ( np.mean(vlist), np.std(vlist) )
        for j in range( shape[-1] ):
            print >> fp, "%i %g" % ( j, vlist[j] )
    elif len(shape)==2:
        for i in range( shape[0] ):
            print >> fp, "@s%i legend \"%s\"" % ( i, legs[i] )
            print >> fp, "# MEAN +- SIG: %g %g" % ( np.mean(vlist[i]), np.std(vlist[i]) )
            for j in range( shape[-1] ):
                print >> fp, "%i %g" % ( j, vlist[i,j] )
            print >> fp, '&'
    else:
        print >> sys.stderr, "= = ERROR: Wrong number of dimensions of data!"
        sys.exit(1)
    fp.close()

def write_logKd_histogram(fn, legs, logKds, xmin, xmax):
    nplots=len(legs)
    fp=open(fn,'w')
    for i in range(nplots):
        print >> fp, '@with g%i' % i
        print >> fp, '@s0 legend "%s"' % legs[i]
        print >> fp, '@s0 line type 3'
        h, edges = np.histogram(logKds[:,i], bins=5*int(xmax-xmin), range=(xmin,xmax))
        for j in range(len(h)):
            print >> fp, "%g %g" % ( edges[j], h[j] )
#            print >> fp, "%g %g" % ( 0.5*(np.sum(edges[j:j+2])), h[j] )
        print >> fp, "%g %g" % (edges[-1], 0)
        print >> fp, "&"
    fp.close()

# = = = = Bound ratio as a simple 2-state reaction.
# = = = = [P]+[Q] <-> [PQ]
# [PQ] Kd = ([Pinit] - [PQ])*([Qinit] - [PQ] )
# x*Kd = (A-x)(B-x)
# x/A  = ( (A+B+Kd) - sqrt( (A+B+Kd)^2 - 4AB ) ) / ( 2*A )
# Assume A is the state that is trated as a ratio, and B is neglected.
def boundconc_2state(inP, inQ, Kd):
    s=inP+inQ+math.fabs(Kd)
    PQ=0.5*(s-math.sqrt(s*s-4*inP*inQ))
    return inP-PQ, inQ-PQ, PQ


# This function assumes that the concentration of species P is fixed,
# and the X-data of the rows are in terms of P:Q ratios.
def load_targets_multifile_modelA( filelist, concA ):
    output=[]
    for i in range(len(filelist)):
        x, y = gs.load_xy(efiles[i])
        series = np.zeros( (len(x),3) )
        for j in range(len(x)):
            series[j] = [ concA, concA*x[j], y[j] ]
        output.append( series )
    return output

# datablock format:
# a list of N-titrations. Each with its own set of [inP],[inQ], and titrations
def load_targets_singlefile_modelA( fn, concA, bGetSig=False ):
    output=[]
    errors=[]
    legs, xlist, ylist, dylist = gs.load_sxydylist(fn)
    if bGetSig:
        for i in range(len(xlist)):
            series = np.zeros( (len(xlist[i]),4) )
            for j in range(len(xlist[i])):
                series[j] = [ concA, concA*xlist[i][j], ylist[i][j], dylist[i][j] ]
            output.append( series )
    else:
        for i in range(len(xlist)):
            series = np.zeros( (len(xlist[i]),3) )
            for j in range(len(xlist[i])):
                series[j] = [ concA, concA*xlist[i][j], ylist[i][j]]
            output.append( series )
    return output, legs


def gaussian_noise_datablock(raw_block, scale=1.0):
    out=[]
    npts=np.sum([ len(x) for x in raw_block])
    delta=np.random.normal(0.0,scale,npts)
    i=0
    for titr in raw_block:
        series=np.zeros( (len(titr),3) )
        for j in range(len(titr)):
            series[j] =[ titr[j][0], titr[j][1], titr[j][2]+titr[j][3]*delta[i] ]
            i+=1
        out.append(series)
    return out

def split_datablock( raw_block ):
    o1=[]
    o2=[]
    for i in raw_block:
        p1=[] ; p2=[]
        for j in i:
            p1.append(j[0:3])
            p2.append([ j[0],j[1],j[3] ] )
        o1.append(np.array(p1)) ; o2.append(np.array(p2))
    return o1, o2


# Specific for the format of the datablock above.
def estimate_initial_parameters_modelA_old( datablock ):
    valA=[]
    valB=[]
    param = [] # The other params here are Kds of individial series.
    for i in range(len(datablock)):
        db=datablock[i]
        loc=np.argmin( db[:,1] )
        valA.append( db[loc,2] )
        loc=np.argmax( db[:,1] )
        valB.append( db[loc,2] )
        param.append( 0.8*np.min(db[:,1])+0.2*np.max(db[:,1]) )

    param.insert(0, np.mean(valB) )
    param.insert(0, np.mean(valA) )
    return param

def obtain_PQ_conc_from_datablock( datablock ):
    PQset=set([])
    for titr in datablock:
        for point in titr:
            conc=( point[0], point[1] )
            if conc not in PQset:
                PQset.add( conc )
    return PQset

# Conc being a 2-tuple ( [P], [Q] )
def remove_conc_from_datablock( datablock, conc):
    out=[]
    for titr in datablock:
        o=[ pt for pt in titr if ( pt[0]!=conc[0] or pt[1]!=conc[1] ) ]
        out.append(np.array(o))
    return out

def write_fitted_modelA(outfn, legends, params, db):
    ncurves=len(db)
    delta, bases, Kds = extract_params_modelA(params, ncurves)
    Kds = np.fabs( Kds )
    fp = open(outfn, 'w')
    write_xvg_header_simexp(fp, 8)
    #write_xvg_legends(fp, "K\\sd\\N =", "\\xm\\f{}M", Kds, 2)
    if len(db[0][0])==3:
        bErr=False
    elif len(db[0][0])==4:
        bErr=True
    s=0
    for i in range(len(Kds)):
#        print >> fp, "@s%i legend \"%s: S\\s\\xD\\f{}\\N = %5.3g \\xm\\f{}M\"" % (s, legends[i], Kds[i])
        print >> fp, "@s%i legend \"%s: K\\sD\\N\\S\\eff.\\N = %5.3g \\xm\\f{}M\"" % (s, legends[i], Kds[i])
        s+=2
    for i in range(len(db)):
        print >> fp, "@type xy"
        for j in range(len(db[i])):
            inP = db[i][j][0] ; inQ = db[i][j][1] ; targ=db[i][j][2]
            P, Q, PQ = boundconc_2state( inP, inQ, Kds[i] )
            model_val = (P/inP)*bases + (PQ/inP)*(bases+delta)
            print >> fp, "%g %g" % (inQ/inP, model_val)
        print >> fp, "&"
        if bErr:
            print >> fp, "@type xydy"
        else:
            print >> fp, "@type xy"
        for j in range(len(db[i])):
            inP = db[i][j][0] ; inQ = db[i][j][1] ; targ=db[i][j][2]
            if bErr:
                print >> fp, "%g %g %g" % (inQ/inP, db[i][j][2], db[i][j][3])
            else:
                print >> fp, "%g %g" % (inQ/inP, db[i][j][2])
        print >> fp, "&"
    fp.close()
#    for i in range(len(db)):
#        for j in range(len(db[i])):
#            inP, inQ, targ = db[i][j]
#            P, Q, PQ = boundconc_2state( inP, inQ, Kds[i] )
#            model_val = (P/inP)*valA + (PQ/inP)*valB
#            print >> fp, "%g %g" % (inQ/inP, model_val)
#        print >> fp, "&"
#        for j in range(len(db[i])):
#            inP, inQ, targ = db[i][j]
#            print >> fp, "%g %g" % (inQ/inP, targ)
#        print >> fp, "&"
#    fp.close()

# Fit function for modelA
def fitfunc_modelA(params, *args):
    ncurves=args[0]
    data=args[1]
    delta, bases, Kds = extract_params_modelA(params, ncurves)
    #valA=params[0]
    #valB=params[1]
    #Kds=params[2:]
    #Kds = np.fabs( Kds )
    chi2=0.0
    count=0
    for i in range(len(data)):
        for j in range(len(data[i])):
            inP, inQ, target = data[i][j]
            #print "...minimising...", inP, inQ, target
            P, Q, PQ = boundconc_2state( inP, inQ, Kds[i] )
            model_val = (P/inP)*bases + (PQ/inP)*(bases+delta)
            #model_val = (P/inP)*valA + (PQ/inP)*valB
            chi2  += ( target - model_val )**2.0
            count += 1
    return 1.0*chi2/count

def extract_params_modelA(params, ncurves):
    delta=params[0]
    bases=params[1]
    Kds=params[1:1+ncurves]
    return delta, bases, Kds

def concat_params_modelA(delta, baseline, Kds):
    return [ delta ] + [ baseline ] + list(Kds)

def estimate_initial_parameters_modelA( datablock ):
    mins = []
    deltas = []
    Kds = [] # The other params here are Kds of individial series.
    for i in range(len(datablock)):
        db=datablock[i]
        loc1=np.argmin( db[:,1] )
        setmin = db[loc1,2]
        loc2=np.argmax( db[:,1] )
        setmax = db[loc2,2]
        mins.append( setmin )
        if loc2 > loc1 :
            #print "1:", setmax-setmin
            deltas.append( setmax-setmin )
        else:
            #print "2:", setmin-setmax
            deltas.append( setmin-setmax )
        Kds.append( 0.8*np.min(db[:,1])+0.2*np.max(db[:,1]) )
    params = [ np.mean(deltas) ] + [ np.mean(mins) ] + Kds
    return params

# = = = = = MODEL B: = = = = =
# We permit the fact that each apo has a different starting concentration.

def extract_params_modelB(params, ncurves):
    delta=params[0]
    bases=params[1:1+ncurves]
    Kds=params[1+ncurves:1+2*ncurves]
    return delta, bases, Kds

def concat_params_modelB(delta, baseline, Kds):
    return [ delta ] + list(baseline) + list(Kds)

# Specific for the format of the datablock above.
def estimate_initial_parameters_modelB( datablock ):
    mins = []
    deltas = []
    Kds = [] # The other params here are Kds of individial series.
    for i in range(len(datablock)):
        db=datablock[i]
        loc1=np.argmin( db[:,1] )
        setmin = db[loc1,2]
        loc2=np.argmax( db[:,1] )
        setmax = db[loc2,2]
        mins.append( setmin )
        if loc2 > loc1 :
            #print "1:", setmax-setmin
            deltas.append( setmax-setmin )
        else:
            #print "2:", setmin-setmax
            deltas.append( setmin-setmax )
        Kds.append( 0.8*np.min(db[:,1])+0.2*np.max(db[:,1]) )
    params = [ np.mean(deltas) ] + mins + Kds
    return params

def write_fitted_modelB(outfn, legends, params, db):
    ncurves=len(db)
    delta, bases, Kds = extract_params_modelB(params, ncurves)
    Kds = np.fabs( Kds )
    fp = open(outfn, 'w')
    write_xvg_header_simexp(fp, ncurves)
    if len(db[0][0])==3:
        bErr=False
    elif len(db[0][0])==4:
        bErr=True
    #write_xvg_legends(fp, "K\\sd\\N =", "\\xm\\f{}M", Kds, 2)
    s=0
    for i in range(len(Kds)):
#        print >> fp, "@s%i legend \"%s: S\\s\\xD\\f{}\\N = %5.3g \\xm\\f{}M\"" % (s, legends[i], Kds[i])
        print >> fp, "@s%i legend \"%s: K\\sD\\N\\S\\eff.\\N = %5.3g \\xm\\f{}M\"" % (s, legends[i], Kds[i])
        s+=2
    for i in range(len(db)):
        print >> fp, "@type xy"
        for j in range(len(db[i])):
            inP = db[i][j][0] ; inQ = db[i][j][1] ; targ=db[i][j][2]
            P, Q, PQ = boundconc_2state( inP, inQ, Kds[i] )
            model_val = (P/inP)*bases[i] + (PQ/inP)*(bases[i]+delta)
            print >> fp, "%g %g" % (inQ/inP, model_val)
        print >> fp, "&"
        if bErr:
            print >> fp, "@type xydy"
        else:
            print >> fp, "@type xy"
        for j in range(len(db[i])):
            inP = db[i][j][0] ; inQ = db[i][j][1] ; targ=db[i][j][2]
            if bErr:
                print >> fp, "%g %g %g" % (inQ/inP, db[i][j][2], db[i][j][3])
            else:
                print >> fp, "%g %g" % (inQ/inP, db[i][j][2])
        print >> fp, "&"
    fp.close()

# Fit function for modelB
def fitfunc_modelB(params, *args):
    ncurves=args[0]
    data=args[1]
    delta, bases, Kds = extract_params_modelB(params, ncurves)
    chi2=0.0
    count=0
    for i in range(len(data)):
        for j in range(len(data[i])):
            inP, inQ, target = data[i][j]
            #print "...minimising...", inP, inQ, target
            P, Q, PQ = boundconc_2state( inP, inQ, Kds[i] )
            model_val = (P/inP)*bases[i] + (PQ/inP)*(bases[i]+delta)
            chi2  += ( target - model_val )**2.0
            count += 1
    return 1.0*chi2/count

# = = = = = Model C: as with model B, but with an additional contribution from unbound ligand
def extract_params_modelC(params, ncurves):
    delta_cplx=params[0]
    #delta_lig=params[1]
    bases=params[1:2+ncurves]
    Kds=params[1+ncurves:1+2*ncurves]
    return delta_cplx, bases, Kds
    #delta_lig=params[1:1+ncurves]
    #bases=params[1+ncurves:1+2*ncurves]
    #Kds=params[1+2*ncurves:1+3*ncurves]
    #return delta_cplx, delta_lig, bases, Kds

# Specific for the format of the datablock above.
def estimate_initial_parameters_modelC( datablock, apo_conc ):
    mins = []
    deltas = []
    #deltas_lig = []
    Kds = [] # The other params here are Kds of individial series.
    for i in range(len(datablock)):
        db=datablock[i]
        loc1=np.argmin( db[:,1] )
        setmin = db[loc1,2]
        loc2=np.argmax( db[:,1] )
        setmax = db[loc2,2]
        mins.append( setmin / apo_conc)
        if loc2 > loc1 :
            print "1:", setmax-setmin
            deltas.append( (setmax-setmin)/apo_conc  )
        else:
            print "2:", setmin-setmax
            deltas.append( (setmin-setmax)/apo_conc )
        Kds.append( 0.8*np.min(db[:,1])+0.2*np.max(db[:,1]) )
        #deltas_lig.append( -0.1*deltas[-1] )

    #params = [ np.mean(deltas) ] + [ np.mean(deltas)*-0.1 ] + mins + Kds
    params = [ np.mean(deltas) ] + mins + Kds
    return params

def write_fitted_modelC(outfn, legends, params, db):
    ncurves=len(db)
    delta_cplx, bases, Kds = extract_params_modelC(params, ncurves)
    Kds = np.fabs( Kds )
    fp = open(outfn, 'w')
    write_xvg_header_simexp(fp, ncurves)
    #write_xvg_legends(fp, "K\\sd\\N =", "\\xm\\f{}M", Kds, 2)
    s=0
    for i in range(len(Kds)):
        print >> fp, "@s%i legend \"%s: K\\sd\\N = %5.3g \\xm\\f{}M\"" % (s, legends[i], Kds[i])
        s+=2
    for i in range(len(db)):
        for j in range(len(db[i])):
            inP, inQ, targ = db[i][j]
            P, Q, PQ = boundconc_2state( inP, inQ, Kds[i] )
            Qfrac = -0.01
            model_val = P*bases[i] + PQ*(bases[i]+delta_cplx) + Qfrac*Q*bases[i]
            print >> fp, "%g %g" % (inQ/inP, model_val)
        print >> fp, "&"
        for j in range(len(db[i])):
            inP, inQ, targ = db[i][j]
            print >> fp, "%g %g" % (inQ/inP, targ)
        print >> fp, "&"
    fp.close()

# Fit function for modelC
def fitfunc_modelC(params, *args):
    ncurves=args[0]
    data=args[1]
    delta_cplx, bases, Kds = extract_params_modelC(params, ncurves)
    chi2=0.0
    count=0
    for i in range(len(data)):
        for j in range(len(data[i])):
            inP, inQ, target = data[i][j]
            #print "...minimising...", inP, inQ, target
            P, Q, PQ = boundconc_2state( inP, inQ, Kds[i] )
            #print "...minimising2...", P, Q, PQ
            Qfrac = -0.01
            model_val = P*bases[i] + PQ*(bases[i]+delta_cplx) + Qfrac*Q*bases[i]
            chi2  += ( target - model_val )**2.0
            count += 1
    return 1.0*chi2/count

# = = = = = MODEL D and E: = = = = =
# Instead of permitting the K_D to vary across individual titrations,
# We will permit the maximum to vary across individual titrations.
# This is the case for, e.g. NMR chemical shift perturbation?
#
# In Model D, we will also fix the base value to be zero,
# a single K_D value shared across all titrations, and
# let the value rises from 0 to delta[i] for each titration.
#
# In Model E, we will permit K_Ds to change per titration.
# This decouples all of the titrations from each other.

def csp_ratio_from_conc( inP, inQ, Kd ):
    # From Williamson 2013, which is basically identical to a direct ratio in the limit of fast exchange.. -_-
    # P, Q, PQ = boundconc_2state( inP, inQ, Kd )
    s = inP + inQ + math.fabs(Kd)
    csp_ratio = ( s - math.sqrt(s*s-4*inP*inQ ) )/(2.0*inP)
    return csp_ratio

def extract_params_modelD(params, ncurves):
    deltas=params[0:ncurves]
    Kd=params[ncurves]
    return deltas, Kd

def concat_params_modelD(deltas, Kd):
    return list( deltas ) + [ Kd ]

# Specific for the format of the datablock above.
def estimate_initial_parameters_modelD( datablock ):
    deltas = []
    Kds = [] # The other params here are Kds of individial series.
    for i in range(len(datablock)):
        db=datablock[i]
        loc1=np.argmin( db[:,1] )
        setmin = db[loc1,2]
        loc2=np.argmax( db[:,1] )
        setmax = db[loc2,2]
        if loc2 > loc1 :
            #print "1:", setmax-setmin
            deltas.append( setmax-setmin )
        else:
            #print "2:", setmin-setmax
            deltas.append( setmin-setmax )
        Kds.append( 0.8*np.min(db[:,1])+0.2*np.max(db[:,1]) )
    params = list(deltas) + [ np.mean(Kds) ]
    return params

def write_fitted_modelD(outfn, legends, params, db, conc_rec):
    ncurves=len(db)
    deltas, Kd = extract_params_modelD(params, ncurves)
    Kd = np.fabs( Kd )
    fp = open(outfn, 'w')
    write_xvg_header_simexp(fp, ncurves)
    if len(db[0][0])==3:
        bErr=False
    elif len(db[0][0])==4:
        bErr=True
    print >> fp, "@subtitle \"K\\sD\\N: %5.3g \\xm\\f{}M ; [rec]: %g \\xm\\f{}M \"" % ( Kd, conc_rec )
    s=0
    for i in range(len(deltas)):
        print >> fp, "@s%i legend \"%s\"" % (s, legends[i])
        s+=2
    for i in range(len(db)):
        print >> fp, "@type xy"
        for j in range(len(db[i])):
            inP = db[i][j][0] ; inQ = db[i][j][1] ; targ=db[i][j][2]
            model_val = csp_ratio_from_conc( inP, inQ, Kd ) * deltas[i]
            #P, Q, PQ = boundconc_2state( inP, inQ, Kd )
            #model_val = (PQ/inP)*deltas[i]
            print >> fp, "%g %g" % (inQ/inP, model_val)
        print >> fp, "&"
        if bErr:
            print >> fp, "@type xydy"
        else:
            print >> fp, "@type xy"
        for j in range(len(db[i])):
            inP = db[i][j][0] ; inQ = db[i][j][1] ; targ=db[i][j][2]
            if bErr:
                print >> fp, "%g %g %g" % (inQ/inP, db[i][j][2], db[i][j][3])
            else:
                print >> fp, "%g %g" % (inQ/inP, db[i][j][2])
        print >> fp, "&"
    fp.close()

def fitfunc_modelD(params, *args):
    ncurves=args[0]
    data=args[1]
    deltas, Kd = extract_params_modelD(params, ncurves)
    chi2=0.0
    count=0
    for i in range(len(data)):
        for j in range(len(data[i])):
            inP, inQ, target = data[i][j]
            #print "...minimising...", inP, inQ, target
            model_val = csp_ratio_from_conc( inP, inQ, Kd ) * deltas[i]
            #P, Q, PQ = boundconc_2state( inP, inQ, Kd )
            #model_val = (PQ/inP)*deltas[i]
            chi2  += ( target - model_val )**2.0
            count += 1
    return 1.0*chi2/count

# = = = = =

def modelfunc_2state(y1, y2, popB):
    return (1-popB)*y1+popB*y2

def modelfit_2state(params, *args):
    f   =params[0]
    popB=params[1]

    ay=args[0] ; ady=args[1] #/ay ; ady=ady/np.average(ady)
    by=args[2] ; bdy=args[3] #/by ; bdy=bdy/np.average(bdy)
    ey=args[4] ; edy=args[5] #/ey ; edy=edy/np.average(edy)
    bWeight=args[6]
    bWeight=args[7]

    npts=len(ay)
    chi2=0.0
    if bWeight:
        if bLog:
            for i in range(npts):
                if ey[i]<0:
                    continue
                sig2  = (edy[i]/ey[i])**2 + (1-popB)*(ady[i]/ay[i])**2 + popB*(bdy[i]/by[i])**2
                chi2 += (math.log(ey[i]) - math.log(f*( (1-popB)*ay[i] + popB*by[i] )) )**2 / sig2
        else:
            for i in range(npts):
                sig2  = edy[i]**2 + (1-popB)*ady[i]**2 + popB*bdy[i]**2
                chi2 += (ey[i] - f*( (1-popB)*ay[i] + popB*by[i] ) )**2 / sig2

    else:
        if bLog:
            for i in range(npts):
                if ey[i]<0:
                    continue
                chi2 += (math.log(ey[i]) - math.log(f*( (1-popB)*ay[i] + popB*by[i] )) )**2
        else:
            for i in range(npts):
                chi2 += (ey[i] - f*( (1-popB)*ay[i] + popB*by[i] ) )**2

    #print "= = = In minimisation function: f=%g popB=%g chi2=%g bWeight=%r" % (f, popB, chi2/npts, bWeight)
    return chi2/npts

scriptname=os.path.basename(__file__)

parser = argparse.ArgumentParser(description='Conduct a fit of effective K_D given the concentrations of receptor:ligand'
                                'and some measure that is a proxy of their relative populations.'
                                'This script operates in two modes:'
                                'In mode 1, the measure for each rec:lig is a full curve, with the assumption that the complex curve is unknown.'
                                'The script assumes that, for each *different* ligand titrated, the change in the receptor apo-curve is identical,'
                                'and thus can be derived from their respective apo curves.'
                                'In mode two, the measure for each rec:lig is a single value, agains with the asusmption that the complexation'
                                'produces an identical change upon the apo curve.'
                                'For both modes, an apo-receptor measure is expected at the beginning of each ligand titration series.'
                                'All data will then be fitted together.',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--curves_complex', type=str, dest='rfile', default='', required=False,
                    help='Mode 1: A file containing three columns: absolute concentrations of the receptor and ligand, then the file name of the curve.'
                    'For example, this is the list of measured combined SAXS curves in a titration. The apo receptor curve should be contained here as the first row.')
parser.add_argument('--curves_lig', type=str, dest='cfile', default='', required=False,
                    help='Mode 1: Optional file that contains two columns: absolute concentration of the ligand, then the file name of the curve.'
                    'This triggers a mode where the unbound ligand contributes to the final measure additively, the unbound ligand signal is otherwise ignored.'
                    'That is I_tot = (c_bound) I_complex + (c_apo-c_bound) I_apo + (c_lig-c_bound) I_lig.'
                    'For example, this is the list of pure ligand SAXS curves in a titration.')
parser.add_argument('-t', '--targ', type=str, dest='tfile', default='', required=False,
                    help='Mode 2: Single file containing one measurement per receptor:ligand combo, in \'&\'-delimited groups.'
                    'For example, this is the Kurtosis of the experimental SAXs curve.'
                    'The X-values represent lig:rec ratios, e.g. 0.0 <-> 1.0')
parser.add_argument('--conc_rec', type=float, dest='concRec', default='50',
                    help='Mode 2: The concentration of the receptor required to interpret the X-column in --targ.')
#parser.add_argument('efiles', metavar='Exp_files', type=str, nargs='+',
#                    help='Files of all experimental intensity curves.')
parser.add_argument('-o', '--outpref', type=str, dest='opref', default='', required=True,
                    help='Output prefix, for recording data such as estimated affinity, model fits, histograms, etc.')
parser.add_argument('--bLog', dest='bLog', action='store_true',
                    help='Use the logarithm of the fitting function.')
parser.add_argument('--bSigma', dest='bSigma', action='store_true',
                    help='Use the third column, uncertainty, as a weight.')
parser.add_argument('--bNMRTitration', action='store_true',
                    help='Use the model for NMR titration, where one assumes all CSPs represent a single Kd curve '
                        'with different deltas and a minimum of zero.')
parser.add_argument('--3body', dest='b3Body', action='store_true',
                    help='Also model a contribution from unbound ligand. Initial estimate will be -10%% of the holo-change.')
parser.add_argument('--errormode', type=str, dest='err_type', default='none',
                    help='Define a model of error-analysis. Options are: none, 1-point, noise2sig, noise1sig.'
                    '1-point means to remove data at one ratio at a time and analyse the difference.'
                    'noise2sig means to add random noise to the mean value, taking the errors of'
                    'said value to be 2-sigma, the 95%% confidence interval.')
parser.add_argument('--ntrials', type=int, dest='ntrials', default=50,
                    help='In Error analyses that require trials, run this many trials. Applicable to, e.g. noise1sig.')
parser.add_argument('--bSmallPerturb', action='store_true',
                    help='Perform an self-correction routine where the expected change in signal should not be greater than the baseline. '
                         'This is outdated and kept for historical purposes.' )

args = parser.parse_args()

opref=args.opref
fileScore1=args.opref+'_score1.xvg'
fileScore2=args.opref+'_score2.xvg'
fileModel=args.opref+'_model.xvg'
fileHist=args.opref+'_hist.xvg'

bLog=args.bLog
bWeight=args.bSigma
bLig=False
bSmallPerturb=args.bSmallPerturb

errMode=args.err_type
if errMode=='noise1sig' or errMode=='noise2sig':
    bReadSig=True
else:
    bReadSig=False

#Determine which mode the python script is being called
mode=0
#Mode 1
if args.rfile != '':
    mode=1
    rfiles = args.rfile
    if args.cfile != '':
        bLig=True
        cfiles=args.cfile
#Mode 2
if args.tfile != '':
    if mode != 0:
        print >> sys.stderr, "= = ERROR: Cannot use both mode 1 and mode 2. Please check your arguments."
        sys.exit(1)
    mode=2
    efiles = args.tfile
    concRec= args.concRec
# Sanity check
if mode == 0:
    print >> sys.stderr, "= = ERROR: The script must operate in either mode 1 or mode 2. Please check help  with -h."

if mode==2:
    # Standard mode for taking a single measure per titration point.
    # = = =
    # Pose data array as 3D entries:
    # X-values in our target are the ratios of receptor ligand.
    data_block, legends = load_targets_singlefile_modelA( efiles, concRec, bReadSig )
    #data_block, err_block = split_datablock( data_block )
    debug_data_block( data_block )

    if args.b3Body:
        # Adopt a model that sets a different baseline for each titration series,
        # but the change in signal upon binding is conserved, as well as the unbound ligand signal.
        # Reminder: Parameters are, 2 deltas, N minimums, and N K_D for each concentration series in datablock.
        # i.e.: 2N+2 parameters to be fillted.
        params = estimate_initial_parameters_modelC( data_block, concRec )
        ncurves = len(data_block)
        print "= = = Initial parameter estimates:"
        print params
        dummy=0
        fminOut = fmin_powell(fitfunc_modelC, params, args=(ncurves, data_block),
                              full_output=True)
        xopt    = fminOut[0]
        funcopt = fminOut[1]
        print "= = Optimised parameters: "
        print xopt
        write_fitted_modelC( fileModel, legends, xopt, data_block )
        print "= = = Written model and targets to %s ." % fileModel
    elif args.bNMRTitration:
    	# Adopt NMR titration mode using chemical shift perturbation.
	# Which corresponds to model-D in this script.
        params = estimate_initial_parameters_modelD( data_block )
        ncurves = len(data_block)
        print "= = = Initial parameter estimates:"
        print params
        fminOut = fmin_powell(fitfunc_modelD, params, args=(ncurves, data_block),
                              full_output=True)
        xopt    = fminOut[0]
        funcopt = fminOut[1]
        print "= = Optimised parameters: "
        print xopt
        write_fitted_modelD( fileModel, legends, xopt, data_block, concRec )
        print "= = = Written model and targets to %s ." % fileModel

    elif True:
        # Some hacking here to switch between model A and B
        # Adopt a model where, for each titration series, the apo-baseline can and will change,
        # but the magnitude of the change is conserved
        # Reminder: Parameters are, 1 shared delta, 1 shared minima, and N K_D for each concentration series in datablock.
        # This is model A
        # i.e.: N_series + 2 parameters
        # Model B gives N minima instead of just one.
        if errMode=='none':
            params = estimate_initial_parameters_modelB( data_block )
            ncurves = len(data_block)
            print "= = = Initial parameter estimates:"
            print params
            fminOut = fmin_powell(fitfunc_modelA, params, args=(ncurves, data_block), full_output=True)
            #fminOut = fmin_powell(fitfunc_modelB, params, args=(ncurves, data_block), full_output=True)
            xopt    = fminOut[0] ; funcopt = fminOut[1]
            print "= = Optimised parameters: "
            print xopt
            write_fitted_modelA( fileModel, legends, xopt, data_block )
            #write_fitted_modelB( fileModel, legends, xopt, data_block )
            print "= = = Written model and targets to %s ." % fileModel
        elif errMode=='1-point':
            params = estimate_initial_parameters_modelB( data_block )
            ncurves = len(data_block)
            conc_list=obtain_PQ_conc_from_datablock( data_block )
            nconcs = len(conc_list)
            print "= = List of concentrations detected:"
            print conc_list
            logKd=[]
            for c in conc_list:
                db=remove_conc_from_datablock(data_block, c)
                fminOut = fmin_powell(fitfunc_modelB, params, args=(ncurves, db),
                                  full_output=True)
                xopt    = fminOut[0] ; funcopt = fminOut[1]
                delta, bases, Kds = extract_params_modelB(xopt, ncurves)
                tmp =  np.clip(np.log10(np.fabs(Kds)),-2.0,5.0) ; print "= = log10-Kd:", tmp
                logKd.append( tmp )
            logKd=np.array(logKd)
            mean = np.mean(logKd,axis=0)
            std  = np.std(logKd,axis=0)
            write_score_file( fileScore1, legends, np.stack( (mean, std), axis=-1 ))
            Kd_mean = np.mean(np.power(10.0,logKd),axis=0)
            Kd_std  = np.std(np.power(10.0,logKd),axis=0)
            write_score_file( fileScore2, legends, np.stack( (Kd_mean, Kd_std), axis=-1 ))
        elif errMode=='noise2sig' or errMode=='noise1sig':
            print "= = Error mode requested: noise-addition "
            if errMode=='noise2sig':
                efact=0.5
            else:
                efact=1.0
            ncurves=len(data_block)
            ntrials=args.ntrials
            deltas=np.zeros(ntrials)
            baselines=np.zeros( ntrials )
            #baselines=np.zeros( (ntrials,ncurves) )
            logKd=np.zeros( (ntrials,ncurves) )
            i=0 ; nRedo=0 ; nRedoTot=0
            while True:
            #for i in range(ntrials):
                print "= = = Conducting trial %i of %i ..." % (i+1, ntrials)
                db=gaussian_noise_datablock(data_block, scale=efact)
                #print "Debug:", db
                params = estimate_initial_parameters_modelA( db )
                #print "Debug:", params
                fminOut = fmin_powell(fitfunc_modelA, params, args=(ncurves, db), full_output=True)
                #params = estimate_initial_parameters_modelB( db )
                #fminOut = fmin_powell(fitfunc_modelB, params, args=(ncurves, db), full_output=True)
                xopt    = fminOut[0] ; funcopt = fminOut[1]

                deltas[i], baselines[i], Kds = extract_params_modelA(xopt, ncurves)
                #deltas[i], baselines[i], Kds = extract_params_modelB(xopt, ncurves)

                #if True:
                #    # Debug to temporary file.
                #    write_fitted_modelB( 'temp-Kds-%i.dat' % i, legends, xopt, db )
                print "= = Optimised delta:", deltas[i]
                print "= = Optimised Baselines:", baselines[i]
                logKd[i] = np.clip(np.log10(np.fabs(Kds)),-2.0,5.0) ; print "= = log10-Kd:", logKd[i]
                # = = = Note: I think this is an outdated check
                if bSmallPerturb and np.any( np.fabs(deltas[i]) > np.fabs(baselines[i]) ):
                    nRedo+=1 ; nRedoTot+=1
                    # Delta is really out of whack. Redo.
                    if nRedo<5:
                        print "= = = ERROR: optimised delta is greater than twice the baseline! This is wrong. Redoing... (Attempts: %i)" % nRedo
                        continue
                    else:
                        print "= = = ABORTING ERROR: Too many redo attempts, bailing! Please check the delta limits in the script and your model."
                        sys.exit(2)
                else:
                    i+=1 ; nRedo=0
                    if i==ntrials:
                        break

            if nRedoTot > 0:
                print "= = = NOTE: A total of %i trials have been redone due to tripping the error signal." % nRedoTot
            logKd_mean = np.mean(logKd,axis=0)
            logKd_std  = np.std(logKd,axis=0)
            Kd_mean = np.mean(np.power(10.0,logKd),axis=0)
            Kd_std  = np.std(np.power(10.0,logKd),axis=0)
            write_distrib_file( opref+'_logKd.xvg',    legends, logKd.T,    header='# Fitted Log_10(Kd)\n# Number of trials: %i' % ntrials)
            write_distrib_file( opref+'_baseline.xvg', legends, baselines.T, header='# Fitted baselines\n# Number of trials: %i' % ntrials)
            write_distrib_file( opref+'_Deltas.xvg', ['Shared-Delta'], deltas, header='# Fitted Deltas\n# Number of trials: %i' % ntrials)
            #write_score_file( fileScore1, legends, np.stack( (logKd_mean, logKd_std), axis=-1 ), header='# Number of trials: %i' % ntrials )
            #write_score_file( fileScore2, legends, np.stack( (Kd_mean, Kd_std), axis=-1 ), header='# Number of trials: %i' % ntrials )
            write_logKd_histogram(fileHist, legends, logKd, xmin=-2.0, xmax=5.0)

            xopt_avg = concat_params_modelA( np.mean(deltas), np.mean(baselines), np.power(10, logKd_mean) )
            write_fitted_modelA(fileModel, legends, xopt, data_block )
            #xopt_avg = concat_params_modelB( np.mean(deltas), np.mean(baselines, axis=0), np.power(10, logKd_mean) )
            #write_fitted_modelB(fileModel, legends, xopt_avg, data_block )

    else:
        # Adopt a simpler model where all apo-baselines are assumed to be the same.
        # So variations given are just errors.

          # Reminder: Parameters are, valA, valB, and 1 K_D for each concentration series in datablock.
        # i.e.: N_series + 2 parameters
        params = estimate_initial_parameters_modelA( data_block )
        print "= = = Initial parameter estimates:"
        print params
        dummy=0
        fminOut = fmin_powell(fitfunc_modelA, params, args=(data_block, dummy),
                              full_output=True)
        xopt    = fminOut[0]
        funcopt = fminOut[1]
        #xopt = xopt*np.sign(xopt)
        print "= = Optimised parameters: "
        print xopt

        write_fitted_modelA( fileModel, xopt, data_block )
        print "= = = Written model and targets to %s ." % fileModel

if mode==1:
    #Standard mode for taking a full curve for each measurement
    data_block, legends = load_targets_multifile_threecol( efiles )
    debug_data_block( data_block )
