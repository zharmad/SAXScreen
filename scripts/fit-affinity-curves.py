import sys, os, math
import argparse
import general_scripts as gs
import numpy as np
from scipy.optimize import fmin_powell
from scipy.optimize import minimize
from datetime import datetime

def _get_current_function():
    import inspect
    frame=inspect.currentframe()
    return inspect.getframeinfo(frame).function

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

def set_outliers_to_nan(arr, nSig):
    arr[np.fabs(arr-np.nanmean(arr)) > nSig*np.nanstd(arr)]=np.nan
    return arr
        
def convert_log_statistics( base, mean, std ):
    err = 0.5*(np.power(base,mean+std)-np.power(base,mean-std))
    return np.power(base,mean), err

def convert_log_xmgrace( base, mean, std):
    m=np.power(base, mean)
    return m, np.power(base,mean+std)-m, m-np.power(base,mean-std)

# = = = = Debugging
def debug_data_block ( db ):
    print "= = = Debug data_block."
    print "...number of conc series :", len(db)
    for i in range(len(db)):
        print "......series %i has %i values." % (i, len(db[i]))
        print "......first entry:", db[i][0]
    return

def write_xvg_header_simexp(fp, nCurves):
    cmax=31
    s=0
    c=0
#   print >> fp, '@with g0'
    print >> fp, '@ view 0.1, 0.1, 0.7, 0.95'
    print >> fp, '@ legend 0.7, 0.95'
    print >> fp, '@ legend char size 0.75'
    print >> fp, '@ legend length 2'
    for i in range(nCurves):
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

def write_fit_results(fn, nCurves, nTrials, legends, meanLogKd, sigLogKd, meanBaseline, sigBaseline, meanDelta, sigDelta, meanQoF, sigQoF, bTestFail):
    fp=open(fn, 'w')
    print >> fp, "= = = Summary of fitting."
    print >> fp, "= = = Date: %s" % datetime.now().strftime("%Y.%m.%d %H:%M:%S") 
    print >> fp, ""
    print >> fp, "mean QoF over %i trials : %g +- %g" % ( nTrials, meanQoF, sigQoF )
    print >> fp, ''
    for i in range(nCurves):
        print >> fp, "Statistics for curve %i - %s :" % (i+1, legends[i])
        print >> fp, "logKd : %g +- %g" % ( meanLogKd[i], sigLogKd[i] )
        print >> fp, "Kd,delta+,delta- : %g +- %g %g" % convert_log_xmgrace(10, meanLogKd[i], sigLogKd[i])
        if meanBaseline != () :
            print >> fp, "baseline : %g +- %g" % ( meanBaseline[i], sigBaseline[i] )
        else:
            print >> fp, "baseline : %g +- %g" % ( meanBaseline, sigBaseline )
        if meanDelta != () :
            print >> fp, "delta : %g +- %g" % ( meanDelta[i], sigDelta[i] )
        else:
            print >> fp, "delta : %g +- %g" % ( meanDelta, sigDelta )

        print >> fp, "Failed quality check? : %s" % bTestFail[i]
        print >> fp, ''

    fp.close()

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
        print >> fp, "# MEAN +- SIG: %g %g" % ( np.nanmean(vlist), np.nanstd(vlist) )
        print >> fp, "# NaN count: %d" % ( len(vlist[np.isnan(vlist)]) )
        for j in range( shape[-1] ):
            if not np.isnan(vlist[j]):
                print >> fp, "%i %g" % ( j, vlist[j] )
            else:
                print >> fp, "#%i %g" % ( j, vlist[j] )
    elif len(shape)==2:
        for i in range( shape[0] ):
            print >> fp, "@s%i legend \"%s\"" % ( i, legs[i] )
            print >> fp, "# MEAN +- SIG: %g %g" % ( np.nanmean(vlist[i]), np.nanstd(vlist[i]) )
            print >> fp, "# NaN count: %d" % ( len(vlist[i][ np.isnan(vlist[i])]) )
            for j in range( shape[-1] ):
                if not np.isnan(vlist[i,j]):
                    print >> fp, "%i %g" % ( j, vlist[i,j] )
                else:
                    print >> fp, "#%i %g" % ( j, vlist[i,j] )
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
def boundconc_2state(inP, inQ, Kd, bFrac=False):
    s=inP+inQ+np.fabs(Kd)
    PQ=0.5*(s-np.sqrt(s*s-4*inP*inQ))
    if bFrac:
        return PQ/inP
    else:
        return inP-PQ, inQ-PQ, PQ

def model_2state(inP, inQ, Kd, baseline, delta):
    P, Q, PQ = boundconc_2state( inP, inQ, Kd )
    return (P/inP)*baseline + (PQ/inP)*(baseline+delta)

def run_global_estimate( titration ):
    """
    Returns ( baseline, delta, Kd ),
    given a 2D-array titration of dimensions ( nPoints x 3or4 )
    Each titration point contains [P], [Q], Target [, dTarget]
    """
    outKd=None ; outBase=None ; outDelta=None
    Ps = titration[:,0] ; Qs = titration[:,1] ; targets = titration[:,2]
    if len(titration[0])==4:
        error = titration[:,3]
    baseInit = targets[np.argmin(Qs)]
    deltaInit = ( np.max(targets)-np.min(targets) )*np.sign( Qs[np.argmax(targets)]-Qs[np.argmin(targets)] )
    KdInit = estimate_Kd( titration ) 
    chiSqMin = 1e99
    for base in np.linspace(0.5,1.5,9)*baseInit :
        for delta in np.linspace(0.5,1.5,9)*deltaInit :
            for Kd in np.logspace(-2,2,13)*KdInit :
                models = model_2state( Ps, Qs, Kd, base, delta) 
                chiSq  = np.mean( np.power(targets-models,2.0) )
                if chiSq < chiSqMin:
                    outKd=Kd
                    outBase=base
                    outDelta=delta
                    chiSqMin=chiSq
    if outKd is None:
        print "= = WARNING: Global parameter estimate has failed to produce chi values less than 1e99!"
        print "    ...returning basic prediction values base %g ; delta %g ; Kd %g" % (baseInit, deltaInit, KdInit)
        print targets
        sys.exit(1)
        return baseInit, deltaInit, KdInit
    else:
        print "= = Global estimate parameters (%g): base %g ; delta %g ; Kd %g" % (chiSqMin, outBase, outDelta, outKd)
        return outBase, outDelta, outKd

# = = = General fast estimate of KD assuming somewhat well-shaped curve?
def estimate_Kd( dbOne ):
    # Assumes dataBlock as a set of titration-points, with each point as [Rec],[Lig], M, dM
    
    # = = Assuming that this is a consstant-[Rec] titration,
    #     obtain the average of three [Lig-values] closest to the midpoint between min(M) max(M),
    #     as a way of estimating the half-maximum point.
    #     This is ~Kd for theoretical binding where Kd > fixed-[R]
    midM = 0.5*( np.min(dbOne[:,2])+np.max(dbOne[:,2]) )
    listDiff = np.fabs(dbOne[:,2]-midM)
    dbSorted = dbOne[ np.argsort(listDiff) ]
    #print "= = Debug: [Lig,M]:", dbOne[:,1:3]
    #print "= = Debug: [Lig] closest to the midpoint:", dbSorted[0:3,1:3]
    #print "= = Debug: Comparison of new Kd estimate with trivial:", np.mean(dbSorted[0:3,1]), 0.8*np.min(dbOne[:,1])+0.2*np.max(dbOne[:,1])
    return np.mean(dbSorted[0:3,1])
    # = = Trivial estimate. Use the observed range [Lig]
    #return 0.8*np.min(dbOne[:,1])+0.2*np.max(dbOne[:,1])

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
def load_targets_singlefile( fn, concA, bGetSig=False ):
    """
    Load the target metric data and convert to a MxNx3or4 list
    each containing:
    [Receptor] [Ligand] Metric dMetric
    M is the number of titrations curves, N is the number of titration points in each titration.
    """
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

def fit_minimize_API(params, *args):
    """
    Interface function between sharedInformationModel and scipy.optimize.
    Calls the relevant function, assuming that the first argument entry is the model.
    """
    fitModel=args[0]
    return fitModel.fit_minimize_internal(params,args[1:])

class sharedInformationModel:
    """
    Overall class for fitting a sets of curves in which
    -  some parameters are shared between curves, while others are curve specific.
    -  permits multiple replicate fits.
    The collection of parameters as stored as a dictionary of various names.
    curveNames and paramNames is stored as a list so that the order can be conserved, in addition to serving as a key.
    """
    def __init__(self):
        self.paramNames=[]
        self.curveNames=[]
        self.bFailedFits=[]
        self.bShared={}
        self.paramList={}
        self._generate_strShared()
        self._generate_paramList()
        self._update_counts()

    def _generate_strShared(self):
        self.strShared={}
        for key in self.bShared:
            if self.bShared[key]:
                self.strShared[key]="shared "
            else:
                self.strShared[key]="per-curve "

    def _generate_paramList(self):
        self.paramList={}
        for key in self.paramNames:
            if self.bShared[key]:
                self.paramList[key]=None
            else:
                self.paramList[key]=[ None for x in self.curveNames ]

    def _update_counts(self):
        self.nCurves=len(self.curveNames)
        self.nParamTypes=len(self.paramNames)
        self.bFailedFits=[ None for c in self.curveNames ]
        c=0
        for i in range(self.nParamTypes):
            if self.bShared[self.paramNames[i]]:
                c+=1
            else:
                c+=self.nCurves
        self.totParams=c

    def _error_not_implemented( functionName ):
        print >> sys.stderr, "= = ERROR: function %s is not implemented! The subclass has forgotten to implement this." % functionName
        sys.exit(1)


    def resolve_index(self, a, index=None, name=None):
        if not index is None:
            if index >= len(a):
                print >> sys.stderr, "= = ERROR: An illegal index request in resolve_index! %i > size of namelist (%i)" % (index, len(a))
                sys.exit(99)
            name=a[index]
        elif not name is None:
            if name not in a:
                print >> sys.stderr, "= = ERROR: Requsted name not found in resolve_index! %s is not in the list of names!" % (name), a
                sys.exit(99)
            index=np.where(self.paramNames==key)[0]
        else:
            print >> sys.stderr, "= = ERROR: resolve_index encountered None arguments, and has no information to retrieve the matching entry!"
            sys.exit(99)
        return index, name

    def report_parameters(self, fp=sys.stdout):
        print >> fp, "= = Reporting fitting parameters."
        for i in range(self.nParamTypes):
            key=self.paramNames[i]
            print >> fp,  "    ... %s %s:" % (self.strShared[key], key),
            print >> fp, self.paramList[key]

    def print_parameters(self, fp=sys.stdout):
        self.report_parameters( fp )

    def write_param_distribution(self, filePrefix, index=None, name=None, header=''):
        index, name = self.resolve_index(self.paramNames, index, name)
        fileName=filePrefix+"_"+name+".xvg"
        fp=open(fileName,'w')
        print >> fp, header
        print >> fp, '@type xy'
        if self.bShared[key]:
            v=self.sampleParamList[key]
            print >> fp, "@s%i legend \"%s\"" % ( 0, name )
            print >> fp, "# MEAN +- SIG: %g %g" % ( np.nanmean(v), np.nanstd(v) )
            print >> fp, "# NaN count: %d" % ( len(v[np.isnan(v)]) )
            for j in range(len(v)):
                if not np.isnan(v[j]):
                    print >> fp, "%i %g" % ( j, v[j] )
                else:
                    print >> fp, "#%i %g" % ( j, v[j] )
        else:
            for i in range( self.nCurves ):
                v=self.sampleParamList[key][...,i]
                print >> fp, "@s%i legend \"%s\"" % ( i, self.curveNames[i] )
                print >> fp, "# MEAN +- SIG: %g %g" % ( np.nanmean(v), np.nanstd(v) )
                print >> fp, "# NaN count: %d" % ( len(v[np.isnan(v)]) )
                for j in range(len(v)):
                    if not np.isnan(v[j]):
                        print >> fp, "%i %g" % ( j, v[j] )
                    else:
                        print >> fp, "#%i %g" % ( j, v[j] )
                print >> fp, '&'

        fp.close()
        return

    def write_fitted_model(self, fileName, db, bUseMean=False):
        bFailArray=self.bFailedFits
        if bUseMean:
            pList=self.meanParamList
        else:
            pList=self.paramList
        fp = open(fileName, 'w')
        write_xvg_header_simexp(fp, self.nCurves)
        s=0
        for c in range(self.nCurves):
            # = = Print legends first.
            n=self.curveNames[c]
            if not self.bShared['logKd']:
                Kd=np.power(10.0, pList['logKd'][c])
                print >> fp, "@s%i legend \"%s: K\\sD\\N = %5.3g \\xm\\f{}M\"" % (s, n, Kd)
            else:
                Kd=np.power(10.0, np.nanmean(pList['logKd'][c]))
                print >> fp, "@s%i legend \"%s: mean K\\sD\\N = %5.3g \\xm\\f{}M\"" % (s, n, Kd)
            if bFailArray[c]:
                print >> fp, "@s%i line linestyle 2" % s
            s+=2
        #= = Print model values, then target values
        for c in range(self.nCurves):
            x = self.obtain_single_xvalues(db[c])
            model = self.obtain_single_model(db[c], c, bUseMean)
            self.print_model_curve(fp, x, model, syntax='xvg')
            self.print_target_curve(fp, db[c], syntax='xvg')
        fp.close()
        return

    def print_model_curve(self, fp, x, model, syntax='raw'):
        _error_not_implemented( _get_current_function() )

    def print_target_curve(self, fp, dbOne, syntax='raw'):
        _error_not_implemented( _get_current_function() )

    def print_curve(self, fp, x, y, dy=None, syntax='raw', gXvg=None, sXvg=None):
        bError=False
        if not dy is None:
            bError=True
        if syntax=='xvg':
            if (not gXvg is None) and (not sXvg is None):
                print >> fp, "@target g%d.s%d" % (gXvg, sXvg)
            if bError:
                print >> fp, "@type xydy"
            else:
                print >> fp, "@type xy"
        if bError:
            for j in range(len(x)):
                print >> fp, "%g %g %g" % (x[j], y[j], dy[j] )
        else:
            for j in range(len(x)):
                print >> fp, "%g %g" % (x[j], y[j])
        if syntax=='xvg':
            print >> fp, "&"
        return

    def estimate_initial_parameters( self, dataBlock ):
        _error_not_implemented( _get_current_function() )

    def extract_params_direct(self):
        out=[]
        for i in range(self.nParamTypes):
            key=self.paramNames[i]
            out.append(self.paramList[key])
        return [x for x in out]

    def extract_params_as_list(self):
        """
        Extracts the class params for use in an external minimization algorithm.
        Paired with store_params_from_list to complete the API.
        """
        out=[]
        for i in range(self.nParamTypes):
            key=self.paramNames[i]
            if self.bShared[key]:
                out.append(self.paramList[key])
            else:
                for j in range(self.nCurves):
                    out.append(self.paramList[key][j])
        return out
    
    def store_params_from_list(self, params):
        if len(params) != self.totParams:
            print >> sys.stderr, "= = ERROR encountered in store_params_from_list. Number of input parameters is not equal to the number of total degrees of freedom in the fitting model! %i vs %i" % (len(params), self.totParams)
            sys.exit(99)
        c=0
        for i in range(self.nParamTypes):
            key=self.paramNames[i]
            if self.bShared[key]:
                self.paramList[key]=params[c]
                c+=1
            else:
                self.paramList[key]=params[c:c+self.nCurves]
                c+=self.nCurves
    
    def obtain_single_xvalues(self, dbOne):
        """
        Default X-value to plot is the ligand-receptor ratio.
        """
        inP = dbOne[:,0] ; inQ=dbOne[:,1]
        return inQ/inP

    def obtain_single_weight(self, dbOne):
        """
        Default weight is unity.
        """
        return 1.0

    def obtain_single_model(self, dbOne, curveIndex=None, curveName=None, bUseMean=False):
        """
        Default model here is the two-state model, given three named parameters:
        - logKd, apoValue, holoValue
        """
        if bUseMean:
            pList=self.meanParamList
        else:
            pList=self.paramList
        i,n=self.resolve_index(self.curveNames, curveIndex, curveName)
        inP = dbOne[:,0] ; inQ=dbOne[:,1]
        P, Q, PQ = boundconc_2state( inP, inQ, np.power(10.0,pList['logKd']) )
        modelVals = (P/inP)*pList['apoValue'] + (PQ/inP)*pList['holoValue']
        return modelVals
    
    def obtain_single_target(self, dbOne):
        """
        the default model assumes that the single datablock is of entries (N,3).
        - N-titration points, with
        - [P] [L] target
        """
        return dbOne[:,2]

    def chiSq_single_model(self, dbOne, curveIndex=None, curveName=None):
        """
        Default chi-squared function is then to combine the single_model operations to match against the target value.
        Without considering errors.
        """
        i,n=self.resolve_index(self.curveNames, curveIndex, curveName)
        weight=self.obtain_single_weight(dbOne)
        model=self.obtain_single_model(dbOne, curveIndex=i)
        target=self.obtain_single_target(dbOne)
        return np.mean( np.power(target-model, 2.0)/weight)

    def chiSq_all(self, dataBlock):
        chiSq=0
        for i in range(self.nCurves):
            chiSq += self.chiSq_single_model(dataBlock[i],curveIndex=i)
        return chiSq/self.nCurves

    def fit_minimize_internal(self, params, args):
        """
        The minimizer function is the chi-square deivation from the target dataSet.
        Note that the target dataSet is itself stored outside the model
        """
        data=args[0]
        self.store_params_from_list(params)
        if 'logKd' in self.paramNames:
            if np.any( self.paramList['logKd'] > globKDClipMax) or np.any(self.paramList['logKd'] < globKDClipMin ):
                return 1e99
        return self.chiSq_all(data)

    # = = = The following sections pertain to sultiple sampling, useful for error analysis.
    #    ...operations are carried out between the standard paramList and the sampleParamList

    def init_sample_variables(self, nSamples):
        self.sampleParamList={}
        self.meanParamList={}
        self.stdParamList={}
        self.nSamples=nSamples
        for i in range(self.nParamTypes):
            key=self.paramNames[i]
            if self.bShared[key]: 
                self.sampleParamList[key]=np.empty(nSamples)
                self.sampleParamList[key][:]=np.nan
                self.meanParamList[key]=np.nan
                self.stdParamList[key]=np.nan
            else:
                self.sampleParamList[key]=np.empty((nSamples,self.nCurves))
                self.sampleParamList[key][:]=np.nan
                self.meanParamList[key]=np.empty(self.nCurves)
                self.meanParamList[key][:]=np.nan
                self.stdParamList[key]=np.copy( self.meanParamList[key] )

    def store_sample_variables(self, s, l=None):
        """
        Optional l arguement allows direct storage from the input list.
        """
        if not l is None:
            self.store_params_from_list( l )
        for key in self.paramNames:
            self.sampleParamList[key][s]=self.paramList[key]

    def reduce_sample(self):
        for i in range(self.nParamTypes):
            key=self.paramNames[i]
            if self.bShared[key]: 
                self.meanParamList[key]=np.nanmean(self.sampleParamList[key])
                self.stdParamList[key]=np.nanstd(self.sampleParamList[key])
            else:
                self.meanParamList[key]=np.nanmean(self.sampleParamList[key],axis=0)
                self.stdParamList[key]=np.nanstd(self.sampleParamList[key],axis=0)

    def conduct_sample_quality_check(self, dataBlock ):
        """
        Sample quality checks are often model-specific. This is left empty.
        """
        _error_not_implemented( _get_current_function() )


    def report_sample_latest_comparison(self, fp=sys.stdout, bNormalKd=False):
        print >> fp, "= = Comparison of sample statistics with latest fit in paramList:" 
        self.reduce_sample()
        for key in self.paramNames:
            if self.bShared[key]:
                print >> fp, "%s %s statistics:" % (self.strShared[key], key ),
                print >> fp, " single fit: %g" % self.paramList[key] ,
                print >> fp, " mean value from sample: %g +- %g" % \
                        (self.meanParamList[key], self.stdParamList[key])
                if bNormalKd and key == 'logKd':
                    # = = Also report standard Kd
                    print >> fp, "%s %s statistics:" % ("converted", "Kd")
                    print >> fp, "normal Kd single fit: %g" % np.power(10,self.paramList[key])
                    meanKd, sigKd = convert_log_statistics(10.0, \
                            self.meanParamList[key], self.stdParamList[key])
                    print >> fp, "normal Kd mean value %g +- %g" % ( meanKd, sigKd )
                    meanKd, plus, minus = convert_log_xmgrace(10.0, \
                            self.meanParamList[key], self.stdParamList[key])
                    print >> fp, "normal Kd mean xmgrace : %g %g %g" % (meanKd, plus, minus)

        for c in range(self.nCurves):
            print >> fp, "= = titration %s (#%d):" % (self.curveNames[c], c+1)
            if not self.bFailedFits[c] is None:
                print >> fp, "= = titration %s Failed quality check? %s" % (self.curveNames[c], self.bFailedFits[c])
            for key in self.paramNames:
                if not self.bShared[key]:
                    print >> fp, "%s %s statistics:" % (self.strShared[key], key ),
                    print >> fp, " single fit: %g" % self.paramList[key][c],
                    print >> fp, " mean value from sample: %g +- %g" % \
                            (self.meanParamList[key][c], self.stdParamList[key][c])
                    if bNormalKd and key == 'logKd':
                        # = = Also report standard Kd
                        print >> fp, "%s %s statistics:" % ("converted", "Kd")
                        print >> fp, "normal Kd single fit: %g" % np.power(10,self.paramList[key][c])
                        meanKd, sigKd = convert_log_statistics(10.0, \
                                self.meanParamList[key][c], self.stdParamList[key][c])
                        print >> fp, "normal Kd mean value %g +- %g" % ( meanKd, sigKd )
                        meanKd, plus, minus = convert_log_xmgrace(10.0, \
                                self.meanParamList[key][c], self.stdParamList[key][c])
                        print >> fp, "normal Kd mean xmgrace : %g %g %g" % (meanKd, plus, minus)


    def report_sample_stats(self, fp=sys.stdout, bNormalKd=False):
        self.reduce_sample()
        for key in self.paramNames:
            print >> fp, "    %s %s:" % (self.strShared[key], key )
            if self.bShared[key]: 
                print >> fp, "    -- mean value from sample: %g +- %g" % \
                        (self.meanParamList[key], self.stdParamList[key])
            else:
                for c in range(self.nCurves):
                    print >> fp, "    -- mean value from sample: %g +- %g" % \
                        (self.meanParamList[key][c], self.stdParamList[key][c])
            if key == 'logKd':
                # = = Also report standard Kd:
                if self.bShared[key]: 
                    meanKd, sigKd = convert_log_statistics(10.0, self.meanParamList[key], self.stdParamList[key])
                    print >> fp, "= = Overall titration raw Kd estimate: %g +- %g" % ( meanKd, sigKd )
                    meanKd, plus, minus = convert_log_xmgrace(10.0, self.meanParamList[key], self.stdParamList[key])
                    print >> fp, "    ... for xmgrace Kd : %g %g %g" % (meanKd, plus, minus)
                else:
                    for c in range(self.nCurves):
                        meanKd, sigKd = convert_log_statistics(10.0, \
                                self.meanParamList[key][c], self.stdParamList[key][c])
                        print >> fp, "= = Titration %s (#%i) raw Kd estimate: %g +- %g" % \
                            ( self.curveNames[c], c+1, meanKd, sigKd )
                        meanKd, plus, minus = convert_log_xmgrace(10.0, \
                                self.meanParamList[key][c], self.stdParamList[key][c])
                        print >> fp, "    ... for xmgrace Kd : %g %g %g" % (meanKd, plus, minus)

    def write_sample_fit_results(self, fileName):
        fp=open(fileName, 'w')
        print >> fp, "= = = Summary of fitting."
        print >> fp, "= = = Date: %s" % datetime.now().strftime("%Y.%m.%d %H:%M:%S") 
        print >> fp, "= = = Using "
        print >> fp, "mean QoF over %i trials : %g +- %g" % ( nTrials, meanQoF, sigQoF )
        print >> fp, ''
        self.write_sample_latest_comparison( fp, bNormalKd=True )
        fp.close() 

def report_parameters(deltas, baselines, Kds, adjective=''):
    print "= = Reporting %s Deltas:" % adjective, deltas
    print "= = Reporting %s baselines:" % adjective, baselines
    print "= = Reporting %s Kds:" % adjective, Kds

# = = = Model C = = =
# Single baseline but multiple deltas to capture different binding modes.
class model_sharedApo(sharedInformationModel):
    """
    Initialisation requires the names of the ligands to be titrated. Can be a dummy list.
    The internal parameters for this model are, in order:
    - apoValue, shared
    - holoValue, unique
    - logKd, unique
    The desired datablock is a list of M ligand titrations:
    - each a 2D-array with any number of titration points.
    - [ [concP, concL, value], [concP, concL, value], ... ]
    """
    def __init__(self, curveNames ):
        self.paramNames=['apoValue','holoValue','logKd']
#        self.paramNames=['holoValue','apoValue','logKd']
        self.bShared={'apoValue': True, 'holoValue': False, 'logKd': False }
        self.curveNames=list(curveNames)
        self._generate_strShared()
        self._generate_paramList()
        self._update_counts()

    def print_model_curve(self, fp, x, model, syntax='raw'):
        self.print_curve(fp, x, model, syntax='xvg')

    def print_target_curve(self, fp, dbOne, syntax='raw'):
        x = self.obtain_single_xvalues(dbOne)
        if dbOne.shape[-1]==4:
            self.print_curve(fp, x, dbOne[:,2], dy=dbOne[:,3], syntax='xvg')
        else:
            self.print_curve(fp, x, dbOne[:,2], syntax='xvg')

    def obtain_single_weight(self, dbOne):
        """                                   
        Divide by the uncertainties, if the entry has them. Else return 1.0.
        """
        if dbOne.shape[-1]==4:
            return dbOne[:,3]
        else:
            return 1.0

    def obtain_single_model(self, dbOne, curveIndex=None, curveName=None, bUseMean=False):
        """
        Default model here is the two-state model, given three named parameters:
        - logKd, apoValue, holoValue
        """
        if bUseMean:
            pList=self.meanParamList
        else:
            pList=self.paramList
        i,n=self.resolve_index(self.curveNames, curveIndex, curveName)
        inP = dbOne[:,0] ; inQ=dbOne[:,1]
        P, Q, PQ = boundconc_2state( inP, inQ, np.power(10.0,pList['logKd'][i]) )
        modelValues = (P/inP)*pList['apoValue'] + (PQ/inP)*pList['holoValue'][i]
        return modelValues

    def estimate_initial_parameters(self, dataBlock, bMinimizer=True):
        apoList=[] 
        for i in range(len(dataBlock)):
            a, h, logKd = self._global_search_one_curve( dataBlock[i], curveIndex=i )
            apoList.append(a)
            self.paramList['holoValue'][i]=h
            self.paramList['logKd'][i]=np.clip(logKd, globKDClipMin, globKDClipMax)
        self.paramList['apoValue']=np.mean(apoList)

        if bMinimizer:
            return self.extract_params_as_list()

    def _estimate_logKd( self, dbOne):
        # = = Assuming that this is a consstant-[Rec] titration,
        #     obtain the average of three [Lig-values] closest to the midpoint between min(M) max(M),
        #     as a way of estimating the half-maximum point.
        #     This is ~Kd for theoretical binding where Kd > fixed-[R]
        midM = 0.5*( np.min(dbOne[:,2])+np.max(dbOne[:,2]) )
        listDiff = np.fabs(dbOne[:,2]-midM)
        dbSorted = dbOne[ np.argsort(listDiff) ]
        return np.log10( np.mean(dbSorted[0:3,1]) )

    def _global_search_one_curve( self, dbOne, curveIndex ):
        logKdOut=None ; apoOut=None ; holoOut=None
        Ps = dbOne[:,0] ; Qs = dbOne[:,1] ; targets = dbOne[:,2]
        if len(dbOne[0])==4:
            error = dbOne[:,3]
        apoInit  = targets[np.argmin(Qs)]
        holoInit = targets[np.argmax(Qs)]
        logKdInit = self._estimate_logKd( dbOne ) 
        chiSqMin = 1e99
        for apo in np.linspace(0.5,1.5,9)*apoInit :
            self.paramList['apoValue']=apo
            for holo in np.linspace(0.5,1.5,9)*holoInit :
                self.paramList['holoValue'][curveIndex]=holo
                for logKd in np.linspace(-2,2,9)*logKdInit :
                    self.paramList['logKd'][curveIndex]=logKd
                    chiSq = self.chiSq_single_model(dbOne, curveIndex)
                    if chiSq < chiSqMin:
                        logKdOut=logKd ; apoOut=apo ; holoOut=holo
                        chiSqMin=chiSq
        if logKdOut is None:
            print "= = WARNING: Global parameter estimate has failed to produce chi values less than 1e99!"
            print "    ...returning basic prediction values apo %g ; holo %g ; Kd %g" % \
                    (apoInit, holoInit, logKdInit)
            return apoInit, holoInit, logKdInit
        else:
            print "= = Global estimate parameters (%g): apo %g ; holo %g ; Kd %g" % \
                    (chiSqMin,  apoOut, holoOut, logKdOut)
            return apoOut, holoOut, logKdOut
    
    def obtain_datablock_apoVariation(self, dataBlock ):
        val=[] ; err=[]
        for titr in dataBlock:
            for pt in titr:
                if pt[1]==0.0:
                    val.append( pt[2] )
                    if len( pt )==4:
                        err.append( pt[3] )
        if len( val )>0:
            ret=np.std( val )
            if len( err )>0:
                ret=np.sqrt( ret**2.0+np.mean(err)**2.0 )
            return ret if ret>0.0 else None
        else:
            return None

    def obtain_datablock_observedRange(self, dataBlock ):
        max=[ np.max( dbOne[...,2] ) for dbOne in dataBlock ]
        min=[ np.min( dbOne[...,2] ) for dbOne in dataBlock ]
        return np.max( max ) - np.min( min )

    def conduct_sample_quality_check(self, dataBlock ):
        """
        In this type of model, we can check for the following:
        - the total perturbation compared to the variation in apo-measurements.
        - the fitted perturbation as a function of total measured difference.
        """
        apoVariation  = self.obtain_datablock_apoVariation( dataBlock )
        observedRange = self.obtain_datablock_observedRange( dataBlock )
        if apoVariation is None:
            print "= = WARNING: cannot conduct apo-variations test as there appears to be only zero/one apo-measurements!"
            apoVariation=0.0

        # = = = Checks for the sample fit. Eliminate outliers for each
        self.reduce_sample()
        # = = = sample shape: nTrials x nCurves ; meanShape: nCurves
        sampleDeltas=self.sampleParamList['holoValue']-self.sampleParamList['apoValue'][:,np.newaxis]
        meanDeltas=self.meanParamList['holoValue']-self.meanParamList['apoValue']
        sigDeltas= np.sqrt( np.square(self.stdParamList['holoValue']) + np.square(self.stdParamList['apoValue']) )
        nanIndices = np.greater( np.fabs(sampleDeltas-meanDeltas), globRejectThreshold*sigDeltas )
        print "= = NOTE: Using criteria %g sigma,  %d outliers will be rejected from %i by %i trials" % \
                ( globRejectThreshold, len(nanIndices[nanIndices]), self.nSamples, self.nCurves )
        self.sampleParamList['holoValue'][nanIndices]=np.nan
        self.sampleParamList['logKd'][nanIndices]=np.nan
        # = = = Update after np.nans have been assigned.
        self.reduce_sample()

        # = = = Checks for the single fit.
        deltas=self.paramList['holoValue']-self.paramList['apoValue']
        for c in range(self.nCurves):
            self.bFailedFits[c]=False 
            delta=deltas[c]
            if np.fabs(delta) < apoVariation :
                print "= = = WARNING: Curve %d does not show significant fitted delta-deviations from variations in apo-measurements!" % (c+1)
                print "    ... | %g | < %g" % ( delta, apoVariation )
                #self.paramList['holoValue'][c] = self.paramList['apoValue']
                #self.paramList['logKd'][c] = globKDClipMax
                self.bFailedFits[c]=True 
            if np.fabs(delta) > globObservedRangeLimitMult*observedRange:
                print "= = = WARNING: Curve %d shows unbelievably large delta compared to the limits of target values! This indicates a bad fit." % (c+1)
                print "    ... | %g | >> %g" % ( delta, observedRange )
                #outDelta[i]   = 0.0
                #meanKd[i] = np.power(10,KDClipMax)
                self.bFailedFits[c] = True

        return self.bFailedFits

# The order of parameters matter for fmin_powell as it searches sequentially along each dimension
def extract_params_modelC(params, nCurves):
    #baseline=params[0]
    #deltas=params[1:nCurves+1]
    deltas=params[0:nCurves]
    baseline=params[nCurves]
    Kds=params[1+nCurves:1+2*nCurves]
    return deltas, baseline, Kds

def concat_params_modelC(deltas, baseline, Kds):
#    return [baseline] + list(deltas) + list(Kds)
    return list(deltas) + [baseline] + list(Kds)

# Specific for the format of the datablock above.
def estimate_initial_parameters_modelC( datablock ):
    baseline = []
    deltas = []
    Kds = [] # The other params here are Kds of individial series.
    for i in range(len(datablock)):
        # Get locations where [Lig] is the lowest or highest
        #oneTitration=datablock[i]
        #Ps = titration[:,0] ; Qs = titration[:,1] ; targets = titration[:,2]
        #baseInit = targets[np.argmin(Qs)]
        #deltaInit = ( np.max(targets)-np.min(targets) )*np.sign( Qs[np.argmax(targets)]-Qs[np.argmin(targets)] )
        #KdInit = estimate_Kd( titration ) 
        base, delta, Kd = run_global_estimate( datablock[i] )
        baseline.append( base ) 
        deltas.append( delta )
        Kds.append( Kd )
    params = concat_params_modelC( deltas, np.mean(baseline), Kds )
    return params

def write_fitted_modelC(outfn, legends, params, db, bFailArray=None):
    nCurves=len(db)
    deltas, baseline, Kds = extract_params_modelC(params, nCurves)
    Kds = np.fabs( Kds )
    report_parameters( deltas, baseline, Kds, 'averaged model')
    if bFailArray is None:
        bFailArray=[ False for x in Kds ]

    fp = open(outfn, 'w')
    write_xvg_header_simexp(fp, nCurves)
    if len(db[0][0])==3:
        bErr=False
    elif len(db[0][0])==4:
        bErr=True
    #write_xvg_legends(fp, "K\\sd\\N =", "\\xm\\f{}M", Kds, 2)
    s=0
    for i in range(len(Kds)):
#        print >> fp, "@s%i legend \"%s: S\\s\\xD\\f{}\\N = %5.3g \\xm\\f{}M\"" % (s, legends[i], Kds[i])
        print >> fp, "@s%i legend \"%s: K\\sD\\N\\S\\eff.\\N = %5.3g \\xm\\f{}M\"" % (s, legends[i], Kds[i])
        if bFailArray[i]:
            print >> fp, "@s%i line linestyle 2" % s
        s+=2
    for i in range(len(db)):
        print >> fp, "@type xy"
        for j in range(len(db[i])):
            inP = db[i][j][0] ; inQ = db[i][j][1] ; targ=db[i][j][2]
            P, Q, PQ = boundconc_2state( inP, inQ, Kds[i] )
            model_val = (P/inP)*baseline + (PQ/inP)*(baseline+deltas[i])
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

# Fit function for modelC
def fitfunc_modelC(params, *args):
    nCurves=args[0]
    data=args[1]
    deltas, baseline, Kds = extract_params_modelC(params, nCurves)
    chi2=0.0
    count=0
    for i in range(len(data)):
        for j in range(len(data[i])):
            inP, inQ, target = data[i][j]
            #print "...minimising...", inP, inQ, target
            P, Q, PQ = boundconc_2state( inP, inQ, Kds[i] )
            model_val = (P/inP)*baseline + (PQ/inP)*(baseline+deltas[i])
            chi2  += ( target - model_val )**2.0
            count += 1
    return 1.0*chi2/count

# = = = = = Model D: as with model B, but with an additional contribution from unbound ligand
def extract_params_modelD(params, nCurves):
    delta_cplx=params[0]
    #delta_lig=params[1]
    bases=params[1:2+nCurves]
    Kds=params[1+nCurves:1+2*nCurves]
    return delta_cplx, bases, Kds
    #delta_lig=params[1:1+nCurves]
    #bases=params[1+nCurves:1+2*nCurves]
    #Kds=params[1+2*nCurves:1+3*nCurves]
    #return delta_cplx, delta_lig, bases, Kds

# Specific for the format of the datablock above.
def estimate_initial_parameters_modelD( datablock, apo_conc ):
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
        Kds.append( estimate_Kd( db ) )
        #deltas_lig.append( -0.1*deltas[-1] )

    #params = [ np.mean(deltas) ] + [ np.mean(deltas)*-0.1 ] + mins + Kds
    params = [ np.mean(deltas) ] + mins + Kds
    return params

def write_fitted_modelD(outfn, legends, params, db):
    nCurves=len(db)
    delta_cplx, bases, Kds = extract_params_modelD(params, nCurves)
    Kds = np.fabs( Kds )
    fp = open(outfn, 'w')
    write_xvg_header_simexp(fp, nCurves)
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

# Fit function for modelD
def fitfunc_modelD(params, *args):
    nCurves=args[0]
    data=args[1]
    delta_cplx, bases, Kds = extract_params_modelD(params, nCurves)
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

def extract_params_modelD(params, nCurves):
    deltas=params[0:nCurves]
    Kd=params[nCurves]
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
        deltas.append( (setmax-setmin)*np.sign(db[loc2,1]-db[loc1,1]) )
        Kds.append( estimate_Kd( db ) )
    params = list(deltas) + [ np.mean(Kds) ]
    return params

def write_fitted_modelD(outfn, legends, params, db, conc_rec):
    nCurves=len(db)
    deltas, Kd = extract_params_modelD(params, nCurves)
    Kd = np.fabs( Kd )
    fp = open(outfn, 'w')
    write_xvg_header_simexp(fp, nCurves)
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
            model_val = boundconc_2state( inP, inQ, Kd, bFrac=True ) * deltas[i]
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
    nCurves=args[0]
    data=args[1]
    deltas, Kd = extract_params_modelD(params, nCurves)
    chi2=0.0
    count=0
    for i in range(len(data)):
        for j in range(len(data[i])):
            inP, inQ, target = data[i][j]
            #print "...minimising...", inP, inQ, target
            model_val = boundconc_2state( inP, inQ, Kd, bFrac=True ) * deltas[i]
            chi2  += ( target - model_val )**2.0
            count += 1
    return 1.0*chi2/count

# = = NMR chemical-shift perturbation based titration.
#     This will utilise a slightly more advanced CSP model than assuming everything is in fast exchange.

#def estimate_initial_parameters_NMR_CSP( dataBlock ):
#def fitfunc_modelCSP( params, args=(nCurves, data_block):
#write_fitted_modelCSP( fileModel, legends, xopt, data_block, concRec )

def extract_params_modelCSP(params, nCurves):
    deltas=params[0:nCurves]
    taus =params[nCurves:2*nCurves]
    Kd=params[2*nCurves]
    return deltas, taus, Kd

def concat_params_modelCSP(deltas, taus, Kd):
    return list( deltas ) + list(taus) + [ Kd ]

def fitfunc_modelCSP(params, *args):
    nCurves=args[0]
    data=args[1]
    deltas, taus, Kd = extract_params_modelCSP(params, nCurves)
    chiSq=0.0 ; count=0
    for i in range(len(data)):
        for j in range(len(data[i])):
            inP, inQ, target = data[i][j]
            boundFrac = boundconc_2state( inP, inQ, Kd, bFrac=True ) 
            modelValue = mix_CSP( 0.0, deltas[i], boundFrac, taus[i] )
            # = = = Will return three values in medium/slow-exchange. Pick the closest one of three.
            chiSq  += np.min( np.power( target - modelValue, 2.0) )
            count += 1
    return 1.0*chiSq/count

def estimate_initial_parameters_modelCSP( dataBlock ):
    deltas = [] ; taus = []
    Kds = [] # The other params here are Kds of individial series.
    for i in range(len(dataBlock)):
        db=dataBlock[i]
        loc1=np.argmin( db[:,1] )
        setmin = db[loc1,2]
        loc2=np.argmax( db[:,1] )
        setmax = db[loc2,2]
        deltas.append( (setmax-setmin)*np.sign(db[loc2,1]-db[loc1,1]) )
        # = = = Assume 0.5 of tau_critical, although I think the units are wrong?
        taus.append( 0.5*np.sqrt(2)/np.pi/deltas[-1] )
        Kds.append( estimate_Kd( db ) )
    params = list(deltas) + list(taus) + [ np.mean(Kds) ]
    return params

def write_fitted_modelCSP( outFile, legends, params, dataBlock, concRec ):
    nCurves=len(dataBlock)
    deltas, taus, Kd = extract_params_modelCSP(params, nCurves)
    Kd = np.fabs( Kd )
    fp = open(outFile, 'w')
    write_xvg_header_simexp(fp, nCurves)
    if len(dataBlock[0][0])==3:
        bErr=False
    elif len(dataBlock[0][0])==4:
        bErr=True
    print >> fp, "@subtitle \"K\\sD\\N: %5.3g \\xm\\f{}M ; [rec]: %g \\xm\\f{}M \"" % ( Kd, concRec )
    s=0
    for i in range(len(deltas)):
        print >> fp, "@s%i legend \"%s\"" % (s, legends[i])
        s+=2
    for i in range(len(dataBlock)):
        # = = Model curve = =
        print >> fp, "@type xy"
        for j in range(len(dataBlock[i])):
            inP = dataBlock[i][j][0] ; inQ = dataBlock[i][j][1] ; targ=dataBlock[i][j][2]
            boundFrac = boundconc_2state( inP, inQ, Kd, bFrac=True )
            modelValue = mix_CSP( 0.0, deltas[i], boundFrac, taus[i] )
            if len(modelValue)>1:
                loc = np.argmin( np.power(modelValue - dataBlock[i][j][2],2.0) )
            else:
                loc = 0
            # = = = Will return three values in medium/slow-exchange. Pick the closest one of three.
            print >> fp, "%g %g" % (inQ/inP, modelValue[loc])
        print >> fp, "&"
        # = = Target Curve = =
        if bErr:
            print >> fp, "@type xydy"
        else:
            print >> fp, "@type xy"
        for j in range(len(dataBlock[i])):
            inP = dataBlock[i][j][0] ; inQ = dataBlock[i][j][1] ; targ=dataBlock[i][j][2]
            if bErr:
                print >> fp, "%g %g %g" % (inQ/inP, dataBlock[i][j][2], dataBlock[i][j][3])
            else:
                print >> fp, "%g %g" % (inQ/inP, dataBlock[i][j][2])
        print >> fp, "&"
    fp.close()

def mix_CSP(wFree, wBound, pBound, tauBound):
    """
    Computes the chemical shift of the peak in a mixture of
    free and bound populations in the following 2-state reaction:
    [ FreeProtein ] + ligand <-> [ BoundProtein ]
    Based on London, J. Magn. Reson., 1993 ;
    who further cites Lyondenn-Bell, Prog. NMR. Spectrosc., 1967.

    Arguments:
    - wFree and wBound are the chemical shifts of the free and bound states
    - pBound is the fraction of bound state population,
      where pFree + pBound = 1
    - tauBound is the lifetime of the bound population,
      is therefore the inverse of k^-1

    Returns the rotts of the cubic solution, which will either be
    one or three chemical shift values depending on the input.

    Derivation:
    According to Eq. 1 (London, 1993), the spectral intensity of the mix is:
            p_f . p_b . (w_f - w_b)^2
    I(w) = --------------------------                                 ,
            [ tau(w_f - w) (w_b - w) ]^2 + [p_f w_f + p_b w_b - w ]^2 
    where p_f + p_b = 1,
    p_f / p_b = tau_f / tau_b,
    thus: tau = tau_b . tau_f / (tau_b + tau_f) = (1 - p_b).tau_b
    Note transition pointsnear tau_c = sqrt(2) / ( pi dv)

    The peak positions are obtained by setting d I(w) / dw = 0

    This implies solving:
      [2 tau^2] w^3 - [3 tau^2(w_f+w_b)] w^2
    + [tau^2( (w_f + w_b)^2 + 2 w_f w_b ) + 1] w
    - [w_f w_b (w_f + w_b) + p_fw_f + p_b w_b]
    = 0 
    We will use np.roots to solbe the degree-3 polynomial
    """
    tauSq =  np.power( (1.0-pBound)*tauBound, 2.0)
    fbSum =  wFree+wBound
    A     =  2.0*tauSq
    B     = -3.0*tauSq*fbSum
    C     =  tauSq*(fbSum*fbSum+2.0*wFree*wBound)+1.0
    D     = -1.0*(fbSum*wFree*wBound+(1.0-pBound)*wFree+pBound*wBound)

    wOut = np.roots([A,B,C,D])
    # = = Return only the real solutions.
    return [ w.real for w in wOut if w.imag == 0 ]

# = = = = = End MNR CSP section.

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
                                'In mode two, the measure for each rec:lig is a single value, agains with the assumption that the complexation'
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
parser.add_argument('--conc_rec', type=float, dest='concRec', default=None,
                    help='Mode 2: The concentration of the receptor required to interpret the X-column in --targ.')
parser.add_argument('--model', type=str, default='A',
                    help='Model selection for shared information between titrations.'
                    'A: Fit single baseline and single delta, representing a single binding mode to the same receptor.'
                    'B: Fit multiple baseline but single delta, representing possibility of some contamination or baseline drift.'
                    'C: Fit single baseline but multiple deltas, representing the possibility of multiple binding modes to the same receptor.')
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
                    '1-point : to remove data at one ratio at a time and analyse the difference. This method does not require uncertainty estimates at each point.'
                    'noise1sig / noise2sig : add random noise to the mean value scaled to the uncertainty,'
                    'where 1sig and 2sig represents taking either the 64%% or 95%% confidence interval for a guassian noise model.')
parser.add_argument('--nTrials', type=int, dest='ntrials', default=50,
                    help='In Error analyses that require trials, run this many trials. Applicable to, e.g., noise1sig.')
parser.add_argument('--rejectThreshold', type=float, dest='outlierThreshold', default=2.0,
                    help='In Error analyses that require trials, reject trials where the results are more than N standard-deviations away from the set mean. Applicable to, e.g., noise1sig.')
parser.add_argument('--strongKdLimit', type=float, default=-2.0, 
                    help='Temporary sensitivity limits to reject fits that don\'t make sense. Overwritten when a constant receptor concentration is given.')
parser.add_argument('--weakKdLimit', type=float, default=5.0,
                    help='Temporary sensitivity limits to reject fits that don\'t make sense. Overwritten when a constant receptor concentration is given.')

args = parser.parse_args()

opref=args.opref
fileScore1=args.opref+'_score1.xvg'
fileScore2=args.opref+'_score2.xvg'
fileModel=args.opref+'_model.xvg'
fileHist=args.opref+'_hist.xvg'

bLog=args.bLog
bWeight=args.bSigma
bLig=False

errMode=args.err_type
if errMode=='noise1sig' or errMode=='noise2sig':
    bReadSig=True
else:
    bReadSig=False

globRejectThreshold=args.outlierThreshold
if not args.concRec is None:
    logConc = np.log10( args.concRec )
    globKDClipMin = logConc-4.0
    globKDClipMax = logConc+4.0
else:
    globKDClipMin = args.strongKdLimit
    globKDClipMax = args.weakKdLimit

globObservedRangeLimitMult = 10.0

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
    # This assumes that the receptor concentration is held constant throughout the titration series.
    # = = =
    # Pose data array as 3D entries [ titratio-series, titration-points, 2or3]
    # X-values in our target are ligand-receptor ratios. Y-values and dY values  M +- dM
    # This is converted to [Rec] [Lig] M dM
    data_block, legends = load_targets_singlefile( efiles, concRec, bReadSig )
    #data_block, err_block = split_datablock( data_block )
    debug_data_block( data_block )

    # = = Begin setup for quality checks.
    #   - First gather the scatter of apo measurements to test if model deviation is too small.
    #   - Then gather the maximum observed range in titration to test model deviation is too large.
#    apoMeasures=[]
#    valueRange=[]
#    for titr in data_block:
#        for pt in titr:
#            if pt[1]==0.0:
#                apoMeasures.append( pt[2] )
#            valueRange.append( pt[2] )
#    if len(apoMeasures) > 1:
#        bTestDeltaSignificance = True
#        apo2sigma = 2.0*np.std( apoMeasures )
#        print "= = Quality check section: from %d repeats, the 2-sigma of apo measurement deviations is %g" % ( len(apoMeasures), apo2sigma )
#        if apo2sigma > 0:
#            print "= = Turning on Delta-Significance testing."
#            bTestDeltaSignificance = True
#    else:
#        bTestDeltaSignificance = False
#        print "= = NOTE: no or only one receptor-apo measurements have been found. Cannot do significance checks."
#    valueBounds = np.max( valueRange ) - np.min( valueRange )
#    del apoMeasures, valueRange
    # = = = End set up for quality checks

    # = = = Big IF block to process different kinds of applications.
    bTestGoodnessOfFit = False
    if args.b3Body:
        # Adopt a model that sets a different baseline for each titration series,
        # but the change in signal upon binding is conserved, as well as the unbound ligand signal.
        # Reminder: Parameters are, 2 deltas, N minimums, and N K_D for each concentration series in datablock.
        # i.e.: 2N+2 parameters to be fillted.
        print >> sys.stderr, "= = = ERROR: 3-body fitting is currently broken."
        sys.exit(1)
        params = estimate_initial_parameters_modelD( data_block, concRec )
        nCurves = len(data_block)
        print "= = = Initial parameter estimates:"
        print params
        dummy=0
        fminOut = fmin_powell(fitfunc_modelD, params, args=(nCurves, data_block),
                              full_output=True)
        xopt    = fminOut[0]
        funcopt = fminOut[1]
        print "= = Optimised parameters: "
        print xopt
        write_fitted_modelD( fileModel, legends, xopt, data_block )
        print "= = = Written model and targets to %s ." % fileModel
    elif args.bNMRTitration:
    	# Adopt NMR titration mode using chemical shift perturbation.
        # This means computing both an on/off rate as well as the overall affinity Kd
        # the range of parameters are: [N*delta_max, N*tauBound, Kd ]
        params = estimate_initial_parameters_modelCSP( data_block )
        nCurves = len(data_block)
        print "= = = Initial parameter estimates:"
        print params
        fminOut = fmin_powell(fitfunc_modelCSP, params, args=(nCurves, data_block),
                              full_output=True)
        xopt    = fminOut[0]
        funcopt = fminOut[1]
        print "= = Optimised parameters: "
        print xopt
        write_fitted_modelCSP( fileModel, legends, xopt, data_block, concRec )
        print "= = = Written model and targets to %s ." % fileModel

    elif True:
        # = = = Define fitting model here!
        fitModel=model_sharedApo(legends)

        if errMode=='none':
            params = fitModel.estimate_initial_parameters( data_block )
            print "= = = Initial parameter estimates:"
            fitModel.report_parameters()
            fminOut = fmin_powell(fit_minimize_API, params, args=(fitModel, data_block), full_output=True)
            xopt    = fminOut[0] ; funcopt = fminOut[1]
            print "= = Optimised parameters:"
            print xopt
            fitModel.store_params_from_list( xopt )
            fitModel.write_fitted_model( fileModel, data_block )
            print "= = = Written model and targets to %s ." % fileModel
            sys.exit()
        elif errMode=='1-point':
            bTestGoodnessOfFit = True
            print "= = Error mode requested: single-point removal"
            params = fitModel.estimate_initial_parameters( data_block )
            print "= = = Initial parameter estimates:"
            fitModel.report_parameters()
            nCurves = fitModel.nCurves
            conc_list=obtain_PQ_conc_from_datablock( data_block )
            nConcs = len(conc_list)
            print "= = List of concentrations detected:"
            print conc_list
            fitModel.init_sample_variables(nConcs)
            fitQuality=np.zeros( nConcs )
            # = = = Estimate uncertainty here by removing a single titration point.
            for i, c in enumerate(conc_list):
                dbLoc=remove_conc_from_datablock(data_block, c)
                print "= = Round %d : removing concentration point %s .." % (i+1, str(c))
                print "= = ... current number of points in each titration set:", [ len(x) for x in dbLoc ]
                fminOut = fmin_powell(fit_minimize_API, params, args=(fitModel, dbLoc), full_output=True)
                fitQuality[i] = fminOut[1]
                #fitModel.store_params_from_list( xopt )
                fitModel.store_sample_variables( i, l = fminOut[0] )
                fitModel.report_parameters()
                #logKd[i] = np.clip(np.log10(np.fabs(Kds)),KDClipMin,KDClipMax) ; print "= = log10-Kd:", logKd[i]
                #print ""

            # = = = Compute the value from the full curve.
            fminOut = fmin_powell(fit_minimize_API, params, args=(fitModel, data_block), full_output=True)
            fitModel.store_params_from_list( fminOut[0] )
            fitModel.report_sample_latest_comparison()
            fitModel.conduct_sample_quality_check( data_block )
            fitModel.write_fitted_model( fileModel, data_block, bUseMean=False )
            # = = = Write Summary data.
            fp = open( opref+'_results.txt', 'w')
            fitModel.report_sample_latest_comparison( fp, bNormalKd=True )
            fp.close()
            for key in fitModel.paramNames:
                fitModel.write_param_distribution( opref, name=key)
            print "= = = Complete!"
            sys.exit()
        elif errMode=='noise2sig' or errMode=='noise1sig':
            bTestGoodnessOfFit = True
            print "= = Error mode requested: noise-addition "
            if errMode=='noise2sig':
                efact=0.5
            else:
                efact=1.0
            params = fitModel.estimate_initial_parameters( data_block )
            print "= = = Initial parameter estimates:"
            fitModel.report_parameters()
            nCurves=fitModel.nCurves
            nTrials=args.ntrials
            fitQuality=np.zeros( nTrials )
            fitModel.init_sample_variables( nTrials )
            bFirst=True
            for i in range(nTrials):
                print "= = = Conducting trial %i of %i ..." % (i+1, nTrials)
                dbLoc=gaussian_noise_datablock(data_block, scale=efact)
                fminOut = fmin_powell(fit_minimize_API, params, args=(fitModel, dbLoc), full_output=True)
                fitQuality[i] = fminOut[1]
                #fitModel.store_params_from_list( xopt )
                fitModel.store_sample_variables( i, l = fminOut[0] )
                fitModel.report_parameters()
                print "    ...trial %i goodness-of-fit: %g" % (i+1, fitQuality[i])

            # = = = Compute the value from the full curve.
            fminOut = fmin_powell(fit_minimize_API, params, args=(fitModel, data_block), full_output=True)
#                    direc=np.fabs(params)*0.1 )
            fitModel.store_params_from_list( fminOut[0] )
            fitModel.report_sample_latest_comparison()
            fitModel.conduct_sample_quality_check( data_block )
            fitModel.write_fitted_model( fileModel, data_block, bUseMean=True )
            # = = = Write Summary data.
            fp = open( opref+'_results.txt', 'w')
            fitModel.report_sample_latest_comparison( fp, bNormalKd=True )
            fp.close()
            for key in fitModel.paramNames:
                fitModel.write_param_distribution( opref, name=key)
            print "= = = Complete!"
            sys.exit()

            # = = TODO OLD. Convert to mean statistics with outlier discards.
            print "= = Conducting mean/standard analysis to discard statistical outliers."
            print "    ...rejecting any trials in which any of deltas, baseline, or, logKd are outside %i standard-deviations" % outlierThreshold

            badDeltaIndices = np.greater(np.fabs(deltas), 10.0*valueBounds)
            print "= = NOTE: %d bad deltas will be rejected across %i by %i trial" % ( len(badDeltaIndices[badDeltaIndices]), nCurves, nTrials )
            print "    ...as they are at last an order of magnitude greater than the observed spread of all measurements."
            deltas[badDeltaIndices]=np.nan ; logKd[badDeltaIndices]=np.nan

            meanDelta = np.nanmean(deltas,axis=0) ; sigDelta = np.nanstd(deltas, axis=0)
            #meanBaseline = np.nanmean(baselines,axis=0) ; sigBaseline = np.nanstd(baselines, axis=0)
            meanLogKd = np.nanmean(logKd,axis=0) ; sigLogKd  = np.nanstd(logKd,axis=0)
            failDeltaIndices = np.greater( np.fabs(deltas-meanDelta), outlierThreshold*sigDelta)
            failLogKdIndices = np.greater( np.fabs(logKd-meanLogKd), outlierThreshold*sigLogKd)
            nanIndices=np.logical_or(failDeltaIndices,failLogKdIndices)
            print "= = NOTE: %d outliers will be also rejected across %i by %i trials" % ( len(nanIndices[nanIndices]), nCurves, nTrials )
            deltas[nanIndices]=np.nan ; logKd[nanIndices]=np.nan
            meanDelta = np.nanmean(deltas,axis=0) ; sigDelta = np.nanstd( deltas, axis=0)
            #meanBaseline = np.nanmean(baselines,axis=0) ; sigBaseline = np.nanstd(baselines, axis=0)
            meanBaseline = np.nanmean(baselines) ; sigBaseline = np.nanstd(baselines)
            meanLogKd = np.nanmean(logKd,axis=0) ; sigLogKd  = np.nanstd(logKd,axis=0)
            meanKd, sigKd = convert_log_statistics( 10.0, meanLogKd, sigLogKd)

        else:
            print >> sys.stderr, "= = ERROR: invalid error mode selected?"
            sys.exit(1)
        # = = = End of error-mode selection.
    # = = = End of application selection

    # = = = output raw fitting files
    if args.bNMRTitration:
        print "= = Fit complete. kD value is %f" % ( np.fabs(xopt[-1]) )
        sys.exit()

    # = = = Should never go below here now
    print "= = NOTE: Main computation has finished."
    print "    ... conducting quality checks ..."
    if bTestGoodnessOfFit:
        #print "    ....noting exceptions where the goodness-of-fit isitself an outlier."
        meanQoF = np.mean(fitQuality) ; sigQoF = np.std(fitQuality)
        print "= = = NB: Quality-of-fit statistic from all trials: %g +- %g" % ( meanQoF, sigQoF )
        outliers = [ x for x in fitQuality if x-meanQoF>3.0*sigQoF ]
        print "   ...and particular outliers:", outliers

