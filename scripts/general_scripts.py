import numpy as np
import sys
import csv
import re
# This function returns a string that contains
# the value of the function with its error formatted
# to the same exponent.
# e.g. f(0.00302, 4.52e-5, 2) = 3.020e-5 +- 0.045e-5
# The precision argument N sets the number of decimal points
# such at most N+1 sig. figs. are in both values

def normalise_vector_array(v):
    """
    Takes an ND-array of vectors of shape (i, j, k, ..., 3)
    and noramlises it along the last axis.
    """
    return v/np.sqrt((v**2).sum(-1))[..., np.newaxis]

def format_float_with_error(val, err, prec=4):
    # NB: int floors anyway.
    exp_val = np.floor(np.log10(np.abs(val)))
    exp_err = np.floor(np.log10(np.abs(err)))

    exp_out = max(exp_val,exp_err)
    sig_out = prec
    #sig_out = prec+int(np.abs(exp_val-exp_err))
    #r = " %.*e%i +- %.*e%i " % ( sig_val, val*10**(exp_out-exp_val), exp_out, sig_err, )
    return "%.*fe%i +- %.*fe%i" % (sig_out, val*10**(-exp_out), exp_out, sig_out, err*10**(-exp_out), exp_out)

def load_matrix(fn):
    x=[]
    bFirst=True
    for l in open(fn):
        if l[0]=="#" or l[0]=="@" or l[0]=="&" or l=="" or l[0]=="\n":
            continue
        line = l.split()
        x.append([float(e) for e in line])
        if bFirst:
            xl=len(line)
            bFirst=False
        else:
            if len(line) != xl:
                print( "= = = WARNING: file read in general_scripts.load_matrix() is not rectangular!", file=sys.stderr )
                print( "= = = ...detected a different row length to the first row read: %d, %d" % (xl, len(line)), file=sys.stderr )

    return np.array(x)

def load_xy(fn):
    x=[]
    y=[]
    bWarned=False
    for l in open(fn):
        if l[0]=="#" or l[0]=="@" or l[0]=="&" or l=="" or l[0]=="\n":
            continue
        lines = l.split()
        #if lines[0].isalnum():
        if re.search('[a-df-zA-DF-Z]', lines[0]):
            if not bWarned:
                bWarned=True
                print( "= = WARNING: Input file %s to general_scripts.load_xy contains uncommented alphanumerics: %s" % (fn, lines[0]), file=sys.stderr )
            continue
        x.append(float(lines[0]))
        y.append(float(lines[1]))
    return np.array(x), np.array(y)

def load_xys(fn):
    x=[]
    y=[]
    bWarned=False
    for l in open(fn):
        if l[0]=="#" or l[0]=="@" or l[0]=="&" or l=="":
            continue
        lines=l.split()
        #if lines[0].isalnum():
        if re.search('[a-df-zA-DF-Z]', lines[0]):
            if not bWarned:
                bWarned=True
                print( "= = WARNING: Input file %s to general_scripts.load_xys contains uncommented alphanumerics: %s" % (fn, lines[0]), file=sys.stderr )
            continue
        vars = [ float(i) for i in lines ]
        x.append(vars[0])
        y.append(vars[1:])
    return np.array(x), np.array(y)

def load_xydy(fn):
    x=[]
    y=[]
    dy=[]
    bWarned=False
    for l in open(fn):
        if l[0]=="#" or l[0]=="@" or l[0]=="&" or l=="" or l=="\n":
            continue
        lines = l.split()
        #if lines[0].isalnum():
        if re.search('[a-df-zA-DF-Z]', lines[0]):
            if not bWarned:
                print( "= = WARNING: Input file %s to general_scripts.load_xydy contains uncommented alphanumeric text: %s" % (fn, lines[0]), file=sys.stderr )
                bWarned=True
            continue
        x.append(float(lines[0]))
        y.append(float(lines[1]))
        if len(lines) <=2:
            print( '= = ERROR: in general_scripts.load_xydy, input file does not have a third column as dy! Aborting..', file=sys.stderr )
            print( '- - line content: '+l, file=sys.stderr )
            sys.exit(1)
        dy.append(float(lines[2]))
    return np.array(x), np.array(y), np.array(dy)

def load_block_as_numpy(fn, ignores='#@', newblock='&', bVerbose=False ):
    """
    load_numpyblock transforms a freeform input file into an numpy array,
    without explicitly checking that all numbers are well-formated.
    1) It expects to ignore any lines starting with a comment character from the given set
    (by default the xmgrace '#' for comments and '@' for commands). Any empty lines
    will by default be also ignored for safety purposes.
    2) It will take the character '&' to denote termination of a 2D-block.

    Special arguments:
    - For newline delimited 3D arrays such as gnuplot, please give the empty string '' to the newblock argument.
    - For files where text are not commented such as synchroton data files, include the string 'alpha' in ingores,
      this will add [a-zA-Z]
    """
    out3D=[] ; out2D=[]
    bEmptyAsNewblock = ( newblock == '' )
    bAlpha = ( 'alpha' in ignores )
    if bAlpha:
        ignores = ignores.replace('alpha','')

    if bVerbose:
        print( "= = load_block_as_numpy is reading file %s ..." % fn )
    for l in open(fn):
        if any( l[0]==char for char in ignores ):
            continue
        if bAlpha and l[0].isalpha():
            continue
        if l[0]=='':
            if bEmptyAsNewblock:
                out3D.append( out2D ) ; out2D = []
            continue
        if any( l[0]==char for char in newblock ):
            # Empty newblocks autofail this condition, AFAIK.
            out3D.append( out2D ) ; out2D = []
            continue
        # Assume to be data line.
        lines = l.split()
        out2D.append( [ float(i) for i in lines ] )

    if bVerbose:
        print( "= = ... read finished." )
    # Now determine dimensionality. Three are three cases.
    if out3D == []:
        # Basic 2D without "&" character.
        return np.array( out2D )
    else:
        # Potential 3D data.
        if out2D != []:
            # Unclean input file where the final 2D-set does no have the required "&"
            out3D.append( out2D ) ; out2D = []
        if len(out3D) == 1:
            # Basic 2D set but terminated with "&" character
            return np.array( out3D[0] )
        else:
            # Full 3D set
            return np.array( out3D )
    # Impossible section.
    return -1

def load_xylist(fn):
    xlist=[] ; ylist=[]
    x=[] ; y=[]
    for l in open(fn):
        if l[0]=="#" or l[0]=="@" or l=="":
            continue
        if l[0]=="&":
            xlist.append(x) ; ylist.append(y)
            x=[] ; y=[]
            continue
        lines = l.split()
        x.append(float(lines[0])) ; y.append(float(lines[1]))
    if x != []:
        xlist.append(x) ; ylist.append(y)

    return xlist, ylist

def load_xydylist(fn):
    xlist=[] ; ylist=[] ; dylist=[]
    x=[] ; y=[] ; dy=[]
    for l in open(fn):
        if l[0]=="#" or l[0]=="@" or l=="":
            continue
        if l[0]=="&":
            xlist.append(x) ; ylist.append(y) ; dylist.append(dy)
            x=[] ; y=[] ; dy=[]
            continue
        lines = l.split()
        x.append(float(lines[0]))
        y.append(float(lines[1]))
        dy.append(float(lines[2]))

    if x != []:
        xlist.append(x) ; ylist.append(y) ; dylist.append(dy)

    return xlist, ylist, dylist

def load_sxydylist(fn, key="legend"):
    leglist=[]
    xlist=[]
    ylist=[]
    dylist=[]
    x=[] ; y=[] ; dy=[]
    for l in open(fn):
        lines = l.split()
        if l=="" or l=="\n":
            continue
        if l[0]=="#" or l[0]=="@":
            if key in l:
                leglist.append(lines[-1].strip('"'))
            continue
        if l[0]=="&":
            xlist.append(x) ; ylist.append(y)
            if len(dy)>0:
                dylist.append(dy)
            x=[] ; y=[] ; dy=[]
            continue
        x.append(float(lines[0]))
        y.append(float(lines[1]))
        if len(lines) > 2:
            dy.append(float(lines[2]))

    if x != []:
        xlist.append(x) ; ylist.append(y) ; dylist.append(dy)

    if dylist != []:
        return leglist, np.array(xlist), np.array(ylist), np.array(dylist)
    else:
        return leglist, np.array(xlist), np.array(ylist), []


# = = = Printing functions = = =

def xmgrace_template( plotType ):
    header="@type %s" % plotType
    footer="&"
    return header, footer

def _print_flexible(filePointer, variable):
    """
    Internal handler for taking several kinds of inputs of var, according to the type of variable.
    Giving the empty string will do nothing.
    """
    if variable=='':
        return
    t=type(variable)
    if t==type('a') or t==type(1):
        print( variable, file=filePointer )
    elif t==type([]) or t==type(()):
        for line in variable:
            print( line, file=filePointer )
    else:
        print( '= = ERROR in _print_flexible: type %s is not accounted for in the function.' % t, file=sys.stderr )
        sys.exit(1)

def print_xy(fn, x, y, dy=[], header="", footer=""):
    fp=open(fn,'w')
    _print_flexible(fp, header)
    if dy==[]:
        for i in range(len(x)):
            print( x[i], y[i], file=fp )
    else:
        for i in range(len(x)):
            print( x[i], y[i], dy[i], file=fp )
    _print_flexible(fp, footer)
    fp.close()

def print_xydy(fn, x, y, dy, header="", footer=""):
    print_xy(fn, x, y, dy, header, footer)

def print_xylist(fn, x, ylist, bCols=False, header=""):
    """
    Array formats: x(nvals) y(nplots,nvals)
    bCols will stack all y contents in the same line, useful for errors.
    """
    fp = open( fn, 'w')
    _print_flexible(fp, header)

    ylist=np.array(ylist)
    shape=ylist.shape
    print( shape )
    if len(shape)==1:
        for j in range(len(x)):
            print( x[j], ylist[j], file=fp )
        print( "&", file=fp )
    elif len(shape)==2:
        nplot=shape[0]
        nvals=shape[1]
        if bCols:
            for j in range(nvals):
                s="%g " % x[j] + " ".join("%g" % ylist[i][j] for i in range(nplot))
                print( s, file=fp )
            print( "&", file=fp )
        else:
            for i in range(nplot):
                for j in range(len(x)):
                    print( x[j], ylist[i][j], file=fp )
                print( "&", file=fp )
    fp.close()

def print_sxylist(fn, legend, x, ylist, header=[]):
    fp = open( fn, 'w')
    for line in header:
        print( "%s" % line, file=fp )

    ylist=np.array(ylist)
    shape=ylist.shape
    nplot=len(ylist)
    s=0
    for i in range(nplot):
        print( "@s%d legend \"%s\"" % (s, legend[i]), file=fp )
        for j in range(len(x)):
            print( x[j], str(ylist[i][j]).strip('[]'), file=fp )
        print( "&", file=fp )
        s+=1
    fp.close()

def print_s3d(fn, legend, arr, cols, header=[]):
    fp = open( fn, 'w')
    for line in header:
        print( "%s" % line, file=fp )

    shape=arr.shape
    ncols=len(cols)
    s=0
    for i in range(shape[0]):
        print( "@s%d legend \"%s\"" % (s, legend[i]), file=fp )
        for j in range(shape[1]):
            s=" ".join("%g" % arr[i,j,cols[k]] for k in range(ncols))
            print( s, file=fp )
        print( "&", file=fp )
        s+=1
    fp.close()

def print_R_hist(fn, hist, edges, header=''):
    fp = open(fn, 'w')
    nbins=hist.shape
    dim=len(nbins)

    if header != '':
        print( '%s' % header, file=fp )
    print( '# DIMENSIONS: %i' % dim, file=fp )
    s="# BINWIDTH: " + " ".join("%g" % ((edges[i][-1]-edges[i][0])/nbins[i]) for i in range(dim) )
    print( s, file=fp )
    s="# NBINS: " + " ".join("%g" % (nbins[i]) for i in range(dim) )
    print( s, file=fp )
    for index, val in np.ndenumerate(hist):
        s=" ".join("%g %g" % (edges[i][index[i]], edges[i][index[i]+1]) for i in range(dim) )
        s=s+" %g" % val
        print(s, file=fp)

def print_gplot_hist(fn, hist, edges, header='', bSphere=False):

    fp = open(fn, 'w')
    nbins=hist.shape
    dim=len(nbins)

    # For Gnuplot, we should plot the average between the two edges.
    # ndemunerate method iterates across all specified dimensions
    # First gather and print( some data as comments to aid interpretation and future scripting. )
    # Use final histogram output as authoritative source!
    if header != '':
        print( '%s' % header, file=fp )
    print( '# DIMENSIONS: %i' % dim, file=fp )
    s="# BINWIDTH: " + " ".join("%g" % ((edges[i][-1]-edges[i][0])/nbins[i]) for i in range(dim) )
    print( s, file=fp )
    s="# NBINS: " + " ".join("%g" % (nbins[i]) for i in range(dim) )
    print( s, file=fp )

    if bSphere:
        if dim != 2:
            print( "= = = ERROR: histogram data is not in 2D, but spherical histogram plotting is requested!", file=sys.stderr )
            sys.exit(1)
        # Handle spherical data by assuming that X is wrapped, and y is extended.
        # Can assume data is only 2D.
        xmin=0.5*(edges[0][0]+edges[0][1])
        ymin=edges[1][0]
        ymax=edges[1][-1]
        for eX in range(nbins[0]):
            xavg=0.5*(edges[0][eX]+edges[0][eX+1])
            # Print polar-caps to complete sphere
            print( '%g %g %g' % (xavg, ymin, hist[eX][0]), file=fp )
            for eY in range(nbins[1]):
                yavg=0.5*(edges[1][eY]+edges[1][eY+1])
                print( '%g %g %g' % (xavg, yavg, hist[eX][eY]), file=fp )
            print( '%g %g %g' % (xavg, ymax, hist[eX][-1]), file=fp )
            print( '', file=fp )
        # Print first line again to complete sphere, with 2-pi deviation just in case.
        print( '%g %g %g' % (xmin+2*np.pi, ymin, hist[0][0]), file=fp )
        for eY in range(nbins[1]):
            yavg=0.5*(edges[1][eY]+edges[1][eY+1])
            print( '%g %g %g' % (xmin+2*np.pi, yavg, hist[0][eY]), file=fp )
        print( '%g %g %g' % (xmin+2*np.pi, ymax, hist[0][-1]), file=fp )
        print( '', file=fp )
    else:
        for index, val in np.ndenumerate(hist):
            s=" ".join("%g" % (0.5*(edges[i][index[i]]+edges[i][index[i]+1])) for i in range(dim) )
            s=s+" %g" % val
            print(s, file=fp)
            if index[-1] == nbins[-1]-1:
                print( '', file=fp )
    fp.close()

def print_gplot_4d(fn, datablock, x, y, z, header=''):
    sh = datablock.shape
    dims = len(sh)
    if dims != 3:
        print( " = = ERROR: general_scripts.print_gplot_3d requires the data to be in 3D.", file=sys.stderr )
        sys.exit(1)

    fp = open(fn, 'w')
    if header != '':
        print( header, file=fp )

    for i in range(sh[0]):
        for j in range(sh[1]):
            for k in range(sh[2]):
                print( "%g %g %g %g" % (x[i], y[j], z[k], datablock[i,j,k]), file=fp )

    fp.close()

def print_numpy_block(fn, data, header='', delim='&', axis=-1):
    """
    Unformatted output for numpy arrays in 2D or 3D.
    axis == -1 means each line spans the last  axis of data.
    axis ==  0 means each line spans the first axis of data.
    """
    if axis != 0 and axis != -1:
        print( "= = ERROR in print_numpy_block, axis argument must be either 0 or -1 !", file=sys.stderr )
        return -1
    shape = data.shape
    dims= len(shape)
    if dims > 3 :
        print( "= = ERROR in print_numpy_block, cannot yet deal with 4+ dimensions in numpy array!", file=sys.stderr )
        return -1

    fp = open(fn, 'w')
    writer = csv.writer(fp, delimiter=' ')
    if header != '':
        print( header, file=fp )
    if dims == 2:
        if axis == -1:
            for i in range(shape[0]):
                writer.writerow(data[i])
        elif axis == 0:
            for i in range(shape[1]):
                writer.writerow(data[:,i])
    elif dims == 3:
        if axis == -1:
            for i in range(shape[0]):
                for j in range(shape[1]):
                    s=" ".join('%g ' % data[i,j,k] for k in range(shape[2]))
                    print(s, file=fp)
                    print('', file=fp)
                    #print( '\n', file=fp )
                print( delim, file=fp )
        elif axis == 0:
            for i in range(shape[1]):
                for j in range(shape[2]):
                    s=" ".join('%g ' % data[k,i,j] for k in range(shape[0]))
                    print(s, file=fp)
                    print('', file=fp)
                    #print( '\n', file=fp )
                print( delim, file=fp )
    # Done.
    fp.close()

def format_header_legend(legends, s_init=0, step=1):
    string=''
    nLegs=len(legends)
    s=s_init
    for i in range(nLegs):
        string=string+'@s%i legend "%s"\n' % (s, legends[i])
        s+=step
    return string

def print_xvg_simple(fn, x, y, dy=[]):
    if dy==[]:
        h="@type xy"
    else:
        h="@type xydy"
    f="&"
    print_xy(fn, x, y, dy=dy, header=h, footer=f)

def print_xvg_matrix(fn, db, header='', bReverse=False, bAmpersand=True):
    fp = open(fn, 'w')
    _print_flexible(fp, header)
    if isinstance( db, list):
        for i in len(db):
            sh = db[i].shape
            if len(sh)==2:
                if bReverse:
                    for j in range(sh[-1]):
                        s=" ".join('%f' % db[i][k,j] for k in range(sh[0]))
                        print(s, file=fp)
                else:
                    for j in range(sh[0]):
                        s=" ".join('%f' % db[i][j,k] for k in range(sh[-1]))
                        print(s, file=fp)
            if bAmpersand:
                print( '&', file=fp )
    elif isinstance( db, np.ndarray):
        sh = db.shape
        if len(sh)==2:
            if bReverse:
                for j in range(sh[-1]):
                    s=" ".join('%f' % db[k,j] for k in range(sh[0]))
                    print(s, file=fp)
            else:
                for j in range(sh[0]):
                    s=" ".join('%f' % db[j,k] for k in range(sh[-1]))
                    print(s, file=fp)
        if bAmpersand:
            print( '&', file=fp )

    fp.close()

def print_xvg_complex(fn, x, dataBlock, header='', bType=False, bTarget=False, graph=0):
    fp = open(fn, 'w')
    sh = dataBlock.shape
    nVals = len(x)
    _print_flexible(fp, header)
    # 3D - multiple plots sharing the same X
    if len(sh) == 3:
        if bType:
            if sh[1] == 2:
                print( '@type xydy', file=fp )
            elif sh[1] == 1:
                print( '@type xy', file=fp )
        series=0
        # Each axis-0 in dataBlock represents a graph that must be plotted.
        for i in range(sh[0]):
            if bTarget:
                print( '@target g%i.s%i' % (graph, series), file=fp )
            for h in range(nVals):
                s='%g ' % x[h] + " ".join('%g' % dataBlock[i,j,h] for j in range(sh[1]))
                print(s, file=fp)
            print( '&', file=fp )
            series+=1
    # 2D - a single plot with a separate X.
    elif len(sh) == 2:
        if bType:
            if sh[0] == 2:
                print( '@type xydy', file=fp )
            elif sh[0] == 1:
                print( '@type xy', file=fp )
        if bTarget:
            print( '@target g%i.s%i' % (graph, series), file=fp )
        for h in range(nVals):
            s='%g ' % x[h] + " ".join('%g' % dataBlock[j,h] for j in range(sh[0]))
            print(s, file=fp)
        print( '&', file=fp )
    # 1D - single plot again with separate X.
    elif len(sh) == 1:
        if bType:
            print( '@type xy', file=fp )
        for h in range(nVals):
            print( '%g %g' % (x[h], dataBlock[h]), file=fp )
        print( '&', file=fp )

    fp.close()
