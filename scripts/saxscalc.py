import numpy as np
import numpy.ma as ma
from scipy.stats import sem, gmean
import math
import random
import sys

def block_average(x, stride, bSEM=False):
    """
    Returns a block-wise average by filling out empty values with np.nan.
    This is used to compute data as blocks in stretches of Delta q = pi / D_max,
    as cited in V_R.
    """
    if len(x) % stride != 0:
        t=np.pad( x, (0, stride - x.size % stride), mode='constant', constant_values=np.NaN).reshape(-1,stride)
    else:
        t=x.reshape(-1,stride)

    if bSEM:
        s  = np.nanmean(t, axis=1)
        ds = sem(t, axis=1, nan_policy='omit').data
        return np.stack((s,ds)).T
    else:
        return np.nanmean(t, axis=1)

def block_gmean(x, stride):
    """
    Returns a block geometric mean for an array of length X.
    Fast if and only if the stride divides the length of X, as gmean does not support np.nan filling.
    """
    l=len(x)
    rem=l%stride ; nBins = l/stride
    if rem==0:
        return gmean(x.reshape(nBins,stride),axis=1)
    else:
        tmp=np.array(x)
        out=np.zeros(nBins+1, dtype=tmp.dtype)
        out[:-1]=gmean(tmp[:-rem].reshape(l/stride,stride),axis=1)
        out[-1]=gmean(tmp[nBins*stride:])
        return out

def block_sum(x, stride):
    if len(x) % stride != 0:
        t=np.pad( x, (0, stride - x.size % stride), mode='constant', constant_values=np.NaN).reshape(-1,stride)
    else:
        t=x.reshape(-1,stride)
    return np.nansum(t, axis=1)

def block_weights(x, stride):
    """
    Reduces the contribution of the last element when there are fewer elements in it to represent its lower contributions.
    """
    l=len(x)
    rem=l%stride
    if rem == 0:
        return np.repeat(1.0,l/stride)
    else:
        ret = np.repeat(1.0,1+l/stride)
        ret[-1] = 1.0*rem/stride
        return ret

def random_cointoss(n):
    return [ (random.randint(0,1)>0) for i in range(n) ]

def run_distribution(boolList, bPrune=True ):
    """
    For a list containing True and False, return the distribution of "runs of equal outcome".
    As an example, the list TTFTTTFFFFTTFFFTF converts to [3,2,3,1,0,0,...]
    The output array length will by default be pruned to be as small as possible,
    so as to save some reporting space. This can be turned off using bPrune=False
    so that he reported array is of the same length as the input.
    """
    l=len(boolList)
    out=np.zeros(l)
    c=0
    for x in boolList:
        if c == 0:
            last=x
        elif x != last:
            last=x
            out[c-1]+=1
            c=0
        c+=1
    if c>0:
        out[c-1]+=1
    if bPrune:
        lastElem=np.argwhere(out>0)[-1][0]
        return out[0:lastElem+1]
    else:
        return out

def schilling_equation(n, c):
    """
    A_n(C) as defined by Schilling, Coll. Math. J., 1990 and Franke et al., Nat. Methods, 2015.
    Won't work above n=30 or so due to recursive limit.
    """
    return sum( [schilling_equation(n-1-j, c) for j in range(0,c+1)]) if n>c else 2**n

def schilling_iterative(n, c):
    """
    A_n(C) reformulated to an iterative and not recursive definition.
    The trick is to define an array of size c+1 and replace elements
    as we climb the ladder from c+1 to n.
    Returns a long-int to preserve the actual fidelity.
    """
    if c<0:
        return 0
    elif c>=n:
        return 2**n
    else:
        # float will be required to deal with integer overflow?
        #arr=np.array([ 2**i for i in range(0,c+1)]).astype(float)
        arr=[ 2**i for i in range(0,c+1) ]
        for i in range(0,n-c-1):
            tmp=sum(arr)
            arr[i%(c+1)]=tmp
        return long( sum(arr) )

def probability_cormap_heads(n, c):
    """
    P(R_n>C) as defined by Franke et al., Nat. Methods, 2015.
    Essentially, "What's the likelihood of more than C heads in a row over n tosses?"
    Can have issues calculating the ratio of two very long integers, therefore hack as strings for magnitude comparison.
    Keep at least 16-digits for magnitude comparison.
    """
    if n < 1024:
        return 1.0-1.0*schilling_iterative(n, c)/2**n
    else:
        numer=str(schilling_iterative(n, c))
        denom=str(2**n)
        lenDEL=min(len(numer),len(denom))-16
        return 1.0-float(numer[:-lenDEL])/float(denom[:-lenDEL])

def probability_cormap_either(n, c):
    """
    As heads but with slight shift. Schilling writes B_n(C) = 2*A_{n-1}(C-1).
    """
    if n < 1024:
        return 1.0-2.0*schilling_iterative(n-1, c-1)/2**n
    else:
        numer=str(schilling_iterative(n-1, c-1))
        denom=str(2**n)
        lenDEL=min(len(numer),len(denom))-16
        return 1.0-2.0*float(numer[:-lenDEL])/float(denom[:-lenDEL])

def cormap_value(x1, x2):
    """
    Implementation of Franke, Nat. Methods, 2015.
    Reports on the agreement of the two curves using the deviations as comparison.
    Looks at the statistical distribution of the boolean ( I_1(q) > I_2(q) ) over the range of q.
    For two experimental data sets of the same measurement, this distribution is expected
    to be like a random coin toss. So, if significant stretches of same values occur, then
    it's unlikely to be equal.
    This agreement is reported as a significant level alpha.
    """
    numPoints=len(x1)
    runs = run_distribution( (x1>x2), bPrune=True )
    maxRun = len(runs)
    print >> sys.stderr, "= = Longest sequence found %i in %i tosses." % (maxRun, numPoints)
    probExceed = probability_cormap_either(numPoints, maxRun)
    print >> sys.stderr, "= = Probability of a single random-sequence exceeding this: %g" % probExceed
    return probExceed

def cormap_matrix(x1, x2):
    """
    Plot the CorMap visualisation to see the checkerboard pattern.
    """
    boolVec=(x1>x2)
    return boolVec.reshape(len(boolVec),1) ^ boolVec

def print_cormap_matrix( fp, mat ):
    sh = mat.shape
    for i in range(sh[0]):
        for j in range(sh[1]):
            print >> fp, int(mat[i,j]),
        print >> fp, ''
#        print >> fp, ' '.join(str(mat[i].astype(int))).strip('[]')


def _dummy():
    """
    Notes for volatility ratio uncertainties:
    dR/R = sqrt( (dIa/Ia)^2 + (dIb/Ib)^2 ) - at an individual point q_i

    geometric mean: dR/R 

    then,
    dVR^2 = Sum_i=1^N (Rj*dRi)^2 + (Ri*dRj)^2 ) / (Ri+Rj)^2 

    """

def volatility_ratio_scaling(x1, x2):
    return np.mean( np.divide(x1,x2) )

def volatility_ratio(x1, x2, stride=1, bElementwise=False, bKeepPartialBin=False, bReweightPartialBin=True):
    """
    Implementation of Hura et. al, Nat. Methods, 2013.
          N_bins-1 |  R(q_i) - R(q_{i+1})    |
    V_R =    Sum   |  -------------------    |
             i=1   | (R(q_i) + R(q_{i+1}))/2 |
    Where R is the geometric-mean value of the ratios found within a bin.
    Partial bins are simply removed unless an argument is made to preserve them,
    whereupon they are down-weighted so as to reduce their impact on the total ratio.
    e.g., the last bin might be only 20% of the size of a full bin. Its contribution is then multipled by 0.2.
    Note that negative values are masked out to prevent instability.
    """
    if np.any(x1<=0) or np.any(x2<=0):
#        print >> sys.stderr, "= = = WARNING: negative values encountered in the volatility ratio calculations. " \
#                "To avoid NaN returns, these values will be masked - which will reduce the speed of calculations."
        x1=np.ma.masked_array(x1,mask=(x1<=0))
        x2=np.ma.masked_array(x2,mask=(x2<=0))
#                "To avoid NaN returns, this will be adjusted to a floor value of %g" % floor
#        x1=np.max(np.vstack((x1,np.repeat(floor,len(x1)))),axis=0)
#        x2=np.max(np.vstack((x2,np.repeat(floor,len(x2)))),axis=0)

    ratio=np.divide(x1,x2)
    #ratio/=np.mean(ratio)
    ratio/=gmean(ratio)
    if stride>1:
        rem   = len(x1)%stride ; nBins = len(x1)/stride
        if not bKeepPartialBin and rem>0:
            ratio = block_gmean(ratio[:-rem],stride)
        else:
            ratio = block_gmean(ratio,stride)
    # = = If we're keeping the partial bin, the length of ratio will be 1 longer. = =
    V_R = 2.0*np.fabs( (ratio[:-1]-ratio[1:])/(ratio[:-1]+ratio[1:]) )

    if bReweightPartialBin and stride>1 and rem > 0 :
        V_R[-1] *= 1.0*(len(x1)%stride)/stride
        if bElementwise:
            return V_R
        else:
            return np.sum(V_R)
    else:
        if bElementwise:
            return V_R
        else:
            return np.sum(V_R)

def log_chi_square(x1, x2, dx1=[], dx2=[], stride=1, bElementwise=False):
    if dx1 == [] and dx2 == []:
        chi2 = np.square(np.log(x1)-np.log(x2))
    elif dx1 != [] and dx2 != []:
        chi2 = np.divide( np.square(np.log(x1)-np.log(x2)), dx1*dx1/x1/x1+dx2*dx2/x2/x2 )
    else:
        print >> sys.stderr, "= = ERROR: the uncertainties of both x1 and x2 must be either given or not. Quitting."
        sys.exit(1)
    chi2/=len(chi2)
    if stride>1:
        chi2 = block_sum(chi2, stride)
    if bElementwise:
        return chi2*len(chi2)
    else:
        return np.sum(chi2)

def chi_square(x1, x2, dx1=[], dx2=[], stride=1, bElementwise=False):
    if dx1 == [] and dx2 == []:
        chi2 = np.square( x1 - x2 )
    elif dx1 != [] and dx2 != []:
        chi2 = np.divide( np.square( (x1-x2) ), dx1*dx1+dx2*dx2 )
    else:
        print >> sys.stderr, "= = ERROR: the uncertainties of both x1 and x2 must be either given or not. Quitting."
        sys.exit(1)
    chi2/=len(chi2)
    if stride>1:
        chi2 = block_sum(chi2, stride)
    if bElementwise:
        return chi2*len(chi2)
    else:
        return np.sum(chi2)

def chi_square_free(x1, x2, stride, nParams=0, dx1=[], dx2=[], nRounds=500, bElementwise=False):
    """
    Implementation of Rambo and Tainer, Nature, 2013.

    Partial bins are down-weighted so as to reduce their impact on the total chi value.
    e.g., the last bin might be only 20% of the size of a full bin. Its contribution is then multipled by 0.2.

    The number of degrees of freedom is decided by the number of bins:
    nDOF ~= len(x1) / stride - nParams
    nParams is the number of free parameters used that reduce these effective numbers of freedoms.
    The final weight is a float value if there is a partial bin at the end.
    """
    chi2=chi_square(x1,x2,dx1,dx2,bElementwise=True)
    l=len(chi2) ; rem=l%stride
    nChunks=l/stride
    if rem==0:
       nDOF=1.0*(nChunks-nParams)
    else:
       nDOF=nChunks-nParams+1.0*rem/stride
    chilist = []
    for k in range(nRounds):
        val=0
        for i in range(nChunks):
            j=random.randint(0,stride-1)
            val += chi2[i*stride+j]
        if rem != 0:
            j=random.randint(0,rem-1)
            val += chi2[nChunks*stride+j]*1.0*rem/stride
            #val += chi2[nChunks*stride+j]
        val /= nDOF
        chilist.append(val)
    list.sort(chilist)
    if bElementwise:
        return chilist
    else:
        return chilist[nRounds/2]
