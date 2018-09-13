import math

'''
Module with methods for doing some statistics.

'''

# Fisher's Exact Test

def getFETpval(a1, a2, b1, b2, left=True):
    """Computes Fisher's exact test based on a
    null-hypothesis distribution specified by the totals, and
    an observed distribution specified by b1 and b2, i.e.
    determines the p-value of b's outcomes 1 and 2.
    The default setting is to use the "left" side of the density
    to determine the p-value.

    Returns p-value."""
    (prob, sless, sright, sleft, slarg)=getFETprob(a1, a2, b1, b2)
    if left:
        return sless # sless
    else:
        return slarg # slarg

def getFET2tail(a1, a2, b1, b2):
    """Computes Fisher's exact test based on a
    null-hypothesis distribution specified by the totals, and
    an observed distribution specified by b1 and b2, i.e.
    determines the two-tailed p-value of b's outcomes 1 and 2.

    Returns p-value."""
    (prob, sless, sright, sleft, slarg)=getFETprob(a1, a2, b1, b2)
    return min(1.0, sleft + sright)

def getFETprob(a1, a2, b1, b2):
    """Computes Fisher's exact test based on a
    null-hypothesis distribution specified by the totals, and
    an observed distribution specified by b1 and b2, i.e.
    determines the probability of b's outcomes 1 and 2.

    Returns an immutable list consisting of the exact
    probability, and assorted p-values (sless, sright, sleft,
    slarg) based on the density."""
    sless  = 0.0
    sright = 0.0
    sleft  = 0.0
    slarg  = 0.0
    n = a1 + a2 + b1 + b2
    row1 = a1 + a2 # the row containing the null hypothesis
    col1 = a1 + b1 # the column containing samples for outcome 1
    max = row1
    if col1 < max:
        max = col1
    min = row1 + col1 - n
    if min < 0:
        min = 0
    if min == max:
        rt = (prob, sless, sright, sleft, slarg) = (1.0,1.0,1.0,1.0,1.0)
        return rt
    prob = hyper0(a1, row1, col1, n)
    sleft = 0.0
    p = hyper(min)

    i = min + 1
    while p < (0.99999999 * prob):
        sleft = sleft + p
        p = hyper(i)
        i = i + 1

    i = i - 1
    if p < (1.00000001 * prob):
        sleft = sleft + p
    else:
        i = i - 1

    sright = 0.0
    p = hyper(max)

    j = max - 1
    while p < (0.99999999 * prob):
        sright = sright + p
        p = hyper(j)
        j = j - 1

    j = j + 1
    if p < (1.00000001 * prob):
        sright = sright + p
    else:
        j = j + 1

    if abs(i - a1) < abs(j - a1):
        sless = sleft
        slarg = 1.0 - sleft + prob
    else:
        sless = 1.0 - sright + prob
        slarg = sright
    return (prob, sless, sright, sleft, slarg)

def lngamm(z):
    # Reference: "Lanczos, C. 'A precision approximation
    # of the gamma function', J. SIAM Numer. Anal., B, 1, 86-96, 1964."
    # Translation of  Alan Miller's FORTRAN-implementation
    # See http://lib.stat.cmu.edu/apstat/245
    x = 0.0
    x = x + 0.1659470187408462e-06/(z+7.0)
    x = x + 0.9934937113930748e-05/(z+6.0)
    x = x - 0.1385710331296526    /(z+5.0)
    x = x + 12.50734324009056     /(z+4.0)
    x = x - 176.6150291498386     /(z+3.0)
    x = x + 771.3234287757674     /(z+2.0)
    x = x - 1259.139216722289     /(z+1.0)
    x = x + 676.5203681218835     /(z)
    x = x + 0.9999999999995183
    return math.log(x)-5.58106146679532777-z+(z-0.5)*math.log(z+6.5)

def lnfact(n):
    if n<=1:
        return 0.0
    return lngamm(n+1.0)

def lnbico(n, k):
    return lnfact(n)-lnfact(k)-lnfact(n-k)

def hyper_323(n11, n1_, n_1, n):
    return math.exp(lnbico(n1_,n11)+lnbico(n-n1_,n_1-n11)-lnbico(n,n_1))

(_sn11, _sn1_, _sn_1, _sn, _sprob) = (0,0,0,0,0.0) # global variables used by hyper0
def hyper0(n11i, n1_i, n_1i, ni):
    global _sn11, _sn1_, _sn_1, _sn, _sprob
    if not ((n1_i|n_1i|ni)!=0):
        if not (n11i % 10 == 0):
            if n11i==_sn11+1:
                _sprob = _sprob * ((_sn1_-_sn11)/float(n11i))*((_sn_1-_sn11)/float(n11i+_sn-_sn1_-_sn_1))
                _sn11 = n11i
                return _sprob
            if n11i==_sn11-1:
                _sprob = _sprob * ((_sn11)/float(_sn1_-n11i))*((_sn11+_sn-_sn1_-_sn_1)/float(_sn_1-n11i))
                _sn11 = n11i
                return _sprob
        _sn11 = n11i
    else:
        _sn11 = n11i
        _sn1_=n1_i
        _sn_1=n_1i
        _sn=ni
    _sprob = hyper_323(_sn11,_sn1_,_sn_1,_sn)
    return _sprob

def hyper(n11):
    return hyper0(n11,0,0,0)

def mean(X):
    sum = 0
    for x in X:
        sum += x
    return sum/len(X)

def meanvar(X):
    """ The mean and variance of the sample. """
    mu = mean(X)
    dev = 0
    for x in X:
        dev += (x - mu) * (x - mu)
    return (mu, dev / len(X))

def getZScore(X, sample):
    (mu, var) = meanvar(X)
    return (sample - mu) / math.sqrt(var)

def getZScores(X):
    (mu, var) = meanvar(X)
    Y = [((x - mu) / math.sqrt(var)) for x in X]
    return Y

def getPearson(X, Y):
    """ Pearson correlation coefficient (r). Note that we are using the standard deviation of the sample, NOT the sample standard deviation (see http://en.wikipedia.org/wiki/Standard_deviation).
    """
    (Xmu, Xvar) = meanvar(X)
    (Ymu, Yvar) = meanvar(Y)
    if len(X) != len(Y):
        raise RuntimeError('Vectors are of uneven length')
    n = len(X)
    sum = 0
    for i in range(n):
        sum += (X[i] * Y[i])
    if n == 0 or Xvar == 0 or Yvar == 0:
        return 0
    return (sum - n * (Xmu * Ymu)) / (n * math.sqrt(Xvar) * math.sqrt(Yvar))

# normal distribution
def error(x):
    """
    Error function
    Cephes Math Library Release 2.8:  June, 2000
    Copyright 1984, 1987, 1988, 1992, 2000 by Stephen L. Moshier
    """
    result = 0.0
    xsq = 0.0
    s = 0.0
    p = 0.0
    q = 0.0

    s = +1
    if x<0:
        s = -1
    x = abs(x)
    if x<0.5:
        xsq = x*x
        p = 0.007547728033418631287834
        p = 0.288805137207594084924010+xsq*p
        p = 14.3383842191748205576712+xsq*p
        p = 38.0140318123903008244444+xsq*p
        p = 3017.82788536507577809226+xsq*p
        p = 7404.07142710151470082064+xsq*p
        p = 80437.3630960840172832162+xsq*p
        q = 0.0
        q = 1.00000000000000000000000+xsq*q
        q = 38.0190713951939403753468+xsq*q
        q = 658.070155459240506326937+xsq*q
        q = 6379.60017324428279487120+xsq*q
        q = 34216.5257924628539769006+xsq*q
        q = 80437.3630960840172826266+xsq*q
        result = s*1.1283791670955125738961589031*x*p/q
        return result
    elif x>=10:
        result = s
        return result
    result = s*(1-errorComplement(x))
    return result

def errorComplement(x):
    """
    Complementary error function
    Cephes Math Library Release 2.8:  June, 2000
    Copyright 1984, 1987, 1988, 1992, 2000 by Stephen L. Moshier
    """
    result = 0.0
    p = 0.0
    q = 0.0

    if x<0.0:
        result = 2.0-errorComplement(-x)
        return result
    elif x<0.5:
        result = 1.0-errorComplement(x)
        return result
    elif x>=10:
        result = 0
        return result
    p = 0.0
    p = 0.5641877825507397413087057563+x*p
    p = 9.675807882987265400604202961+x*p
    p = 77.08161730368428609781633646+x*p
    p = 368.5196154710010637133875746+x*p
    p = 1143.262070703886173606073338+x*p
    p = 2320.439590251635247384768711+x*p
    p = 2898.0293292167655611275846+x*p
    p = 1826.3348842295112592168999+x*p
    q = 1.0
    q = 17.14980943627607849376131193+x*q
    q = 137.1255960500622202878443578+x*q
    q = 661.7361207107653469211984771+x*q
    q = 2094.384367789539593790281779+x*q
    q = 4429.612803883682726711528526+x*q
    q = 6089.5424232724435504633068+x*q
    q = 4958.82756472114071495438422+x*q
    q = 1826.3348842295112595576438+x*q
    result = math.exp(-(x*x))*p/q
    return result

def f(x):
    """
    Normal distribution function
    Returns the area under the Gaussian probability density
    function, integrated from minus infinity to x
    Cephes Math Library Release 2.8:  June, 2000
    Copyright 1984, 1987, 1988, 1992, 2000 by Stephen L. Moshier
    """

    result = 0.0

    result = 0.5*(error(x/1.41421356237309504880)+1)
    return result

def inverseError(e):
    """
    Inverse of the error function
    Cephes Math Library Release 2.8:  June, 2000
    Copyright 1984, 1987, 1988, 1992, 2000 by Stephen L. Moshier
    """
    result = 0.0

    result = inverse(0.5*(e+1))/math.sqrt(2)
    return result

def inverse(y0):
    """
    Inverse of Normal distribution function
    Returns the argument, x, for which the area under the
    Gaussian probability density function (integrated from
    minus infinity to x) is equal to y.

    For small arguments 0 < y < exp(-2), the program computes
    z = sqrt( -2.0 * log(y) );  then the approximation is
    x = z - log(z)/z  - (1/z) P(1/z) / Q(1/z).
    There are two rational functions P/Q, one for 0 < y < exp(-32)
    and the other for y up to exp(-2).  For larger arguments,
    w = y - 0.5, and  x/sqrt(2pi) = w + w**3 R(w**2)/S(w**2)).

    Cephes Math Library Release 2.8:  June, 2000
    Copyright 1984, 1987, 1988, 1992, 2000 by Stephen L. Moshier
    """
    result = 0.0
    expm2 = 0.0
    s2pi = 0.0
    x = 0.0
    y = 0.0
    z = 0.0
    y2 = 0.0
    x0 = 0.0
    x1 = 0.0
    code = 0 # int
    p0 = 0.0
    q0 = 0.0
    p1 = 0.0
    q1 = 0.0
    p2 = 0.0
    q2 = 0.0

    MAX_VALUE = 1.e23
    expm2 = 0.13533528323661269189
    s2pi = 2.50662827463100050242
    if y0<=0:
        result = -MAX_VALUE
        return result
    elif y0>=1:
        result = MAX_VALUE
        return result
    code = 1
    y = y0
    if y>1.0-expm2:
        y = 1.0-y
        code = 0
    if y>expm2:
        y = y-0.5
        y2 = y*y
        p0 = -59.9633501014107895267
        p0 = 98.0010754185999661536+y2*p0
        p0 = -56.6762857469070293439+y2*p0
        p0 = 13.9312609387279679503+y2*p0
        p0 = -1.23916583867381258016+y2*p0
        q0 = 1.0
        q0 = 1.95448858338141759834+y2*q0
        q0 = 4.67627912898881538453+y2*q0
        q0 = 86.3602421390890590575+y2*q0
        q0 = -225.462687854119370527+y2*q0
        q0 = 200.260212380060660359+y2*q0
        q0 = -82.0372256168333339912+y2*q0
        q0 = 15.9056225126211695515+y2*q0
        q0 = -1.18331621121330003142+y2*q0
        x = y+y*y2*p0/q0
        x = x*s2pi
        result = x
        return result
    x = math.sqrt(-(2.0*math.log(y)))
    x0 = x-math.log(x)/x
    z = 1.0/x
    if x<8.0:
        p1 = 4.05544892305962419923
        p1 = 31.5251094599893866154+z*p1
        p1 = 57.1628192246421288162+z*p1
        p1 = 44.0805073893200834700+z*p1
        p1 = 14.6849561928858024014+z*p1
        p1 = 2.18663306850790267539+z*p1
        p1 = -(1.40256079171354495875*0.1)+z*p1
        p1 = -(3.50424626827848203418*0.01)+z*p1
        p1 = -(8.57456785154685413611*0.0001)+z*p1
        q1 = 1.0
        q1 = 15.7799883256466749731+z*q1
        q1 = 45.3907635128879210584+z*q1
        q1 = 41.3172038254672030440+z*q1
        q1 = 15.0425385692907503408+z*q1
        q1 = 2.50464946208309415979+z*q1
        q1 = -(1.42182922854787788574*0.1)+z*q1
        q1 = -(3.80806407691578277194*0.01)+z*q1
        q1 = -(9.33259480895457427372*0.0001)+z*q1
        x1 = z*p1/q1
    else:
        p2 = 3.23774891776946035970
        p2 = 6.91522889068984211695+z*p2
        p2 = 3.93881025292474443415+z*p2
        p2 = 1.33303460815807542389+z*p2
        p2 = 2.01485389549179081538*0.1+z*p2
        p2 = 1.23716634817820021358*0.01+z*p2
        p2 = 3.01581553508235416007*0.0001+z*p2
        p2 = 2.65806974686737550832*0.000001+z*p2
        p2 = 6.23974539184983293730*0.000000001+z*p2
        q2 = 1.0
        q2 = 6.02427039364742014255+z*q2
        q2 = 3.67983563856160859403+z*q2
        q2 = 1.37702099489081330271+z*q2
        q2 = 2.16236993594496635890*0.1+z*q2
        q2 = 1.34204006088543189037*0.01+z*q2
        q2 = 3.28014464682127739104*0.0001+z*q2
        q2 = 2.89247864745380683936*0.000001+z*q2
        q2 = 6.79019408009981274425*0.000000001+z*q2
        x1 = z*p2/q2
    x = x0-x1
    if code!=0:
        x = -x
    result = x
    return result

def getRSpval(a, b):
    """
    Compute the Wilcoxon rank sum test (aka the Mann-Whitney U-test), return the p-value
    The approximation is based on the normal distribution and is reliable
    when sample sets are of size 5 or larger.
    The default is based on the area of the left side of the Gaussian, relative the
    estimated z-value.
    NULL: a==b ONE-SIDED: a<b (default=left), for ONE-SIDED: b<a (right) use 1-returned value.
    For a two-tailed, double the p-value.
    Implemented by Mikael Boden
    """

    # create a new list consisting of the two sample sets that can be sorted
    lst=[]
    for elem in a:
        lst.append([elem, +1, 0])
    for elem in b:
        lst.append([elem, -1, 0])
    # ok sort it
    lst.sort(lambda p, q: cmp(p[0], q[0]))

    # let's go through it and edit each rank
    rank=0
    na=0
    nb=0 # the number of points in each set (A & B)
    same=[] # a dynamic list to keep track of elements with same measurement
    measurement=lst[0][0]
    for row in lst:
        if row[1]==+1: # belongs to class 'a'
            na=na+1
        else:
            nb=nb+1
        if (measurement!=row[0]): # here's an entry that differed from the previous...
            # before moving on to handling the new element we need to sort out the "old" same list
            firstInGroup=rank+1-len(same)
            lastInGroup=rank
            average=float(lastInGroup-firstInGroup)/2.0
            for srow in same:
                srow[2]=firstInGroup+average
            same=[]
            measurement=row[0]
        same.append(row)
        rank=rank+1
    # the last batch of entries is handled outside the loop...
    firstInGroup=rank+1-len(same)
    lastInGroup=rank
    average=float(lastInGroup-firstInGroup)/2.0
    for srow in same:
        srow[2]=firstInGroup+average

    n=na+nb      # the total number of measurements
    ta_obs=0     # sum of na ranks in group A
    tb_obs=0     # sum of nb ranks in group B
    # sum the ranks (replace the measurements)
    for entry in lst:
        if entry[1]==+1: # class 'a'
            ta_obs+=entry[2]
        else:
            tb_obs+=entry[2]

    tab=ta_obs+tb_obs                     # sum of n ranks in groups A and B combined
    sd=math.sqrt((na*nb*(n+1.0))/12.0)    # the standard deviation is the same in both sets
    ta_null=na*(n+1.0)/2.0                # the sum of the "null" case
    tb_null=nb*(n+1.0)/2.0                # the sum of the "null" case
    ta_max=na*nb+(na*(na+1.0))/2.0        # the max sum set A can take
    tb_max=na*nb+(nb*(nb+1.0))/2.0        # the max sum set B can take
    ua=ta_max-ta_obs                      # the "U" value for A which is the mirror of ...
    ub=tb_max-tb_obs                      # the "U" value for B (we only need one)
    ua_null=ta_max-ta_null                # the U value for the null case
    ub_null=tb_max-tb_null
    if ta_obs>ta_null:                    # a "continuity correction" for A
        da=-0.5
    else:
        da=+0.5
    if tb_obs>tb_null:                    # a "continuity correction" for B
        db=-0.5
    else:
        db=+0.5
    za=((ta_obs-ta_null)+da)/sd           # the z value for A which is the mirror of ...
    zb=((tb_obs-tb_null)+db)/sd           # the z value for B (we only need one)
    p=f(za)                        # figure out the area of the normal distribution
    u=ua;                                 # remember one of the U values
    return p                              # the p-value: null is that a==b, one-sided (a has lower values)


def getPointBiserialCorr(group1, group2):
    """
    The point biserial correlation coefficient (rpb) is a correlation coefficient used when one variable (e.g. Y) is dichotomous,
    with continuous data divided into two groups (group1 and group2 here).
    group1 corresponds to "greater", group to "lesser", i.e. 1 and 0 respectively
    See https://en.wikipedia.org/wiki/Point-biserial_correlation_coefficient
    """
    n1 = len(group1)
    n0 = len(group2)
    if n1 < 1 or n0 < 1:
        raise RuntimeError('At least one group is empty')
    n  = n1 + n0
    M1 = sum(group1) / float(n1)
    M0 = sum(group2) / float(n0)
    M  = (M1 * n1 + M0 * n0) / float(n)
    all = []
    all.extend(group1)
    all.extend(group2)
    sn = math.sqrt(sum([(x_i - M)**2 for x_i in all]) / float(n))
    return (M1 - M0) / sn * math.sqrt((n1 * n0) / float(n**2))