from math import log, exp
import sys


MAXIT = 100
EPS = 3.0e-7
FPMIN = 1.0e-300
gamma_c = [76.18009172947146,
           -86.50532032941677,
           24.01409824083091,
           -1.23173957245,
           0.1208650973866179e-2,
           -0.5395239384953e-5]


def log_binomial_ncdf(N, k, p):
    """
    Log of one minus the cumulative distribution function of the binomial dist.

    The binomial density gives the probability of k successes in N independent
    trials each with probability p of success.
    """
    if (k==0):
        return 0
    else:
        return log_betai(k, N-k+1, p)


def betai (a, b, x):
    """
    Incomplete beta function
    """
    if (x<0 or x>1): die("Bad x=`" + str(x) + "' in routine betai")
    if (x==0 or x==1):
        bt = 0
    else:
        bt = exp(gammaln(a+b)-gammaln(a)-gammaln(b)+a*log(x)+b*log(1-x))

    thresh = (a+1)/(a+b+2.0)
    if (x<thresh):
        return(bt*betacf(a,b,x)/a)
    else:
        return(1.0-bt*betacf(b,a,1.0-x)/b)


def log_betai(a, b, x):
    """
    log incomplete beta function
    """
    if (x<0 or x>1): die("Bad x=`" + str(x) + "' in routine betai")
    if (x==0 or x==1):
        log_bt = -1e300           # log(0)
    else:
        log_bt = gammaln(a+b)-gammaln(a)-gammaln(b)+a*log(x)+b*log(1.0-x)

    thresh = (a+1.0)/(a+b+2.0)
    if (x<thresh):
        return(log_bt + log(betacf(a,b,x)/a))
    else:
        return(log(1.0 - exp(log_bt)*betacf(b,a,1.0-x)/b))


def betacf(a, b, x):
    """
    used by betai
    """
    qab = a+b
    qap = a+1.0
    qam = a-1.0
    c = 1.0
    d = 1.0-qab*x/qap

    if (abs(d) < FPMIN): d = FPMIN
    d = 1.0/d
    h = d

    for m in range(1, MAXIT+1):
        m2 = 2.0*m
        aa = m*(b-m)*x/((qam+m2)*(a+m2))
        d=1.0+aa*d
        if (abs(d) < FPMIN): d=FPMIN
        c=1.0+aa/c
        if (abs(c) < FPMIN): c=FPMIN
        d = 1.0/d
        h *= d*c
        aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2))

        d=1.0+aa*d
        if (abs(d) < FPMIN): d=FPMIN
        c=1.0+aa/c
        if (abs(c) < FPMIN): c=FPMIN
        d = 1.0/d

        delta = d*c
        h *= delta
        if (abs(delta-1.0) < EPS): break

    if (m > MAXIT):  print(("a or b too big or MAXIT too small "
                                           "in betacf"), file=sys.stderr)
    return h


def gammaln(x):
    """
    Compute log gamma function
    """
    xx = x
    s = 1.000000000190015
    for i in range(0, 6):
        xx += 1
        s += gamma_c[i]/xx

    res = ((x+0.5) * log(x+5.5)) - (x+5.5) + log(2.5066282746310005*s/x)
    if (res >= 0):
        return res
    else:
        return 0					# avoid roundoff error


def die(string):
    print(string, file=sys.stderr)

