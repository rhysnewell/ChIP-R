#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 14:23:48 2017

@author: rhys newell
"""
#
from chipr import bed
from chipr import ival
import numpy as np
from functools import reduce
import operator
import math
import random
import copy
from chipr import multipletesting
import scipy.stats


def rankmetric(elem):
    return elem[1]


def orderentries(vec, reverse=True):
    return sorted(range(len(vec)), key=vec.__getitem__, reverse=reverse)


def prod(factors):
    return reduce(operator.mul, factors, 1)


def rankentries(vec, reverse=True, duphandling='average', random_seed=0.5):
    n = len(vec)
    vec1 = orderentries(vec, reverse=reverse)  # Ordered vector of entries
    vec2 = [vec[r] for r in vec1]  # List of ranks
    sumr = 0  # sum of ranks
    dups_amount = []
    duplicates = 0  # Count of duplicate ranks
    narray = [0] * n  # Empty array

    for i in range(n):
        sumr += i
        duplicates += 1

        if i == n - 1 or vec2[i] != vec2[i + 1]:
            dups_amount.append(duplicates)

            if duphandling == 'randomize':
                srank = i + 1 - duplicates
                rand_order = randomorder(srank, duplicates, random_seed=random_seed)
                it_rand = iter(rand_order)

                for j in range(i - duplicates + 1, i + 1):
                    narray[vec1[j]] = next(it_rand)
            elif duphandling =='average':
                averank = sumr / float(duplicates) + 1

                for j in range(i - duplicates + 1, i + 1):
                    narray[vec1[j]] = averank

            sumr = 0
            duplicates = 0

    return narray


def randomorder(i, n, random_seed=0.5):
    new_order = []

    for j in range(n):
        new_order.append(i + j)

    random.shuffle(new_order, lambda: random_seed)

    return new_order


def compare_with_ties(x, y):
    if x < y:
        return -1
    elif x > y:
        return 1
    else:
        return random.randint(0, 1) * 2 - 1


# Helper function that returns the overall distance of an interval
def intervaldist(entry):
    dist = entry.chromEnd - entry.chromStart
    return dist


def reduceEntries(bedf, metric, specifyMax):

    if metric.lower().startswith('signalvalue'):
        met = 'signalValue'
    elif metric.lower().startswith('pvalue'):
        met = 'pValue'
    elif metric.lower().startswith('qvalue'):
        met = 'qValue'
    else:
        met = 'signalValue'

    newbedf = []
    for i in bedf:
        ents = []
        sig = []
        for e in i:
            ents.append(e)
            sig.append(getattr(e, met))
        sortedents = [x for _, x in sorted(zip(sig, ents), key=lambda pair: pair[0], reverse=True)]
        if specifyMax is not None:
            b = bed.BedFile(sortedents[0:specifyMax], 'Peaks')
        else:
            b = bed.BedFile(sortedents, 'Peaks')
        newbedf.append(b)

    return newbedf


def rankBed(bedf, rankmethod='signalvalue'):
    '''
    Function for ranking entries for multiple replicates
    Entries are ranked independently between replicates

    @param bedf         a list of all replicates in bedfile format

    @param rankmethod   determines what metric is used to rank the entries
                        'all', 'signalvalue', 'pvalue', 'qvalue'
                        all: will only work when value sother than -1 are provided for each entry

    @param duphandling  determines how to handles duplicates
                        'randomize' randomizes the order of the ranks
                        'average' gives the average of ranks
    '''
    newbedf = []
    for rep in bedf:
        svals = []
        # rep = union(bf, minentries=1)[0]
        zcount = 0
        logged = False
        for peak in rep:
            try:
                if rankmethod.lower().startswith('signalvalue'):
                    svals.append(peak.signalValue)
                elif rankmethod.lower().startswith('pvalue'):
                    svals.append(peak.pValue)
                    if peak.pValue > 1:
                        logged = True
                elif rankmethod.lower().startswith('qvalue'):
                    if peak.qValue > 1:
                        logged = True
            except AttributeError:
                zcount += 1
                if logged:
                    svals.append(0)
                else:
                    svals.append(1)
        # zero_up = np.arange(len(rep) - zcount, len(rep), 1).tolist()
        # zmean = np.mean(zero_up)
        if max(svals) > 1:
            ranks = scipy.stats.rankdata(svals, method='average')
            max_rank = max(ranks)
            ranks = [max_rank - i for i in ranks]
        else:
            ranks = scipy.stats.rankdata(svals, method='average')
        for ent, rank in zip(rep, ranks):
            if ent != 0:
                ent.addOption(rank=rank)
            else:
                pass
        # newbed = bed.BedFile(sortedBed)
        newbedf.append(rep)
    return newbedf

def rankreps(bedf, minentries=None, rankmethod='signalvalue', duphandling='average', random_seed=0.5, specifyMax='maximum'):

    if specifyMax is None:
        specifyMax = 'maximum'
    elif isNumeric(specifyMax):
        specifyMax = int(specifyMax)
        print("Filtering BedFiles to have "+str(specifyMax)+" entries...")
        bedf = reduceEntries(bedf, rankmethod, specifyMax)
        bedf = list(iterflatten(bedf))
        print('Done!')
    elif specifyMax.lower().startswith('maximum'):
        bedf = list(iterflatten(bedf))
    elif specifyMax.lower().startswith('minimum'):

        lens = []
        for i in bedf:
            lens.append(len(i))

        specifyMax = min(lens)
        print("Filtering BedFiles to have " + str(specifyMax) + " entries...")
        bedf = reduceEntries(bedf, rankmethod, specifyMax)
        bedf = list(iterflatten(bedf))
        print('Done!')

    unions = union(bedf, minentries)
    rankedUni = rankBed(unions[1], rankmethod)
    rp = [1]*len(unions[0])

    for idx, rep in enumerate(rankedUni):
        zeros = rep.count(0)
        # Averages from n-zeros to n
        zero_up = np.arange(len(rep)-zeros, len(rep), 1).tolist()
        # All pseudo-peaks in all replicates receive rank n
        # zero_up = [len(rep)]*zeros
        zmean = np.mean(zero_up)
        if math.isnan(zmean):
            zmean = 1
        print(str(zeros) + ' Pseudo-peaks used for replicate number: ' + str(idx + 1))
        for i, ent in enumerate(rep):
            try:
                rp[i] *= ent.rank
            except AttributeError:
                rp[i] *= zmean

    return unions, rp, rankedUni


def rankEntropy(unions):
    """
    Entropy

    function for determining the entropy of the assigned ranks in an attempt to determine a 'surprisal' value

    @param unions - Directly takes the output from the rankreps function

    @return entropy_values - A list of entropy values for the for the replicates

    """
    entries = unions[1]
    entropy_values = []
    min_rank = 1
    max_rank = len(entries[0])
    for e in range(len(entries[0])):
        entry_ranks = []
        for rep in range(len(entries)):
            entry_ranks.append(entries[rep][e][1])
        if max(entry_ranks) - min(entry_ranks) != 0:
            z = [((i - min(entry_ranks))/(max(entry_ranks) - min(entry_ranks)))+1 for i in entry_ranks]
        else:
            z = [1 for i in entry_ranks]
        n = sum(z)
        # print(n)
        probs = [j/n for j in z]
        ent = -sum([k*math.log(k) for k in probs])
        entropy_values.append(ent)

    return entropy_values


def entropy(pvals):
    n = len(pvals)
    maxent = math.log2(n)
    vals = []
    sums = []
    # sumPs = 0
    for idx, i in enumerate(pvals):
        vi = i*math.log10(i)
        vals.append(vi)

        if -sum(vals) >= maxent:
            print(i)
            print(idx)
            break

    return vals, sums

"""Plan to remove this implementation. Causing peaks to spill over chromosomes boundaries"""
def extend(unions, mindist=100):
    print('Extending intervals')
    extended = []

    for ind, i in enumerate(unions[0]):
        dist = intervaldist(i[0])
        if dist < mindist:
            i[0].chromStart = i[0].chromStart - min(abs(i[1]-i[0].chromStart), 50)
            i[0].chromEnd = i[0].chromEnd + min(abs(i[2]-i[0].chromEnd), 50)
        extended.append(i[0])
    return extended, unions[1]


# Function to calculate the product of assigned ranks
# Finds product of "up" and "down" ranks
def product(vec):
    # products_down = []
    products_up = []
    k = len(vec)
    ranks = []
    for i in range(len(vec[0])):
        # ranks_down = []
        ranks_up = []
        for e in range(len(vec)):
            # ranks_down.append(vec[e][i][1])
            ranks_up.append(vec[e][i][1])
        # products_down.append(prod(ranks_down))
        products_up.append(prod(ranks_up))
        ranks.append(ranks_up)

    return products_up, ranks


# Function to calculate the monotonic transformation of the rank values
def monotonetransformation(vec):
    vals = []
    k = len(vec)
    logn = math.log(len(vec[0]) + 1, 10)

    for i in range(len(vec[0])):
        ranks = []
        for e in range(len(vec)):
            ranks.append(vec[e][i][1])
        logranks = [math.log(x / (len(vec[0]) + 1)) for x in ranks]
        vals.append(-(-k * logn + sum(logranks)))

    return vals


def rankprod(reps, method='all', holdout=None):
    '''
    rankprod

    Function for calculating RP values based on ranks applied to bed entries

    reps      array of bed entries and their union from the rankReps function
              If reps come from rankReps then input reps as reps[1]
              to avoid inputting the union list

    method    Specify whether to use all replicates for the product or holdout some
              e.g. 'all', 'holdout'

    holdout   specify how many replicates to withdraw from the analysis,
              removes replicates from the end of the list
    '''
    if method == 'all':
        obsproducts = product(reps[1])
        # RP = [x**(1/len(reps[1])) for x in obsproducts]
        if len(obsproducts[0]) == len(reps[0]):
            return reps, obsproducts[0], obsproducts[1]
        else:
            return [], obsproducts[0], obsproducts[1]

    # Perform rankprod on only some of the reps
    # The withheld reps can be compared against results to observe consistency
    # Automatically holds out items at end of the list of reps
    elif method == 'holdout':
        if holdout is None:
            holdout = 1
        k = len(reps[1])
        obsproducts = product(reps[1][0:(k - holdout)])

        if len(obsproducts[0]) == len(reps[0]):
            return reps, obsproducts[0], obsproducts[1]
        else:
            print("Holdout has resulted in uneven entry lists, please review your input data")


def isNumeric(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def switch(x):
    return {'a': 1, 'b': 2}.get(x, 9)


def iterflatten(root):  # Function flatten list of lists into single list, not limited by depth
    if isinstance(root, (list, tuple)):
        for element in root:
            for e in iterflatten(element):
                yield e
    else:
        yield root


def rankprodbounds(rho, n, k, Delta):
    '''
    rankprodbounds

    Description

    This function computes bounds on the p-value for rank products.


    Arguments

    rho     a vector of integers corresponding to the rank products for which one wishes to
            compute the p-value. These are NOT the geometric means of the rank product values.

    n       the number of peaks.

    k       the number of replicates.

    Delta   a character string indicating whether an upper bound ('upper'), lower bound
            ('lower'), or geometric approximation ('geometric') should be computed.


    Returns

    A vector of p-values, one for each rank product.

    Details

    The exact p-value is guaranteed to be in between the lower and the upper bound. The
    geometric mean of the two bounds can be used as an approximation. Each bound is a piecewise
    continuous function of the rank product. The different pieces each have an analytic form,
    the parameters of which can be computed recursively.

    Note

    This implementation closely follows the description in Heskes, Eisinga, Breitling:
    "A fast algorithm for determining bounds and accurate approximate p-values of the
    rank product statistic for replicate experiments", BMC Bioinformatics, referred to
    below to as HEB.
    '''
    # Input Handling
    if any(rho) > ((n ** k) + 1) or any(rho) < 1:
        return 'rho out of bounds'
    if isNumeric(Delta) is False:
        if Delta == 'geometric':
            temp1 = rankprodbounds(rho, n, k, 'upper')
            temp2 = rankprodbounds(rho, n, k, 'lower')
            # pexact = [np.sqrt(x * y) for x, y in zip(temp1[1], temp2[1])]
            pvalue = [np.sqrt(x * y) for x, y in zip(temp1, temp2)]  # geometric mean of upper and lower bound
            return pvalue
        elif Delta == 'upper':
            Delta = 1  # computing upper and lower bound
        elif Delta == 'lower':
            Delta = 0

    # Compute Intervals that contain the rank products
    logn = math.log(n)
    logrho = [math.log(x) for x in rho]
    logrhon = [-(y / logn) for y in logrho]
    allj = [math.ceil(z + k) for z in logrhon]  # Index specifying the interval that contains rho
    minj = min(allj)  # Lowest interval index
    maxj = max(allj)  # Highest interval index

    # Initialize Parameters
    param = np.empty(shape=[k + 1, maxj + 1], dtype=object)
    for i in range(k + 1):
        for j in range(int(maxj + 1)):
            param[i][j] = [[] for i in range(3)]
            param[i][j].append([0.0])
            param[i][j].append([0.0])
            # param is a matrix of lists; each element of param is a list with values for the parameters
    # 0 through 4, which correspond to the parameters alpha through epsilon in HEB;
    # specifically, param[i+1][j+1][0] corresponds to alpha_{i,j} in HEB, etc, where the offset
    # of 1 is introduced to be able to represent, for example, alpha_{0,0};
    # 0, 1, and 2 can be vectors (with possibly different lengths for different i and j),
    # 3 and 4 are scalars

    # COMPUTE PARAMETERS

    for j in range(minj, maxj + 1):
        param = updateparam(param, n, k, j, Delta)

    # call to the function updateparam which recursively computes all parameters that are needed
    # to calculate the p-value for a rank product rho that lies in the interval with index j

    # Compute Rank Products Given Parameters

    k1 = k + 1
    P = [0] * len(rho)
    G = [0] * len(rho)  # G is a vector of the same length as rho, for each rho bounding the number of rank products
    for j in range(int(minj), int(maxj + 1)):
        j1 = j + 1
        iii = [index for index, item in enumerate(allj) if item == j]  # Indices of all rps that fall in interval j
        thisrho = [rho[x] for x in iii]
        thisparam = copy.deepcopy(param[k1 - 1][j1 - 1])
        thisG = copy.deepcopy(thisparam[4])
        if j != 0:
            nrho = len(thisrho)
            nterms = len(thisparam[0])

            rhoparam = [i * thisparam[3][0] for i in thisrho]
            rhoparam = list(iterflatten(rhoparam))

            thisG = [thisG[0] + i for i in rhoparam]
            thisG = list(iterflatten(thisG))

            mat_thisparam2 = np.asarray(thisparam[2]).reshape((1, len(thisparam[2])))
            mat_thisrho = np.asarray(thisrho).reshape((1, len(thisrho)))

            # print(mat_thisparam2.shape)
            # print(mat_thisrho.shape)
            d1 = np.dot(mat_thisparam2.T, mat_thisrho)
            log_thisrho = []
            for x in thisrho:
                if x != 0:
                    log_thisrho.append(math.log(x))

            log_thisrho = np.asarray(log_thisrho).reshape((1, len(log_thisrho)))
            logn_thisparam1 = [logn * (k - j + x) for x in thisparam[1]]
            logn_thisparam1 = list(iterflatten(logn_thisparam1))

            mat_thisparam1 = np.asarray(logn_thisparam1).reshape((1, len(logn_thisparam1)))
            d2 = np.tile(log_thisrho, (nterms, 1)) - np.tile(mat_thisparam1, (nrho, 1)).T
            mat_thisparam0 = np.asarray(thisparam[0]).reshape((1, len(thisparam[0]))).T
            d3 = np.tile(mat_thisparam0, nrho)
            d2_expd3 = np.power(d2, d3)

            mats = np.multiply(d1, d2_expd3)
            mat_sum = np.sum(mats, axis=0)
            mat_sum = mat_sum.tolist()
            mat_sum = list(iterflatten(mat_sum))
            thisG = [x + y for x, y in zip(thisG, mat_sum)]

        it_thisG = iter(thisG)
        # it_matsum = iter(mat_sum)
        for i in iii:
            try:
                G[i] = next(it_thisG)
            except StopIteration:
                continue
            # P[i] = next(it_matsum)
        G = list(iterflatten(G))
        # P = list(iterflatten(P))

    # pexact = [i / (n ** k) for i in P]
    pvalue = [i / (n ** k) for i in G]

    return pvalue


def updateparam(param, n, k, j, Delta):
    '''

    updateparam

    Description

    This subroutine updates the current set of parameters to make sure that the parameters
    corresponding to k replicates and the j'th interval are included.

    Arguments

    param   a matrix of lists, where each element of param is a list with values for the
            parameters a through e; these parameters specify the functional form of the bound;
            a, b, and c are all vectors of unknown length, d and e are scalars.
    n       the number of peaks
    k       the number of replicates for which we need to compute the corresponding parameters.
    j       the index of the interval for which we need to compute the corresponding parameters.
    Delta   0 for the lower bound and 1 for the upper bound.

    Value

    A possibly updated set of parameters, at least including those corresponding to (k,j).

    Details

    This subroutine make sure that the parameters corresponding to k replicates and a rank product
    within the j'th interval are included. If they already are (because calculated before), it
    does not compute anything. Otherwise, it recursively computes all parameters
    that are needed to arrive at the parameters for (k,j).

    Note

    This implementation closely follows HEB, in particular equations (9) through (11).
    '''
    k1 = k + 1
    j1 = j + 1

    try:
        param[k1 - 1][j1 - 1]
    except IndexError:
        param = np.append(param, np.empty(shape=[1, j1], dtype=object), 0)  # empty, so needs calculating

        for i in range(j1 + 1):
            param[i][j1 - 1] = [[] for i in range(3)]
            param[i][j1 - 1].append([0.0])
            param[i][j1 - 1].append([0.0])

    if param[k1 - 1][j1 - 1][4][0] == 0.0:  # empty, so needs calculating
        if j == 0:
            param[k1 - 1][j1 - 1][4][0] = n ** k
            param[k1 - 1][j1 - 1][3][0] = 0.0
        else:
            k0 = k1 - 1
            j0 = j1 - 1
            param = updateparam(param, n, k - 1, j - 1, Delta)
            # checking that the parameters for (k-1,j-1) that are needed to compute the
            # parameters for (k,j) are indeed available; if not, they are themselves computed
            param00 = copy.deepcopy(param[k0 - 1][j0 - 1])
            try:
                newa0 = [x + 1 for x in param00[0]]
            except TypeError:
                newa0 = []
            newa0 = list(iterflatten(newa0))
            newb0 = copy.deepcopy(param00[1])
            newb0 = list(iterflatten(newb0))
            try:
                newc0 = [x / y for x, y in zip(param00[2], newa0)]
            except TypeError:
                newc0 = []
            newc0 = list(iterflatten(newc0))
            param11 = copy.deepcopy(param00)
            if k == j:  # Calculates param (a to e) when k = j
                param11[4] = [(1.0 - Delta) * (1.0 - x) for x in param00[4]]
                param11[4] = list(iterflatten(param11[4]))

                param11[3] = [Delta * x + y for x, y in zip(param00[3], param00[4])]
                param11[3] = list(iterflatten(param11[3]))

                param11[0] = [1.0, copy.deepcopy(param00[0]), copy.deepcopy(newa0)]
                param11[0] = list(iterflatten(param11[0]))

                param11[1] = [0.0, copy.deepcopy(param00[1]), copy.deepcopy(newb0)]
                param11[1] = list(iterflatten(param11[1]))

                Dparam002 = [Delta * x for x in param00[2]]
                param11[2] = [copy.deepcopy(param00[3]), copy.deepcopy(Dparam002), copy.deepcopy(newc0)]
                param11[2] = list(iterflatten(param11[2]))

            else:  # Calculates parameters (a to e)
                param = updateparam(param, n, k - 1, j, Delta)
                param01 = copy.deepcopy(param[k0 - 1][j1 - 1])

                logn = math.log(n)
                lognnkj = (k - j) * math.log(n)

                newa1 = [x + 1.0 for x in param01[0]]
                newa1 = list(iterflatten(newa1))

                newa = [copy.deepcopy(newa0), copy.deepcopy(newa1)]
                newa = list(iterflatten(newa))

                newb = [copy.deepcopy(newb0), copy.deepcopy(param01[1])]
                newb = list(iterflatten(newb))

                paramnewa = [-x / y for x, y in zip(param01[2], newa1)]
                newc = [copy.deepcopy(newc0), copy.deepcopy(paramnewa)]
                newc = list(iterflatten(newc))

                diffparam = [x - y for x, y in zip(param00[4], param01[4])]
                diffparam = list(iterflatten(diffparam))

                Dparam = [(Delta - 1) * x for x in diffparam]
                Dparam = list(iterflatten(Dparam))

                nparam = [n * param01[4][0]]
                param11[4] = [x + y for x, y in zip(nparam, Dparam)]
                param11[4] = list(iterflatten(param11[4]))

                lognminb = [[-1 * x * logn for x in param00[1]], [(1 - x) * logn for x in param01[1]]]
                lognminb = list(iterflatten(lognminb))

                logexp = [x ** y for x, y in zip(lognminb, newa)]
                logexp = list(iterflatten(logexp))

                newcxlog = [x * y for x, y in zip(newc, logexp)]
                newcxlog = list(iterflatten(newcxlog))

                paraDiv = (param00[4][0] - param01[4][0]) / np.exp(lognnkj)
                Dnparam013 = (1 - Delta) * param01[3][0] / n

                param11[3] = [Delta * x + Dnparam013 + paraDiv - sum(newcxlog) for x in param00[3]]
                param11[3] = list(iterflatten(param11[3]))

                param11[0] = [1.0, 1.0, copy.deepcopy(param00[0]), copy.deepcopy(param01[0]), copy.deepcopy(newa)]
                param11[0] = list(iterflatten(param11[0]))

                param11[1] = [0.0, 1.0, copy.deepcopy(param00[1]), copy.deepcopy(param01[1]), copy.deepcopy(newb)]
                param11[1] = list(iterflatten(param11[1]))

                Dparamc = [Delta * x for x in param00[2]]
                Dnewparamc = [(1 - Delta) * x / n for x in param01[2]]

                negparam013 = [-x for x in param01[3]]
                param11[2] = [copy.deepcopy(param00[3]), copy.deepcopy(negparam013), copy.deepcopy(Dparamc),
                              copy.deepcopy(Dnewparamc), copy.deepcopy(newc)]
                param11[2] = list(iterflatten(param11[2]))
            param[k1 - 1][j1 - 1] = makeunique(param11)

    return param


def makeunique(param):
    '''
    makeunique

    Description

    This subroutine updates the parameters for a specific number of replicates and interval
    such that it contains only unique combinations of the parameters a and b.

    Arguments

    param   a single list with values for the parameters a through e; these parameters
            specify the functional form of the bound; a, b, and c are all vectors of
            unknown length, d and e are scalars.

    Value

    A possibly updated and then more concise set of parameters containing only unique
    combinations of the parameters a and b.

    Details

    While updating the vectors a and b, one may end up with the exact same combinations of
    a and b. Given the functional form of the bound, the representation can then be made more
    concise by simply adding the corresponding elements of c.
    '''

    p1 = [x for x in param[0]]
    p2 = [x for x in param[1]]
    ab = list(zip(p1, p2))

    abunique = []  # The list of unique sets from pairing p1 and p2
    for i in ab:
        if i not in abunique:
            abunique.append(i)

    parr = np.array(abunique)
    ab = np.array(ab)
    nunique = len(abunique)
    param[0] = parr[:, 0].T.tolist()
    param[1] = parr[:, 1].T.tolist()
    newc = [0] * nunique

    for i in range(nunique):
        iii1 = np.where(ab[:, 0] == parr[i, 0])[0].tolist()
        iii2 = np.where(ab[:, 1] == parr[i, 1])[0].tolist()
        iii3 = list(set(iii1).intersection(iii2))
        l = []
        for j in iii3:
            l.append(param[2][j])
        newc[i] = sum(iterflatten(l))
    param[2] = newc

    return param


def matmult(a, b):

    zip_b = zip(*b)
    zip_b = list(zip_b)

    return [[sum(ele_a * ele_b for ele_a, ele_b in zip(row_a, col_b))
             for col_b in zip_b] for row_a in a]


def clip(interval, N, min_N, minsize=0):
    '''''
    clip

    Determining the set of intervals that define the N-fold intersection between the input intervals
    Note that zero-length intervals are NOT adding to the intersection count

    Interval: list of intervals

    '''''

    if min_N < 1:
        return "clip: N must be 1 or greater"
    ret = []
    mins = []
    maxs = []
    for iv in interval:
        mins.append(iv.min)
        maxs.append(iv.max)
    mins = sorted(mins)
    maxs = sorted(maxs)
    avmin = sum(mins)//len(mins)
    avmax = sum(maxs)//len(maxs)
    idx_min = 0
    idx_max = 0
    cnt = 0
    current = None

    while idx_min < len(mins) and idx_max < len(maxs):
        if mins[idx_min] < maxs[idx_max]:
            cnt += 1
            if current is None and cnt >= min_N:
                try:
                    if mins[idx_min+1] < maxs[idx_max]:
                        current = ival.Interval(mins[idx_min], mins[idx_min + 1])
                        if (current.max - current.min) >= minsize:
                            ret.append(current)
                    else:
                        current = ival.Interval(mins[idx_min], maxs[idx_max])
                        if (current.max - current.min) >= minsize:
                            ret.append(current)
                        current = None
                except IndexError:
                    current = ival.Interval(mins[idx_min], maxs[idx_max])
                    if (current.max - current.min) >= minsize:
                        ret.append(current)
                    current = None

            elif current is not None and cnt >= min_N:
                if current.max < mins[idx_min]:
                    current = ival.Interval(current.max, mins[idx_min])
                    if (current.max - current.min) >= minsize:
                        ret.append(current)
            idx_min += 1

        elif mins[idx_min] > maxs[idx_max]:

            if current is not None and cnt == min_N:
                if current.max < mins[idx_min]:
                    current = ival.Interval(current.max, mins[idx_min])
                    if (current.max - current.min) >= minsize:
                        ret.append(current)
                elif current.max < maxs[idx_max]:
                    current = ival.Interval(current.max, maxs[idx_max])
                    if (current.max - current.min) >= minsize:
                        ret.append(current)
            elif current is not None and cnt >= min_N:
                if current.max < mins[idx_min]:
                    current = ival.Interval(current.max, mins[idx_min])
                    if (current.max - current.min) >= minsize:
                        ret.append(current)
                elif current.max < maxs[idx_max]:
                    current = ival.Interval(current.max, maxs[idx_max])
                    if (current.max - current.min) >= minsize:
                        ret.append(current)
            cnt -= 1
            idx_max += 1

        else:
            if current is None and cnt >= min_N:
                try:
                    if mins[idx_min+1] < maxs[idx_max]:
                        current = ival.Interval(mins[idx_min], mins[idx_min + 1])
                        if (current.max - current.min) >= minsize:
                            ret.append(current)
                    else:
                        current = ival.Interval(mins[idx_min], maxs[idx_max])
                        if (current.max - current.min) >= minsize:
                            ret.append(current)
                        current = None
                except IndexError:
                    current = ival.Interval(mins[idx_min], maxs[idx_max])
                    if (current.max - current.min) >= minsize:
                        ret.append(current)
                    current = None


            elif current is not None and cnt >= min_N:
                if current.max < mins[idx_min]:
                    current = ival.Interval(current.max, mins[idx_min])
                    if (current.max - current.min) >= minsize:
                        ret.append(current)
            if current is not None and cnt == min_N:
                if current.max < mins[idx_min]:
                    current = ival.Interval(current.max, mins[idx_min])
                    if (current.max - current.min) >= minsize:
                        ret.append(current)
                elif current.max < maxs[idx_max]:
                    current = ival.Interval(current.max, maxs[idx_max])
                    if (current.max - current.min) >= minsize:
                        ret.append(current)
            elif current is not None and cnt >= min_N:
                if current.max < mins[idx_min]:
                    current = ival.Interval(current.max, mins[idx_min])
                    if (current.max - current.min) >= minsize:
                        ret.append(current)
                elif current.max < maxs[idx_max]:
                    current = ival.Interval(current.max, maxs[idx_max])
                    if (current.max - current.min) >= minsize:
                        ret.append(current)
            idx_max += 1
            idx_min += 1

    if idx_min < len(mins):
        return "Invalid intervals"
    if idx_max < len(maxs) and current is not None:
        while idx_max < len(maxs):
            if current.max < maxs[idx_max] and cnt >= min_N:
                current = ival.Interval(current.max, maxs[idx_max])
                if (current.max - current.min) >= minsize:
                    ret.append(current)
            cnt -= 1
            idx_max += 1
    newret = []
    for idx, r in enumerate(ret):
        if (r.max - r.min) >= minsize:
            newret.append(r)

    return newret, avmin, avmax


def merge(chrom, start, end, stash, min_N, N):
    '''
    merge

    Method for merging multiple entries into one,
    all falling within the given interval (which defines the extremes of the entries).
    With min_N set to 0 or 1, the method simply flattens the entries,
    resulting in exactly one entry with its score set to the maximum score (or number of entries, if score is not set).
    With a greater min_N, the method may return 0, 1 or more entries.

    chrom   the chromosome identifier

    start   the start point of the default interval to make up the entry

    end     the end point of the default interval to make up the entry

    stash   the entries to be merged

    min_N   minimum number of entries for defining an entry

    '''
    entries = []
    if min_N > len(stash):  # No point in the interval can meet the minimum number of overlapping entries
        # clipped = clip(stash, 1)
        return entries
    strand_aware = True  # Only false if one or more of the entries is false
    cnt_pos_strand = 0
    cnt_neg_strand = 0
    strand_reverse = False
    maxscore = 0
    for e in stash:
        if e.usestrand is False:
            strand_aware = False
        else:
            if e.usestrand is True:
                cnt_pos_strand += 1
                cnt_neg_strand += 0
            else:
                cnt_pos_strand += 0
                cnt_neg_strand += 1
    if strand_aware is True:
        strand_reverse = (cnt_neg_strand > cnt_pos_strand)
    nOnStrand = 0
    for e in stash:
        if strand_aware is True:
            if e.strand != strand_reverse:
                try:
                    maxscore = max(maxscore, e.score)
                except AttributeError:
                    maxscore = max(maxscore, e.signalValue)
                nOnStrand += 1
            else:
                try:
                    maxscore = max(maxscore, e.score)
                except AttributeError:
                    maxscore = max(maxscore, e.signalValue)
                nOnStrand += 1
    if min_N < 1:
        if strand_aware is True:
            entry = bed.BedEntry(chrom, start, end)
            entry.addOption(strand=strand_reverse, signalValue=stash[0].signalValue)
        else:
            entry = bed.BedEntry(chrom, start, end)
            entry.addOption(signalValue=stash[0].signalValue)
        if maxscore == 0:
            entry.addOption(score=nOnStrand)
        else:
            entry.addOption(score=maxscore)
        # clipped = clip(entry, 1)
        entries = entry
        return entries, entry.chromStart, entry.chromEnd
    else:
        ilist = []
        signalSum = 0
        depth = 0
        for support in stash:
            ilist.append(ival.Interval(support.chromStart, support.chromEnd))
            signalSum += support.signalValue
            depth += 1
        clipped = clip(ilist, N, min_N)
        iset = clipped[0]

        for iv in iset:

            if strand_aware:
                entry = bed.BedEntry(chrom, iv.min, iv.max)
                entry.strand = strand_reverse
                entry.addOption(signalValue=signalSum, depth=depth)
            else:
                entry = bed.BedEntry(chrom, iv.min, iv.max)
                entry.addOption(signalValue=signalSum, depth=depth)
            if maxscore == 0:
                entry.addOption(score=nOnStrand)
                entry.addOption(signalValue=signalSum, depth=depth)
            else:
                entry.addOption(score=maxscore)
                entry.addOption(signalValue=signalSum, depth=depth)
            entries.append(entry)
        return entries, clipped[1], clipped[2]


def compareTo(this_min, this_max, that_min, that_max):
    if this_min < that_min:
        return -1
    elif this_min > that_min:
        return 1
    elif this_max < that_max:
        return -1
    elif this_max > that_max:
        return 1
    else:
        return 0


def start(vec):

    """
    Return the minimum start point of a vector of entries
    """
    vec = [x for x in vec if x != 0 and x is not None]
    low = vec[0].chromStart
    startent = vec[0]
    for i in vec:
        if i.chromStart < low:
            low = i.chromStart
            startent = i
    return startent


def union(bedfiles, minentries=2, maxdist=2):
    '''
    union

    Create a set of BED entries that represent the union of BED entries in the given BED files.

    bedfiles     the BED files

    minentries   the minimum number of entries that need to fall into the interval to create a union

    maxdist      Specifies the maximum distance between entries before they are seen as no longer connected

    return the BED entries, to be used to construct a new BED file
    '''

    bedfiles = list(iterflatten(bedfiles))
    unions = []
    idx = -1
    maxzeros = len(bedfiles) - minentries
    if maxzeros < 0:
        maxzeros = 0
    entries = [[] for i in range(len(bedfiles))]
    edx = -1
    chroms = []

    for bedfile in bedfiles:
        for i in bedfile.chroms:
            if i not in chroms:
                chroms.append(i)

    for chrom in chroms:
        genarr = [0] * len(bedfiles)  # the generators, one for each chrom
        entryarr = [0] * len(bedfiles)  # the entry array, one for each chrom
        entryarr_prev = [0] * len(bedfiles)  # the entry array for the previous entries in entryarr
        stored = [0]*len(bedfiles)
        index = 0

        for bedfile in bedfiles:
            genarr[index] = bedfile.generate(chrom)
            try:
                entryarr[index] = next(genarr[index], 0)
            except StopIteration:
                print("no more chroms")
            index += 1

        for i in range(len(entryarr)):
            entryarr_prev[i] = entryarr[i]

        if 0 in entryarr:  # Account for mismatched chromosomes between bedfiles
            print("mismatch for " + str(chrom) + ", it has been removed.")
            if entryarr.count(0) > maxzeros:
                continue

        buildme = None
        indexs = []  # A history of all indexs of entryarr that appear in stashed
        stashed = []  # where we keep entries that will be merged, once entries are flattened into one
        while True:

            # index of current "head" BED entry
            # always look for the "first" entry in terms of their order
            first = -1

            for i in range(len(entryarr)):
                if entryarr[i] is not None and entryarr[i] != 0:
                    if first == -1:
                        first = i
                    else:
                        if compareTo(entryarr[first].chromStart, entryarr[first].chromEnd, entryarr[i].chromStart,
                                     entryarr[i].chromEnd) > 0:
                            first = i

            if first != -1:
                chosen = entryarr[first]
                if buildme is None:
                    buildme = ival.Interval(chosen.chromStart, chosen.chromEnd)
                    stashed.append(chosen)
                    indexs.append(first)
                else:  # Checks to see if buildme is complete
                    chosen_ival = ival.Interval(chosen.chromStart, chosen.chromEnd)
                    if ival.isect(buildme, chosen_ival):  # buildme is not complete, so continue extending
                        buildme = ival.Interval(buildme.min, max(buildme.max, chosen.chromEnd))
                        stashed.append(chosen)
                        indexs.append(first)

                    elif buildme.dist(chosen_ival, centre2centre=True) < maxdist:

                        # buildme not complete, because the next entry is not beyond the max distance
                        buildme = ival.Interval(buildme.min, max(buildme.max, chosen.chromEnd))
                        stashed.append(chosen)
                        indexs.append(first)

                    else:
                        # buildme is complete, stash the buildme interval (if sufficiently supported),
                        # then create a new buildme
                        s_ind = set(indexs)
                        if len(stashed) >= minentries and len(s_ind) >= minentries and (buildme.max - buildme.min) > 0:
                            merged = merge(chrom, buildme.min, buildme.max, stashed, minentries, N=len(bedfiles))
                            extentries = merged[0]
                            try:
                                for x in range(len(extentries)):
                                    if merged[0][x].chromEnd - merged[0][x].chromStart > 0:
                                        unions.append(merged[0][x])
                                        for i in range(len(entryarr_prev)):
                                            if entryarr_prev[i] in stashed:
                                                check = entryarr_prev[i].getInterval()
                                                mergcheck = merged[0][x].getInterval()
                                                if check.isectStrict(mergcheck):
                                                    stored[i] = copy.deepcopy(entryarr_prev[i])
                                                else:
                                                    stored[i] = 0

                                            elif entryarr[i] in stashed:
                                                check = entryarr[i].getInterval()
                                                mergcheck = merged[0][x].getInterval()
                                                if check.isectStrict(mergcheck):
                                                    stored[i] = copy.deepcopy(entryarr[i])
                                                else:
                                                    stored[i] = 0

                                            else:
                                                stored[i] = 0
                                        for i in range(len(stored)):
                                            entries[i].append(stored[i])

                            except TypeError:
                                if merged[0].chromEnd - merged[0].chromStart > 0:

                                    unions.append(merged[0])
                                    for i in range(len(entryarr_prev)):
                                        if entryarr_prev[i] in stashed:
                                            check = entryarr_prev[i].getInterval()
                                            mergcheck = merged[0].getInterval()
                                            if check.isectStrict(mergcheck):
                                                stored[i] = copy.deepcopy(entryarr_prev[i])
                                            else:
                                                stored[i] = 0
                                        elif entryarr[i] in stashed:
                                            check = entryarr[i].getInterval()
                                            mergcheck = merged[0].getInterval()
                                            if check.isectStrict(mergcheck):
                                                stored[i] = copy.deepcopy(entryarr[i])
                                            else:
                                                stored[i] = 0
                                        else:
                                            stored[i] = 0
                                    for i in range(len(stored)):
                                        entries[i].append(stored[i])


                            idx += 1
                        buildme = ival.Interval(chosen.chromStart, chosen.chromEnd)
                        stashed = []
                        stashed.append(chosen)
                        indexs = []
                        indexs.append(first)
                entryarr_prev[first] = entryarr[first]
                entryarr[first] = next(genarr[first], 0)

            if first == -1:
                break

        # there is probably a final entry, which also needs to be stashed away
        s_ind = set(indexs)
        if buildme is not None and len(stashed) >= minentries and len(s_ind) >= minentries and (buildme.max - buildme.min) > 0:
            merged = merge(chrom, buildme.min, buildme.max, stashed, minentries, N=len(bedfiles))
            extentries = merged[0]
            try:
                for x in range(len(extentries)):
                    if merged[0][x].chromEnd - merged[0][x].chromStart > 0:

                        unions.append(merged[0][x])
                        for i in range(len(entryarr_prev)):
                            if entryarr_prev[i] in stashed:
                                check = entryarr_prev[i].getInterval()
                                mergcheck = merged[0][x].getInterval()
                                if check.isectStrict(mergcheck):
                                    stored[i] = copy.deepcopy(entryarr_prev[i])
                                else:
                                    stored[i] = 0

                            elif entryarr[i] in stashed:
                                check = entryarr[i].getInterval()
                                mergcheck = merged[0][x].getInterval()
                                if check.isectStrict(mergcheck):
                                    stored[i] = copy.deepcopy(entryarr[i])
                                else:
                                    stored[i] = 0

                            else:
                                stored[i] = 0
                        for i in range(len(stored)):
                            entries[i].append(stored[i])

            except TypeError:
                if merged[0].chromEnd - merged[0].chromStart > 0:

                    unions.append(merged[0])
                    for i in range(len(entryarr_prev)):
                        if entryarr_prev[i] in stashed:
                            check = entryarr_prev[i].getInterval()
                            mergcheck = merged[0].getInterval()
                            if check.isectStrict(mergcheck):
                                stored[i] = copy.deepcopy(entryarr_prev[i])
                            else:
                                stored[i] = 0
                        elif entryarr[i] in stashed:
                            check = entryarr[i].getInterval()
                            mergcheck = merged[0].getInterval()
                            if check.isectStrict(mergcheck):
                                stored[i] = copy.deepcopy(entryarr[i])
                            else:
                                stored[i] = 0
                        else:
                            stored[i] = 0
                    for i in range(len(stored)):
                        entries[i].append(stored[i])

            idx += 1

    return unions, entries


def overlaps(ent1, ent2, buffer=0):
    """
    Overlaps:
    defines whether two bed entries are overlapping with each other
    @param buffer: sets a buffer zone around entry1 extending it's space by that many bps
    """
    if (ent1.chromStart-buffer) <= ent2.chromStart <= (ent1.chromEnd+buffer):
        return True
    elif (ent1.chromStart-buffer) <= ent2.chromEnd <= (ent1.chromEnd+buffer):
        return True
    else:
        return False


def collapse(bedf, idr=False):
    """
    collapse:
    Function that collapses the final union bedfile AFTER the rank product analysis has been performed. If there are multiple
    Bed Entries overlapping one another, they will be merged preventing an overabundance of bed entries created by the
    union function.
    """
    allOl = []
    pvals = []
    if idr is True:
        for chrom in bedf.chroms:
            chromarr = bedf.generate(chrom)
            ol1 = []
            ol2 = []
            ol3 = []
            ol4 = []
            ol5 = []

            for entry in chromarr:

                if len(ol3) == 0:
                    cnt = 1
                    ol1.append(entry.chromStart)
                    ol2.append(entry.chromEnd)
                    ol3.append(entry)
                    try:
                        ol4.append(entry.score)
                    except AttributeError:
                        ol4.append(0)
                    ol5.append(entry.signalValue)

                else:
                    for o in ol3:
                        if overlaps(o, entry):
                            ol1.append(entry.chromStart)
                            ol2.append(entry.chromEnd)
                            ol3.append(entry)
                            try:
                                ol4.append(entry.score)
                            except AttributeError:
                                ol4.append(0)
                            ol5.append(entry.signalValue)

                            break
                        else:
                            cnt += 1

                if cnt > len(ol3):
                    if len(ol3) > 1:
                        mer = bed.BedEntry(chrom, min(ol1), max(ol2))
                        newScore = max(ol4)  # Calculate combined p-value
                        newSignal = sum(ol5)
                        mer.addOption(name='TBD',
                                      score=newScore,
                                      strand='.',
                                      signalValue=newSignal,
                                      pValue=-1,
                                      qValue=-1)
                        allOl.append(mer)

                    else:
                        ol3[0].addOption(name='TBD',
                                         score=max(ol4),
                                         strand='.',
                                         signalValue=sum(ol5),
                                         pValue=-1,
                                         qValue=-1)
                        allOl.append(ol3[0])

                    ol1 = []
                    ol2 = []
                    ol3 = []
                    ol4 = []
                    ol5 = []

                    cnt = 1
                    ol1.append(entry.chromStart)
                    ol2.append(entry.chromEnd)
                    ol3.append(entry)
                    try:
                        ol4.append(entry.score)
                    except AttributeError:
                        ol4.append(0)
                    ol5.append(entry.signalValue)

        return allOl
    else:
        for chrom in bedf.chroms:
            chromarr = bedf.generate(chrom)
            ol1 = []
            ol2 = []
            ol3 = []
            ol4 = []
            ol5 = []
            ol6 = []

            for entry in chromarr:

                if len(ol3) == 0:
                    cnt = 1

                    ol1.append(entry.chromStart)
                    ol2.append(entry.chromEnd)
                    ol3.append(entry)
                    try:
                        ol4.append(entry.pValue)
                    except AttributeError:
                        ol4.append(-1)
                    try:
                        ol5.append(entry.qValue)
                    except AttributeError:
                        ol5.append(-1)
                    try:
                        ol6.append(entry.signalValue)
                    except AttributeError:
                        ol6.append(-1)
                else:
                    for o in ol3:
                        if overlaps(o, entry):
                            ol1.append(entry.chromStart)
                            ol2.append(entry.chromEnd)
                            ol3.append(entry)
                            try:
                                ol4.append(entry.pValue)
                            except AttributeError:
                                ol4.append(-1)
                            try:
                                ol5.append(entry.qValue)
                            except AttributeError:
                                ol5.append(-1)
                            try:
                                ol6.append(entry.signalValue)
                            except AttributeError:
                                ol6.append(-1)
                            break
                        else:
                            cnt += 1
                if cnt > len(ol3):
                    if len(ol3) > 1:
                        mer = bed.BedEntry(chrom, min(ol1), max(ol2))
                        newp = scipy.stats.combine_pvalues(ol4, 'stouffer')[1] #Calculate combined p-value
                        pvals.append(newp)
                        newq = scipy.stats.combine_pvalues(ol5, 'stouffer')[1]
                        # qvals.append(newq)

                        mer.addOption(name='TBD',
                                      score=min([abs(int(125 * math.log2(newp))), 1000]),
                                      strand='.',
                                      signalValue=sum(ol6),
                                      pValue=newp,
                                      qValue=newq)
                        allOl.append(mer)

                    else:
                        try:
                            ol3[0].addOption(name='TBD',
                                             score=min([abs(int(125 * math.log2(ol4[0]))), 1000]),
                                             strand='.',
                                             signalValue=sum(ol6),
                                             pValue=ol4[0],
                                             qValue=ol5[0])
                            allOl.append(ol3[0])
                            pvals.append(ol4[0])
                        except ValueError:
                            ol3[0].addOption(name='TBD',
                                             score=min([abs(int(125 * math.log2(1))), 1000]),
                                             strand='.',
                                             signalValue=sum(ol6),
                                             pValue=-1,
                                             qValue=-1)
                            allOl.append(ol3[0])
                            pvals.append(-1)
                    ol1 = []
                    ol2 = []
                    ol3 = []
                    ol4 = []
                    ol5 = []
                    ol6 = []
                    cnt = 1
                    ol1.append(entry.chromStart)
                    ol2.append(entry.chromEnd)
                    ol3.append(entry)
                    try:
                        ol4.append(entry.pValue)
                    except AttributeError:
                        ol4.append(-1)
                    try:
                        ol5.append(entry.qValue)
                    except AttributeError:
                        ol5.append(-1)
                    try:
                        ol6.append(entry.signalValue)
                    except AttributeError:
                        ol6.append(-1)

        return allOl, pvals


def thresholdCalc(p, k=6):

    """
    thresholdCalc
    Function for calculating a suitable cutoff threshold for a list of P-values.
    Finds intersection point between the binomial cumulative distribution function and the rank product probabilities

    params:
    @p - A list of P-values from the rank product distribution

    """
    n = len(p)
    ps = sorted(p)
    bs = []
    for i, v in enumerate(ps):
        b = scipy.stats.binom.cdf(i, n, v)
        bs.append(b)

    idx = np.argwhere(np.diff(np.sign(np.array(ps) - np.array(bs))) != 0).reshape(-1) + 0
    pits = [ps[i] for i in idx]

    return bs, ps, pits, idx


def connect_entries(bedf, reps, default_min_peak=20, broadpeaks=False):
    """
    connect_entries
    Modified version of union function that takes only one bedfile input and connects entries that are within 1 base pair
    of each other and have a significantly close enough rank product p-value.

    Prevents over fragmentation of peaks.
    """
    chroms = bedf.chroms.keys()

    unions = []
    for chrom in chroms:
        genarr = bedf.generate(chrom)
        entryarr = next(genarr, 0)
        buildme = None
        indexs = []  # A history of all indexs of entryarr that appear in stashed
        stashed = []  # where we keep entries that will be merged, once entries are flattened into one
        while True:

            # index of current "head" BED entry
            # always look for the "first" entry in terms of their order
            if entryarr != 0:
                chosen = entryarr
            if buildme is None:

                buildme = ival.Interval(chosen.chromStart, chosen.chromEnd)
                buildme_vars = vars(chosen)

                stashed.append(chosen)

            else:  # Checks to see if buildme is complete
                chosen_ival = ival.Interval(chosen.chromStart, chosen.chromEnd)
                chosen_vars = chosen.pValue

                if buildme.dist(chosen_ival, centre2centre=False) <= 1:
                    # Test if either interval is shorter than expected
                    lowerbound = calcLowerBound(reps, chosen, default_min_peak, broadpeaks)

                    if lowerbound >= chosen_ival.max - chosen_ival.min \
                            or lowerbound >= buildme.max - buildme.min:
                        # If so, test if they should be joined
                        if round(buildme_vars['pValue'], 1) == round(chosen_vars, 1):  # Test if peaks are similar
                            buildme = ival.Interval(buildme.min, max(buildme.max, chosen.chromEnd))
                            if buildme_vars['pValue'] > chosen_vars:
                                buildme_vars = vars(chosen)
                            stashed.append(chosen)
                        else:  # Treat peaks as independent entities

                            # buildme is complete, stash the buildme interval (if sufficiently supported),
                            # then create a new buildme
                            lowerbound = calcLowerBound(reps, chosen, default_min_peak, broadpeaks)

                            if lowerbound <= buildme.max - buildme.min:
                                merged = bed.BedEntry(chrom, buildme.min, buildme.max)
                                del buildme_vars['chrom'], buildme_vars['chromStart'], buildme_vars['chromEnd']
                                merged.addOption(**buildme_vars)
                                unions.append(merged)

                            buildme = ival.Interval(chosen.chromStart, chosen.chromEnd)
                            buildme_vars = vars(chosen)

                    elif round(buildme_vars['pValue'], 3) == round(chosen_vars, 3): # Test if peaks are similar

                        buildme = ival.Interval(buildme.min, max(buildme.max, chosen.chromEnd))
                        if buildme_vars['pValue'] > chosen_vars:
                            buildme_vars = vars(chosen)
                        stashed.append(chosen)
                        # entryarr = next(genarr, 0)
                    else: # Treat peaks as independent entities

                        # buildme is complete, stash the buildme interval (if sufficiently supported),
                        # then create a new buildme
                        lowerbound = calcLowerBound(reps, chosen, default_min_peak, broadpeaks)

                        if lowerbound <= buildme.max - buildme.min:
                            merged = bed.BedEntry(chrom, buildme.min, buildme.max)
                            del buildme_vars['chrom'], buildme_vars['chromStart'], buildme_vars['chromEnd']
                            merged.addOption(**buildme_vars)
                            unions.append(merged)

                        buildme = ival.Interval(chosen.chromStart, chosen.chromEnd)
                        buildme_vars = vars(chosen)

                else:

                    # buildme is complete, stash the buildme interval (if sufficiently supported),
                    # then create a new buildme
                    lowerbound = calcLowerBound(reps, chosen, default_min_peak, broadpeaks)

                    if lowerbound <= buildme.max - buildme.min:
                        merged = bed.BedEntry(chrom, buildme.min, buildme.max)
                        del buildme_vars['chrom'], buildme_vars['chromStart'], buildme_vars['chromEnd']
                        merged.addOption(**buildme_vars)
                        unions.append(merged)

                    buildme = ival.Interval(chosen.chromStart, chosen.chromEnd)
                    buildme_vars = vars(chosen)

                    stashed = []
                    stashed.append(chosen)

                # entryarr_prev = entryarr
                entryarr = next(genarr, 0)

                if entryarr == 0:
                    break

        # there is probably a final entry, which also needs to be stashed away
        s_ind = set(indexs)
        if buildme is not None:
            lowerbound = calcLowerBound(reps, chosen, default_min_peak, broadpeaks)

            if lowerbound <= buildme.max - buildme.min:
                merged = bed.BedEntry(chrom, buildme.min, buildme.max)
                del buildme_vars['chrom'], buildme_vars['chromStart'], buildme_vars['chromEnd']
                merged.addOption(**buildme_vars)
                unions.append(merged)

    unions = bed.BedFile(unions, 'IDR')

    return unions


def calcLowerBound(reps, chosen, default_min_peak=20, broadpeaks=False):
    closest_widths = []
    for rep in reps:
        closest_peaks = rep.getClosest(chosen)
        if closest_peaks is not None:
            for closest in closest_peaks:
                closest_widths.append(closest.chromEnd - closest.chromStart)
    if len(closest_widths) != 0 and not broadpeaks:
        mean_width = np.mean(closest_widths)
        lowerbound = mean_width - mean_width / 2
        upperbound = mean_width + mean_width / 2
    elif broadpeaks:
        lowerbound = default_min_peak
    else:
        lowerbound = default_min_peak
        upperbound = 450

    return lowerbound


def performrankprod(bedf, minentries=2, rankmethod="signalvalue", specifyMax=None,
                    duphandling='average', random_seed=0.5,
                    alpha=0.05,
                    filename="bedfile_unions.bed",
                    default_min_peak=20,
                    print_pvals=True,
                    fragment=False):
    print("Using high fragmentation:", fragment)

    # First create intersection and rank the entries in each replicate and return the rankproduct values
    ranks = rankreps(bedf, minentries, rankmethod, duphandling, random_seed, specifyMax)

    # Calculate rank product for each entry that contributes to a union entry
    # Calculate the pvalues of the rank product values
    print('Calculating rank product probabilities...')
    rpb_up = rankprodbounds(ranks[1], len(ranks[1]), len(bedf), 'geometric')
    print('Calculating binomial threshold...')
    # Calculate rpb and binomial intersection point
    Pks = thresholdCalc(rpb_up, k=len(bedf)-(minentries-1))
    if len(Pks[2]) != 0:
        binomAlpha = round(min(Pks[2]), 3)
    else:
        print('No binomial convergence, defaulting to 0.1')
        binomAlpha = 0.1

    # Perform multiple hypothesis testing correction upon the pvals
    fdr = multipletesting.fdrcorrection(rpb_up)

    # Determine whether to remove entries that are called significant
    print('Cleaning up output...')
    for i, v in enumerate(ranks[0][0]):
        p = rpb_up[i]

        if p != 0.0:
            ranks[0][0][i].addOption(name='TBD',
                                     score=min([abs(int(125 * math.log2(rpb_up[i]))), 1000]),
                                     strand='.',
                                     pValue=rpb_up[i],
                                     qValue=fdr[1][i])
        else:
            ranks[0][0][i].addOption(name='TBD',
                                     score=1000,
                                     strand='.',
                                     pValue=2.5e-20,
                                     qValue=2.5e-20)
    collapsed = bed.BedFile(ranks[0][0], 'IDR')
    len1 = len(collapsed)

    collapsed = connect_entries(collapsed, bedf, default_min_peak, fragment)

    t3cnt = 0
    t1_unions = []
    # t2_unions = []
    t3_unions = []
    pvals = []

    for i, v in enumerate(collapsed):
        """
        Determine what tier each union falls under
        Tier 1 - Union with significance score <= alpha ~~~~~~ OLD ~~~~~~~

        primary - Union with significance score <= binomAlpha.
                 These are unions that the binomial threshold calculation has deemed to be significant.

        secondary - Union that does not meet requirements for previous Tiers.
                 Should not be discarded as its peaks still appear in majority of replicates.
        """
        pvals.append(v.pValue)
        if v.pValue <= binomAlpha:
            v.addOption(name='primary_peak_' + str(i),
                        strand='.')
            t1_unions.append(v)
        # elif v.pValue <= binomAlpha:
        #     v.addOption(name='T2_peak_' + str(i),
        #                 strand='.')
        #     t2_unions.append(v)
        else:
            t3cnt += 1
            v.addOption(name='secondary_peak_' + str(i),
                        strand='.')
            t3_unions.append(v)

    sortedUnions = [x for _, x in sorted(zip(pvals, collapsed), key = lambda pair:pair[0])]
    print(round((len(t1_unions)/len(collapsed))*100, 2), "% Primary peaks")
    print(round((t3cnt/len(collapsed))*100, 2), "% Secondary peaks")

    if filename is not None:
        bed.writeBedFile(sortedUnions, filename + "_all.bed", format="Peaks")
        bed.writeBedFile(t1_unions, filename + "_optimal.bed", format="Peaks")
        mets = bed.BedFile(sortedUnions).getMetrics()
        with open(filename+'_log.txt', 'w') as f:
            f.write('Primary Peaks: '+str(len(t1_unions))+'\n'+
                    'Secondary Peaks: '+str(t3cnt)+'\n'+
                    'Total Peaks:'+str(t3cnt+len(t1_unions))+'\n')
            f.write('\n Chromosome'+'\t mean peak size\t standard deviation: \n')
            for k, v in mets.items():
                f.write(k+'\t'+str(v[0])+'\t'+str(v[1])+'\n')

    return collapsed, Pks, rpb_up, fdr