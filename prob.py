'''
Module for classes and functions that are representing and processing basic probabilities.
Uses and depends on "Alphabet" that is used to define discrete random variables.
'''
import random
from sym import *
from copy import deepcopy
import math

#################################################################################################
# Generic utility functions
#################################################################################################

def _getMeTuple(alphas, str):
    """ Handy function that resolves what entries that are being referred to in the case
    of written wildcards etc.
    Example y = _getMeTuple([DNA_Alphabet, Protein_Alphabet], '*R') gives y = (None, 'R')
    alphas: the alphabets
    str: the string that specifies entries (may include '*' and '-' signifying any symbol) """
    assert len(str) == len(alphas), "Entry invalid"
    if not type(str) is tuple:
        list = []
        for ndx in range(len(alphas)):
            if str[ndx] == '*' or str[ndx] == '-':
                list.append(None)
            else:
                list.append(str[ndx])
        return tuple(list)
    else:
        return str

#################################################################################################
# Distrib class
#################################################################################################

class Distrib():
    """ A class for a discrete probability distribution, defined over a specified "Alphabet"
        TODO: Fix pseudo counts
              Exclude from counts, specify in constructor,
              include only when computing probabilities by standard formula (n_a + pseudo_a * N^(1/2)) / (N + N^(1/2))
              Exclude from filesaves, include with filereads (optional)
    """
    def __init__(self, alpha, pseudo = 0.0):
        """ Construct a new distribution for a specified alphabet, using an optional pseudo-count.
        alpha: alphabet
        pseudo: either a single "count" that applies to all symbols, OR a distribution/dictionary with counts.
        """
        self.pseudo = pseudo or 0.0
        self.alpha = alpha
        self.cnt = [0.0 for _ in alpha]
        try: # assume pseudo is a dictionary or a Distrib itself
            self.tot = 0
            symndx = 0
            for sym in alpha:
                cnt = float(pseudo[sym])
                self.cnt[symndx] = cnt
                self.tot = self.tot + cnt
                symndx += 1
        except TypeError: # assume pseudo is a single count for each symbol
            self.cnt = [float(self.pseudo) for _ in alpha]
            self.tot = float(self.pseudo) * len(alpha) # track total counts (for efficiency)

    def observe(self, sym, cntme = 1.0):
        """ Make an observation of a symbol
        sym: symbol that is being observed
        cntme: number/weight of observation (default is 1)
        """
        ndx = self.alpha.symbols.index(sym)
        self.cnt[ndx] = self.cnt[ndx] + cntme
        self.tot = self.tot + cntme
        return

    def reset(self):
        """ Re-set the counts of this distribution. Pseudo-counts are re-applied. """
        try:
            self.tot = 0
            symndx = 0
            for sym in self.alpha: # assume it is a Distribution
                cnt = float(self.pseudo[sym])
                self.cnt[symndx] = cnt
                self.tot = self.tot + cnt
                symndx += 1
        except TypeError: # assume pseudo is a single count for each symbol
            self.cnt = [float(self.pseudo) for _ in self.alpha]
            self.tot = float(self.pseudo) * len(self.alpha) # track total counts (for efficiency)

    def reduce(self, new_alpha):
        """ Create new distribution from self, using (smaller) alphabet new_alpha. """
        d = Distrib(new_alpha, self.pseudo)
        for sym in new_alpha:
            d.observe(sym, self.cnt[self.alpha.index(sym)])
        return d

    def count(self, sym = None):
        """ Return the absolute count(s) of the distribution
            or the count for a specified symbol. """
        if sym != None:
            ndx = self.alpha.symbols.index(sym)
            return self.cnt[ndx]
        else:
            d = {}
            index = 0
            for a in self.alpha:
                d[a] = self.cnt[index]
                index += 1
            return d

    def add(self, distrib):
        """ Add the counts for the provided distribution to the present. """
        for i in range(len(self.cnt)):
            cnt = distrib.count(self.alpha[i])
            self.cnt[i] += cnt
            self.tot += cnt

    def subtract(self, distrib):
        """ Subtract the counts for the provided distribution from the present. """
        for i in range(len(self.cnt)):
            cnt = distrib.count(self.alpha[i])
            self.cnt[i] -= cnt
            self.tot -= cnt

    def getSymbols(self):
        return self.alpha.symbols

    def __getitem__(self, sym):
        """ Retrieve the probability of a symbol (ascertained by counts incl pseudo-counts) """
        if self.tot > 0.0:
            return self.count(sym) / self.tot
        else:
            return 1.0 / len(self.alpha) # uniform

    def prob(self, sym = None):
        """ Retrieve the probability of a symbol OR the probabilities of all symbols
        (listed in order of the alphabet index). """
        if sym != None:
            return self.__getitem__(sym)
        elif self.tot > 0:
            return [ s / self.tot for s in self.cnt ]
        else:
            return [ 1.0 / len(self.alpha) for _ in self.cnt ]

    def __iter__(self):
        return self.alpha

    def __str__(self):
        """ Return a readable representation of the distribution """
        str = '< '
        for s in self.alpha:
            str += (s + ("=%4.2f " % self[s]))
        return str + ' >'

    def swap(self, sym1, sym2):
        """ Swap the entries for specified symbols. Useful for reverse complement etc.
            Note that changes are made to the current instance. Use swapxcopy if you
            want to leave this instance intact. """
        sym1ndx = self.alpha.index(sym1)
        sym2ndx = self.alpha.index(sym2)
        tmpcnt = self.cnt[sym1ndx]
        self.cnt[sym1ndx] = self.cnt[sym2ndx]
        self.cnt[sym2ndx] = tmpcnt

    def swapxcopy(self, sym1, sym2):
        """ Create a new instance with swapped entries for specified symbols.
            Useful for reverse complement etc.
            Note that changes are NOT made to the current instance.
            Use swap if you want to modify this instance. """
        newdist = Distrib(self.alpha, self.count())
        newdist.swap(sym1, sym2)
        return newdist

    def writeDistrib(self, filename = None):
        """ Write the distribution to a file or string.
            Note that the total number of counts is also saved, e.g.
            * 1000 """
        str = ''
        for s in self.alpha:
            str += (s + ("\t%f\n" % self[s]))
        str += "*\t%d\n" % self.tot
        if filename != None:
            fh = open(filename, 'w')
            fh.write(str)
            fh.close()
        return str

    def generate(self):
        """ Generate and return a symbol from the distribution using assigned probabilities. """
        alpha = self.alpha
        p = random.random() # get a random value between 0 and 1
        q = 0.0
        for sym in alpha: # pick a symbol with a frequency proportional to its probability
            q = q + self[sym]
            if p < q:
                return sym
        return alpha[len(alpha)]

    def getmax(self):
        """ Generate the symbol with the largest probability. """
        maxprob = 0.0
        maxsym = None
        for sym in self.alpha:
            if self[sym] > maxprob or maxprob == 0.0:
                maxsym = sym
                maxprob = self[sym]
        return maxsym

    def getsort(self):
        """ Return the list of symbols, in order of their probability. """
        symlist = [sym for (sym, _) in self.getProbsort()]
        return symlist

    def getProbsort(self):
        """ Return the list of symbol-probability pairs, in order of their probability. """
        s = [(sym, self.prob(sym)) for sym in self.alpha]
        ss = sorted(s, key=lambda y: y[1], reverse=True)
        return ss

    def divergence(self, distrib2):
        """ Calculate the Kullback-Leibler divergence between two discrete distributions.
            Note that when self.prob(x) is 0, the divergence for x is 0.
            When distrib2.prob(x) is 0, it is replaced by 0.0001.
        """
        assert self.alpha == distrib2.alpha
        sum = 0.0
        base = len(self.alpha)
        for sym in self.alpha:
            if self[sym] > 0:
                if distrib2[sym] > 0:
                    sum += math.log(self[sym] / distrib2[sym]) * self[sym]
                else:
                    sum += math.log(self[sym] / 0.0001) * self[sym]
        return sum

    def entropy(self):
        """ Calculate the information (Shannon) entropy of the distribution.
            Note that the base is the size of the alphabet, so maximum entropy is by definition 1.
            Also note that if the probability is exactly zero, it is replaced by a small value to
            avoid numerical issues with the logarithm. """
        sum = 0.0
        base = len(self.alpha)
        for sym in self.alpha:
            p = self.__getitem__(sym)
            if p == 0:
                p = 0.000001
            sum +=  p * math.log(p, base)
        return -sum

def writeDistribs(distribs, filename):
    """ Write a list/set of distributions to a single file. """
    str = ''
    k = 0
    for d in distribs:
        str += "[%d]\n%s" % (k, d.writeDistrib())
        k += 1
    fh = open(filename, 'w')
    fh.write(str)
    fh.close()

def _readDistrib(linelist):
    """ Extract distribution from a pre-processed list if strings. """
    symstr = ''
    d = {}
    for line in linelist:
        line = line.strip()
        if len(line) == 0 or line.startswith('#'):
            continue
        sections = line.split()
        sym, value = sections[0:2]
        if len(sym) == 1:
            if sym != '*':
                symstr += sym
        else:
            raise RuntimeError("Invalid symbol in distribution: " + sym)
        try:
            d[sym] = float(value)
        except ValueError:
            raise RuntimeError("Invalid value in distribution for symbol " + sym + ": " + value)
    if len(d) == 0:
        return None
    alpha = Alphabet(symstr)
    if '*' in list(d.keys()): # tot provided
        for sym in d:
            if sym != '*':
                d[sym] = d[sym] * d['*']
    distrib = Distrib(alpha, d)
    return distrib

def readDistribs(filename):
    """ Load a list of distributions from file.
    Note that if a row contains '* <number>' then it is assumed that each probability
    associated with the specific distribution is based on <number> counts. """
    fh = open(filename)
    string = fh.read()
    distlist = []
    linelist = []
    for line in string.splitlines():
        line = line.strip()
        if line.startswith('['):
            if len(linelist) != 0:
                distlist.append(_readDistrib(linelist))
            linelist = []
        elif len(line) == 0 or line.startswith('#'):
            pass # comment or blank line --> ignore
        else:
            linelist.append(line)
    # end for-loop, reading the file
    if len(linelist) != 0:
        distlist.append(_readDistrib(linelist))
    fh.close()
    return distlist

def readDistrib(filename):
    """ Load a distribution from file.
    Note that if a row contains '* <number>' then it is assumed that each probability
    is based on <number> counts. """
    dlist = readDistribs(filename)
    if len(dlist) > 0:  # if at least one distribution was in the file...
        return dlist[0] # return the first

import re

def _readMultiCount(linelist, format = 'JASPAR'):
    ncol = 0
    symcount = {}
    if format == 'JASPAR2010':
        for line in linelist:
            line = line.strip()
            if len(line) > 0:
                name = line.split()[0]
                counts = []
                for txt in re.findall(r'\w+', line):
                    try:
                        y = float(txt)
                        counts.append(y)
                    except ValueError:
                        pass # ignore non-numeric entries
                if len(counts) != ncol and ncol != 0:
                    raise RuntimeError('Invalid row in file: ' + line)
                ncol = len(counts)
                if len(name) == 1: # proper symbol
                    symcount[name] = counts
        alpha = Alphabet(''.join(list(symcount.keys())))
        distribs = []
        for col in range(ncol):
            d = dict([(sym, symcount[sym][col]) for sym in symcount])
            distribs.append(Distrib(alpha, d))
    elif format == 'JASPAR':
        alpha_str = 'ACGT'
        alpha = Alphabet(alpha_str)
        cnt = 0
        for sym in alpha_str:
            line = linelist[cnt].strip()
            counts = []
            for txt in re.findall(r'\w+', line):
                try:
                    y = float(txt)
                    counts.append(y)
                except ValueError:
                    pass # ignore non-numeric entries
            if len(counts) != ncol and ncol != 0:
                raise RuntimeError('Invalid row in file: ' + line)
            ncol = len(counts)
            symcount[sym] = counts
            cnt += 1
        distribs = []
        for col in range(ncol):
            d = dict([(sym, symcount[sym][col]) for sym in symcount])
            distribs.append(Distrib(alpha, d))
    else:
        raise RuntimeError('Unsupported format: ' + format)
    return distribs

def readMultiCounts(filename, format = 'JASPAR'):
    """ Read a file of raw counts for multiple distributions over the same set of symbols
        for (possibly) multiple (named) entries.
        filename: name of file
        format: format of file, default is 'JASPAR' exemplified below
        >MA0001.1 SEP4
        0    3    79    40    66    48    65    11    65    0
        94    75    4    3    1    2    5    2    3    3
        1    0    3    4    1    0    5    3    28    88
        2    19    11    50    29    47    22    81    1    6
        returns a dictionary of Distrib's, key:ed by entry name (e.g. MA001.1)
    """
    fh = open(filename)
    linelist = []
    entryname = ''
    entries = {}
    for row in fh:
        row = row.strip()
        if len(row) < 1: continue
        if row.startswith('>'):
            if len(linelist) > 0:
                entries[entryname] = _readMultiCount(linelist, format=format)
                linelist = []
            entryname = row[1:].split()[0]
        else:
            linelist.append(row)
    if len(linelist) > 0:
        entries[entryname] = _readMultiCount(linelist, format=format)
    fh.close()
    return entries

def readMultiCount(filename, format = 'JASPAR'):
    """ Read a file of raw counts for multiple distributions over the same set of symbols.
        filename: name of file
        format: format of file, default is 'JASPAR' exemplified below
        0    3    79    40    66    48    65    11    65    0
        94    75    4    3    1    2    5    2    3    3
        1    0    3    4    1    0    5    3    28    88
        2    19    11    50    29    47    22    81    1    6
        returns a list of Distrib's
    """
    d = readMultiCounts(filename, format=format)
    if len(d) > 0:
        return list(d.values())[0]

#################################################################################################
# Joint class
#################################################################################################

class Joint(object):
    """ A joint probability class.
        The JP is represented as a distribution over n-tuples where n is the number of variables.
        Variables can be for any defined alphabet. The size of each alphabet determine the
        number of entries in the table (with probs that add up to 1.0) """

    def __init__(self, alphas):
        """ A distribution of n-tuples.
        alphas: Alphabet(s) over which the distribution is defined
        """
        if type(alphas) is Alphabet:
            self.alphas = tuple( [alphas] )
        elif type(alphas) is tuple:
            self.alphas = alphas
        else:
            self.alphas = tuple( alphas )
        self.store = TupleStore(self.alphas)
        self.totalCnt = 0

    def getN(self):
        """ Retrieve the number of distributions/random variables. """
        return len(self.alphas)

    def __iter__(self):
        return self.store.__iter__()

    def reset(self):
        """ Re-set the counts of this joint distribution. Pseudo-counts are re-applied. """
        for entry in self.store:
            self.store[entry] = None
        self.totalCnt = 0

    def observe(self, key, cnt = 1):
        """ Make an observation of a tuple/key
        key: tuple that is being observed
        cnt: number/weight of observation (default is 1)
        """
        key = _getMeTuple(self.alphas, key)
        if not None in key:
            score = self.store[key]
            if (score == None):
                score = 0
            self.totalCnt += cnt
            self.store[key] = score + cnt
        else: # there are wildcards in the key
            allkeys = [mykey for mykey in self.store.getAll(key)]
            mycnt = float(cnt)/float(len(allkeys))
            self.totalCnt += cnt
            for mykey in allkeys:
                score = self.store[mykey]
                if (score == None):
                    score = 0
                self.store[mykey] = score + mycnt
        return

    def count(self, key):
        """ Return the absolute count that is used for the joint probability table. """
        key = _getMeTuple(self.alphas, key)
        score = self.store[key]
        if (score == None):
            score = 0.0
            for match in self.store.getAll(key):
                y = self.store[match]
                if y != None:
                    score += y
        return score

    def __getitem__(self, key):
        """ Determine and return the probability of a specified expression of the n-tuple
        which can involve "wildcards"
        Note that no assumptions are made regarding independence. """
        key = _getMeTuple(self.alphas, key)
        score = self.store[key]
        if (score == None):
            score = 0.0
            for match in self.store.getAll(key):
                y = self.store[match]
                if y != None:
                    score += y
        if self.totalCnt == 0:
            return 0.0
        return float(score) / float(self.totalCnt)

    def __str__(self):
        """ Return a textual representation of the JP. """
        str = '< '
        if self.totalCnt == 0.0:
            return str + 'None >'
        for s in self.store:
            if self[s] == None:
                y = 0.0
            else:
                y = self[s]
            str += (''.join(s) + ("=%4.2f " % y))
        return str + ' >'

    def items(self, sort = False):
        """ In a dictionary-like way return all entries as a list of 2-tuples (key, prob).
        If sort is True, entries are sorted in descending order of probability.
        Note that this function should NOT be used for big (>5 variables) tables."""
        if self.totalCnt == 0.0:
            return []
        ret = []
        for s in self.store:
            if self[s] != None:
                ret.append((s, self[s]))
        if sort:
            return sorted(ret, key=lambda v: v[1], reverse=True)
        return ret


class IndepJoint(Joint):

    def __init__(self, alphas, pseudo = 0.0):
        """ A distribution of n-tuples.
        All positions are assumed to be independent.
        alphas: Alphabet(s) over which the distribution is defined
        """
        self.pseudo = pseudo
        if type(alphas) is Alphabet:
            self.alphas = tuple( [alphas] )
        elif type(alphas) is tuple:
            self.alphas = alphas
        else:
            self.alphas = tuple( alphas )
        self.store = [Distrib(alpha, pseudo) for alpha in self.alphas]

    def getN(self):
        """ Retrieve the number of distributions/random variables. """
        return len(self.alphas)

    def __iter__(self):
        return TupleStore(self.alphas).__iter__()

    def reset(self):
        """ Re-set the counts of each distribution. Pseudo-counts are re-applied. """
        self.store = [Distrib(alpha, self.pseudo) for alpha in self.alphas]

    def observe(self, key, cnt = 1, countGaps = True):
        """ Make an observation of a tuple/key
        key: tuple that is being observed
        cnt: number/weight of observation (default is 1)
        """
        assert len(key) == len(self.store), "Number of symbols must agree with the number of positions"
        for i in range(len(self.store)):
            subkey = key[i]
            if subkey == '-' and countGaps == False:
                continue
            if subkey == '*' or subkey == '-':
                for sym in self.alphas[i]:
                    score = self.store[i][sym]
                    if (score == None):
                        score = 0
                    self.store[i].observe(sym, float(cnt)/float(len(self.alphas[i])))
            else:
                score = self.store[i][subkey]
                if (score == None):
                    score = 0
                self.store[i].observe(subkey, cnt)

    def __getitem__(self, key):
        """ Determine and return the probability of a specified expression of the n-tuple
        which can involve "wildcards"
        Note that variables are assumed to be independent. """
        assert len(key) == len(self.store), "Number of symbols must agree with the number of positions"
        prob = 1.0
        for i in range(len(self.store)):
            mykey = key[i]
            if mykey == '*' or mykey == '-':
                pass # same as multiplying with 1.0 (all symbols possible)
            else:
                prob *= self.store[i][mykey]
        return prob

    def get(self, sym, pos):
        """ Retrieve the probability of a specific symbol at a specified position. """
        mystore = self.store[pos]
        return mystore[sym]

    def getColumn(self, column, count = False):
        """ Retrieve all the probabilities (or counts) for a specified position.
            Returns values as a dictionary, with symbol as key."""
        d = {}
        for a in self.alphas[column]:
            if count: # absolute count
                d[a] = self.store[column].count(a)
            else: # probability
                d[a] = self.store[column][a]
        return d

    def getRow(self, sym, count = False):
        """ Retrieve the probabilities (or counts) for a specific symbol over all columns/positions.
            Returns a list of values in the order of the variables/alphabets supplied to the constructor. """
        d = []
        for store in self.store:
            if count: # absolute count
                d.append(store.count(sym))
            else: # probability
                d.append(store[sym])
        return d

    def getMatrix(self, count = False):
        """ Retrieve the full matrix of probabilities (or counts) """
        d = {}
        for a in self.alphas[0]:
            d[a] = self.getRow(a, count)
        return d

    def displayMatrix(self, count = False):
        """ Pretty-print matrix """
        print((" \t%s" % (''.join("\t%5d" % (i + 1) for i in range(len(self.alphas))))))
        for a in self.alphas[0]:
            if count:
                print(("%s\t%s" % (a, ''.join("\t%5d" % (y) for y in self.getRow(a, True)))))
            else:
                print(("%s\t%s" % (a, ''.join("\t%5.3f" % (y) for y in self.getRow(a)))))

    def __str__(self):
        """ Text representation of the table. Note that size is an issue so big tables
        will not be retrieved and displayed. """
        if self.alphas > 5:
            return '< ... too large to process ... >'
        tstore = TupleStore(self.alphas)
        str = '< '
        for key in tstore:
            p = 1.0
            for i in range(len(self.store)):
                value = self.store[i][key[i]]
                if value != None and value != 0.0:
                    p *= value
                else:
                    p = 0;
                    break;
            str += (''.join(key) + ("=%4.2f " % p))
        return str + ' >'

    def items(self, sort = False):
        """ In a dictionary-like way return all entries as a list of 2-tuples (key, prob).
        If sort is True, entries are sorted in descending order of probability.
        Note that this function should NOT be used for big (>5 variables) tables."""
        tstore = TupleStore(self.alphas)
        ret = []
        for key in tstore:
            p = 1.0
            for i in range(len(self.store)):
                value = self.store[i][key[i]]
                if value != None and value != 0.0:
                    p *= value
                else:
                    p = 0;
                    break;
            if p > 0.0:
                ret.append((key, p))
        if sort:
            return sorted(ret, key=lambda v: v[1], reverse=True)
        return ret

class NaiveBayes():
    """ NaiveBayes implements a classifier: a model defined over a class variable
        and conditional on a list of discrete feature variables.
        Note that feature variables are assumed to be independent. """

    def __init__(self, inputs, output, pseudo_input = 0.0, pseudo_output = 0.0):
        """ Initialise a classifier.
            inputs: list of alphabets that define the values that input variables can take.
            output: alphabet that defines the possible values the output variable takes
            pseudo_input: pseudo-count used for each input variable (default is 0.0)
            pseudo_output: pseudo-count used for the output variable (default is 0.0) """
        if type(inputs) is Alphabet:
            self.inputs = tuple( [inputs] )
        elif type(inputs) is tuple:
            self.inputs = inputs
        else:
            self.inputs = tuple( inputs )
        self.condprobs = {}   # store conditional probabilities as a dictionary (class is key)
        for outsym in output: # GIVEN the class
            # for each input variable initialise a conditional probability
            self.condprobs[outsym] = [ Distrib(input, pseudo_input) for input in self.inputs ]
        self.classprob = Distrib(output, pseudo_output) # the class prior

    def observe(self, inpseq, outsym):
        """ Record an observation of an input sequence of feature values that belongs to a class.
            inpseq: sequence/list of feature values, e.g. 'ATG'
            outsym: the class assigned to these feature values. """
        condprob = self.condprobs[outsym]
        for i in range(len(inpseq)):
            condprob[i].observe(inpseq[i])
        self.classprob.observe(outsym)

    def __getitem__(self, key):
        """ Determine and return the class probability GIVEN a specified n-tuple of feature values
        The class probability is given as an instance of Distrib. """
        out = Distrib(self.classprob.alpha)
        for outsym in self.classprob.getSymbols():
            condprob = self.condprobs[outsym]
            prob = self.classprob[outsym]
            for i in range(len(key)):
                prob *= condprob[i][key[i]] or 0.0
            out.observe(outsym, prob)
        return out
