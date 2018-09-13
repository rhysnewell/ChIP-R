"""
Module *** sequence ***

This module depends on the following modules

sym -- defines an alphabet
prob -- defines structures to hold probabilities (prob also depends on sym)

This module incorporates classes for

Sequence -- names and defines a sequence of symbols; computes various transformations and pairwise alignments
Alignment -- defines a multiple sequence alignment; computes stats for use in substitution matrices
SubstMatrix -- substitution matrix class to support alignment methods
Regexp -- defines patterns as regular expressions for textual pattern matching in sequences
PWM -- defines a weight matrix that can score any site in actual sequences

Incorporates methods for loading and saving files relevant to the above (e.g. FASTA, ALN, substitution matrices)
and methods for retrieving relevant data from web services

This code has been adapted to Python 3.5 in 2017

This code has gone through many updates and has benefited from kind contributions of course participants.
Please keep suggestions coming!
Email: m.boden@uq.edu.au
"""

import string, sys, re, math, os, array
import numpy
from webservice import *
from sym import *
from prob import *

# Sequence ------------------****

class Sequence(object):
    """ A biological sequence. Stores the sequence itself (as a compact array), 
    the alphabet (i.e., type of sequence it is), and optionally a name and further 
    information. """
    
    sequence = None # The array of symbols that make up the sequence 
    alphabet = None # The alphabet from which symbols come
    name =     None # The name (identifier) of a sequence
    info =     None # Other information (free text; e.g. annotations)
    length =   None # The number of symbols that the sequence is composed of
    gappy =    None # True if the sequence has "gaps", i.e. positions that represent deletions relative another sequence
    
    def __init__(self, sequence, alphabet = None, name = '', info = '', gappy = False):
        """ Create a sequence with the sequence data. Specifying the alphabet,
        name and other information about the sequence are all optional.
        The sequence data is immutable (stored as a string).
        Example:
        >>> myseq = Sequence('MVSAKKVPAIAMSFGVSF')
        will create a sequence with no name, and assign one of the predefined
        alphabets on the basis of what symbols were used.
        >>> myseq.alphabet.symbols
        will output the standard protein alphabet:
        ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q',
        'R', 'S', 'T', 'V', 'W', 'Y'] """
        
        self.sequence = sequence
        
        # Assign an alphabet
        # If no alphabet is provided, attempts to identify the alphabet from sequence
        self.alphabet = None
        if not alphabet is None:
            for sym in self.sequence:
                if not sym in alphabet and (sym != '-' or not gappy):  # error check: bail out
                    raise RuntimeError('Invalid symbol: %c in sequence %s' % (sym, name))
            self.alphabet = alphabet
        else:
            for alphaName in preferredOrder:
                alpha = predefAlphabets[alphaName]
                valid = True
                for sym in self.sequence:
                    if not sym in alpha and (sym != '-' or not gappy):  
                        valid = False
                        break
                if valid:
                    self.alphabet = alpha
                    break
            if self.alphabet is None:
                raise RuntimeError('Could not identify alphabet from sequence: %s' % name)
        
        # Store other information
        self.name = name
        self.info = info
        self.length = len(self.sequence)
        self.gappy = gappy
        
    def __len__(self):
        """ Defines what the "len" operator returns for an instance of Sequence, e.g.
        >>> seq = Sequence('ACGGTAGGA', DNA_Alphabet)
        >>> print (len(seq))
        9
        """
        return len(self.sequence)

    def __str__(self):
        """ Defines what should be printed when the print statement is used on a Sequence instance """
        str = self.name + ': '
        for sym in self:
            str += sym
        return str
    
    def __iter__(self):
        """ Defines how a Sequence should be "iterated", i.e. what its elements are, e.g.
        >>> seq = Sequence('AGGAT', DNA_Alphabet)
        >>> for sym in seq:
                print (sym)
        will print A, G, G, A, T (each on a separate row)
        """ 
        tsyms = tuple(self.sequence)
        return tsyms.__iter__()
    
    def __contains__(self, item):
        """ Defines what is returned when the "in" operator is used on a Sequence, e.g.
        >>> seq = Sequence('ACGGTAGGA', DNA_Alphabet)
        >>> print ('T' in seq)
        True
            which is equivalent to 
        >>> print (seq.__contains__('T'))
        True
        >>> print ('X' in seq)
        False
        """ 
        for sym in self.sequence:
            if sym == item:
                return True
        return False
        
    def __getitem__(self, ndx):
        """ Retrieve a specified index (or a "slice" of indices) of the sequence data.
            Calling self.__getitem__(3) is equivalent to self[3] 
        """
        if type(ndx) is slice:
            return ''.join(self.sequence[ndx])
        else:
            return self.sequence[ndx]
        
    def writeFasta(self):
        """ Write one sequence in FASTA format to a string and return it. """
        fasta = '>' + self.name + ' ' + self.info + '\n'
        data = ''.join(self.sequence)
        nlines = int(math.ceil((len(self.sequence) - 1) / 60 + 1))
        for i in range(nlines):
            lineofseq = ''.join(data[i*60 : (i+1)*60]) + '\n'
            fasta += lineofseq
        return fasta
    
    def count(self, findme = None):
        """ Get the number of occurrences of specified symbol findme OR
            if findme = None, return a dictionary of counts of all symbols in alphabet """
        if findme != None:
            cnt = 0
            for sym in self.sequence:
                if findme == sym:
                    cnt = cnt + 1
            return cnt
        else:
            symbolCounts = {}
            for symbol in self.alphabet:
                symbolCounts[symbol] = self.count(symbol)
            return symbolCounts

    def find(self, findme):
        """ Find the position of the specified symbol or sub-sequence """
        return ''.join(self.sequence).find(findme)

"""
Below are some useful methods for loading data from strings and files.
Recognize the FASTA format (nothing fancy).
"""
def readFasta(string, alphabet = None, ignore = False, gappy = False):
    """ Read the given string as FASTA formatted data and return the list of
        sequences contained within it.
        If alphabet is specified, use it, if None (default) then guess it.
        If ignore is False, errors cause the method to fail.
        If ignore is True, errors will disregard sequence.
        If gappy is False (default), sequence cannot contain gaps,
        if True gaps are accepted and included in the resulting sequences."""
    seqlist = []    # list of sequences contained in the string
    seqname = None  # name of *current* sequence
    seqinfo = None
    seqdata = []    # sequence data for *current* sequence
    for line in string.splitlines():    # read every line
        if len(line) == 0:              # ignore empty lines
            continue
        if line[0] == '>':  # start of new sequence
            if seqname:     # check if we've got one current
                try:
                    current = Sequence(seqdata, alphabet, seqname, seqinfo, gappy)
                    seqlist.append(current)
                except RuntimeError as errmsg:
                    if not ignore:
                        raise RuntimeError(errmsg)
            # now collect data about the new sequence
            seqinfo = line[1:].split() # skip first char
            if len(seqinfo) > 0:
                try:
                    parsed = parseDefline(seqinfo[0])
                    seqname = parsed[0]
                    seqinfo = line[1:]
                except IndexError as errmsg:
                    if not ignore:
                        raise RuntimeError(errmsg)
            else:
                seqname = ''
                seqinfo = ''
            seqdata = []
        else:               # we assume this is (more) data for current
            cleanline = line.split()
            for thisline in cleanline:
                seqdata.extend(tuple(thisline.strip('*')))
    # we're done reading the file, but the last sequence remains
    if seqname:
        try:
            lastseq = Sequence(seqdata, alphabet, seqname, seqinfo, gappy)
            seqlist.append(lastseq)
        except RuntimeError as errmsg:
            if not ignore:
                raise RuntimeError(errmsg)
    return seqlist

def parseDefline(string):
    """ Parse the FASTA defline (see http://en.wikipedia.org/wiki/FASTA_format)
        GenBank, EMBL, etc                gi|gi-number|gb|accession|locus
        SWISS-PROT, TrEMBL                sp|accession|name
        ...
        Return a tuple with
        [0] primary search key, e.g. UniProt accession, Genbank GI
        [1] secondary search key, e.g. UniProt name, Genbank accession
        [2] source, e.g. 'sp' (SwissProt/UniProt), 'tr' (TrEMBL), 'gb' (Genbank)
    """
    if len(string) == 0: return ('', '', '', '')
    s = string.split()[0]
    if re.match("^sp\|[A-Z][A-Z0-9]{5}\|\S+", s):            arg = s.split('|');  return (arg[1], arg[2], arg[0], '')
    elif re.match("^tr\|[A-Z][A-Z0-9]{5}\|\S+", s):          arg = s.split('|');  return (arg[1], arg[2], arg[0], '')
    elif re.match("^gi\|[0-9]*\|\S+\|\S+", s):               arg = s.split('|');  return (arg[1], arg[3], arg[0], arg[2])
    elif re.match("gb\|\S+\|\S+", s):                        arg = s.split('|');  return (arg[1], arg[2], arg[0], '')
    elif re.match("emb\|\S+\|\S+", s):                       arg = s.split('|');  return (arg[1], arg[2], arg[0], '')
    elif re.match("^refseq\|\S+\|\S+", s):                   arg = s.split('|');  return (arg[1], arg[2], arg[0], '')
    else: return (s, '', '', '')

def readFastaFile(filename, alphabet = None, ignore = False, gappy = False):
    """ Read the given FASTA formatted file and return the list of sequences
        contained within it. Note that if alphabet is NOT specified, it will take a
        separate guess for each sequence.
        If ignore is False, errors cause the method to fail.
        If ignore is True, errors will disregard sequence.
        If gappy is False (default), sequence cannot contain gaps,
        if True gaps are accepted and included in the resulting sequences."""
    fh = open(filename)
    seqlist = []
    batch = '' # a batch of rows including one or more complete FASTA entries
    rowcnt = 0
    for row in fh:
        row = row.strip()
        if len(row) > 0:
            if row.startswith('>') and rowcnt > 0:
                more = readFasta(batch, alphabet, ignore, gappy)
                if len(more) > 0:
                    seqlist.extend(more)
                batch = ''
                rowcnt = 0
            batch += row + '\n'
            rowcnt += 1
    if len(batch) > 0:
        more = readFasta(batch, alphabet, ignore, gappy)
        if len(more) > 0:
            seqlist.extend(more)
    fh.close()
    return seqlist

def writeFastaFile(filename, seqs):
    """ Write the specified sequences to a FASTA file. """
    fh = open(filename, 'w')
    for seq in seqs:
        fh.write(seq.writeFasta())
    fh.close()

def getMarkov(seqs, order = 0):
    """ Retrieve the Markov stats for a set of sequences. """
    myseqs = seqs
    if seqs is Sequence:
        myseqs = list([seqs])
    myalpha = None
    for seq in myseqs:
        if myalpha == None:
            myalpha = seq.alphabet
        else:
            if seq.alphabet != myalpha:
                raise RuntimeError('Sequence ' + seq.name + ' uses an invalid alphabet ')
    jp = Joint([myalpha for _ in range(order + 1)])
    for seq in myseqs:
        for i in range(len(seq) - order):
            sub = seq[i:i + order + 1]
            jp.observe(sub)
    return jp

def getCount(seqs, findme = None):
    if findme != None:
        cnt = 0
        for seq in seqs:
            cnt += seq.count(findme)
        return cnt
    else:
        if len(seqs) > 0:
            alpha = seqs[0].alphabet
            patcnt = {}
            for a in alpha:
                patcnt[a] = getCount(seqs, a)
        return patcnt

# Alignment ------------------

class Alignment():
    """ A sequence alignment class. Stores two or more sequences of equal length where
    one symbol is gap '-' 
    Example usage:
    >>> seqs = [Sequence('THIS-LI-NE-', Protein_Alphabet, gappy = True), Sequence('--ISALIGNED', Protein_Alphabet, gappy = True)]
    >>> print (Alignment(seqs))
     THIS-LI-NE-
     --ISALIGNED """

    alignlen = None
    seqs = None
    alphabet = None

    def __init__(self, seqs):
        self.alignlen = -1
        self.seqs = seqs
        self.alphabet = None
        for s in seqs:
            if self.alignlen == -1:
                self.alignlen = len(s)
            elif self.alignlen != len(s):
                raise RuntimeError("Alignment invalid: different lengths")
            if self.alphabet != None and self.alphabet != s.alphabet:
                raise RuntimeError("Alignment invalid: different alphabets")
            self.alphabet = s.alphabet

    def getnamelen(self):
        namelen = 0
        for seq in self.seqs:
            namelen = max(len(seq.name), namelen)
        return namelen

    def __len__(self):
        """ Defines what the "len" operator returns for an instance of Alignment, e.g.
        >>> seqs = [Sequence('THIS-LI-NE', Protein_Alphabet, gappy = True), Sequence('--ISALIGNED', Protein_Alphabet, gappy = True)]
        >>> aln = Alignment(seqs)
        >>> print(len(aln))
        2
        """
        return len(self.seqs)

    def getSize(self):
        """ Returns the size of an alignment in terms of number of columns """
        return self.alignlen

    def __str__(self):
        string = ''
        namelen = self.getnamelen()
        for seq in self.seqs:
            string += seq.name.ljust(namelen+1)
            for sym in seq:
                string += sym
            string += '\n'
        return string

    def __getitem__(self, ndx):
        return self.seqs[ndx]

    def writeClustal(self, filename = None):
        """ Write the alignment to a string or file using the Clustal file format. """
        symbolsPerLine = 60
        maxNameLength =  self.getnamelen() + 1
        string = ''
        wholeRows = self.alignlen / symbolsPerLine
        for i in range(int(wholeRows)):
            for j in range(len(self.seqs)):
                string += self.seqs[j].name.ljust(maxNameLength) + ' '
                string += self.seqs[j][i*symbolsPerLine:(i+1)*symbolsPerLine] + '\n'
            string += '\n'
        # Possible last row
        lastRowLength = self.alignlen - wholeRows*symbolsPerLine
        if lastRowLength > 0:
            for j in range(len(self.seqs)):
                if maxNameLength > 0:
                    string += self.seqs[j].name.ljust(maxNameLength) + ' '
                string += self.seqs[j][-lastRowLength:] + '\n'
        if filename != None:
            fh = open(filename, 'w')
            fh.write('CLUSTAL W (1.83) multiple sequence alignment\n\n\n') # fake header so that clustal believes it
            fh.write(string)
            fh.close()
            return
        return string

    def getProfile(self, pseudo = 0.0, countGaps = True):
        """ Determine the probability matrix from the alignment, assuming
        that each position is independent of all others. """
        p = IndepJoint([self.alphabet for _ in range(self.alignlen)], pseudo)
        for seq in self.seqs:
            p.observe(seq, 1, countGaps = countGaps)
        return p

    def getConsensus(self):
        """ Construct a consensus sequence. """
        syms = []
        for col in range(self.alignlen):
            d = Distrib(self.alphabet)
            for seq in self.seqs:
                if seq[col] in self.alphabet:
                    d.observe(seq[col])
            syms.append(d.getmax())
        return Sequence(syms)

    def getConsensusForColumn(self, colidx):
        symcnt = {}
        for seq in self.seqs:
            mysym = seq[colidx]
            try:
                symcnt[mysym] += 1
            except:
                symcnt[mysym] = 1
        consensus = None
        maxcnt = 0
        for mysym in symcnt:
            if symcnt[mysym] > maxcnt:
                maxcnt = symcnt[mysym]
                consensus = mysym
        return consensus

    def displayConsensus(self, theta1 = 0.2, theta2 = 0.05, lowercase = True):
        """ Display a table with rows for each alignment column, showing
            column index, entropy, number of gaps, and symbols in order of decreasing probability.
            theta1 is the threshold for displaying symbols in upper case,
            theta2 is the threshold for showing symbols at all, and in lower case. """
        print(("Alignment of %d sequences, with %d columns" % (len(self.seqs), self.alignlen)))
        print(("Column\tEntropy\tGaps\tProb\tConserv\tSymbols (Up>=%.2f;Low>=%.2f)\n" % (theta1, theta2)))
        for col in range(self.alignlen):
            d = Distrib(self.alphabet)
            gaps = 0
            for seq in self.seqs:
                if seq[col] in self.alphabet:
                    d.observe(seq[col])
                else:
                    gaps += 1
            print(((col + 1), "\t%5.3f" % d.entropy(), "\t%4d\t" % gaps,))
            symprobs = d.getProbsort()
            (_, maxprob) = symprobs[0]
            if maxprob >= theta1:
                print(("%d\tTRUE\t" % int(maxprob * 100),))
            else:
                print(("%d\t\t" % int(maxprob * 100),))
            for (sym, prob) in symprobs:
                if prob >= theta1:
                    print((sym, "%d%%" % int(prob * 100),))
                elif prob >= theta2 and lowercase:
                    print((sym.lower(), "%d%%" % int(prob * 100),))
                elif prob >= theta2:
                    print((sym, "%d%%" % int(prob * 100),))
            print()

    def saveConsensus(self, myseq, filename, theta1 = 0.2, theta2 = 0.05, lowercase = True, compact = False):
        """ Display a table with rows for each alignment column, showing
            column index, entropy, number of gaps, and symbols in order of decreasing probability.
            theta1 is the threshold for displaying symbols in upper case,
            theta2 is the threshold for showing symbols at all, and in lower case. """
        filename = ''.join(e for e in filename if e.isalnum() or e == '_' or e == '.')
        f = open(filename, 'w')
        f.write("Alignment of %d sequences, with %d columns\n" % (len(self.seqs), self.alignlen))
        if compact:
            f.write("Column\tConserv\tVariab\tAll (Up>=%.2f;Low>=%.2f)\n" % (theta1, theta2))
        else:
            f.write("Column\tProb\tConserv\tSymbols (Up>=%.2f;Low>=%.2f)\n" % (theta1, theta2))
        countrow = 0
        for col in range(self.alignlen):
            countrow += 1
            if myseq[col] == '-':
                continue
            alist = list(self.alphabet)
            alist.append('-')
            gapalphabet = Alphabet(alist)
            d_gap = Distrib(gapalphabet)
            d_nogap = Distrib(self.alphabet)
            for seq in self.seqs:
                if seq[col] in gapalphabet:
                    d_gap.observe(seq[col])
                if seq[col] in self.alphabet:
                    d_nogap.observe(seq[col])
            f.write("%d\t" % (col + 1))
            symprobs_nogap = d_nogap.getProbsort()
            symprobs_gap = d_gap.getProbsort()
            (maxsym, maxprob) = symprobs_nogap[0]
            if compact:
                if maxprob >= theta1:
                    f.write("%c\t" % maxsym)
                else:
                    f.write("\t")
                    for (sym, prob) in symprobs_gap:
                        if prob >= theta2 and lowercase:
                            f.write("%c" % sym.lower())
                        elif prob >= theta2:
                            f.write("%c" % sym)
                f.write("\t")
            else:
                if maxprob >= theta1:
                    f.write("%d\t" % int(maxprob * 100))
                else:
                    f.write("%d\t\t" % int(maxprob * 100))
            for (sym, prob) in symprobs_gap:
                if prob >= theta1:
                    f.write("%c %d%% " % (sym, int(prob * 100)))
                elif prob >= theta2 and lowercase:
                    f.write("%c %d%% " % (sym.lower(), int(prob * 100)))
                elif prob >= theta2:
                    f.write("%c %d%% " % (sym, int(prob * 100)))
            f.write('\n')
        f.close()

    def calcBackground(self):
        """ Count the proportion of each amino acid's occurrence in the
            alignment, and return as a probability distribution. """
        p = Distrib(self.alphabet)
        for seq in self.seqs:
            for sym in seq:
                if sym in self.alphabet: # ignore "gaps"
                    p.observe(sym)
        return p

    def calcSubstMatrix(self, background = None):
        """ Return a substitutionMatrix whose fg are based on this un-gapped
        multiple sequence alignment. Scores are given in half-bits. """
        # Get a list of the amino acids
        aminoAcids = self.alphabet.symbols
        columns = self.alignlen                   # Length of sequences in alignment
        numSeqs = len(self.seqs)                  # Number of sequences in alignment
        seqPairs = (numSeqs* (numSeqs - 1) ) / 2  # Number of pairs of sequences in ungapped alignment
        aaPairs = seqPairs * columns              # Number of pairs of amino acids in ungapped alignment
        # For each pair of amino acids, calculate the proportion of all aligned
        # amino acids in this alignment which are made up of that pair
        # (i.e., q[ab] = fab / aaPairs, where fab is the number of times
        #  a and b are aligned in this alignment)
        # See page 122 in Understanding Bioinformatics.
        q = {}
        for i in range( len(aminoAcids) ):
            a = aminoAcids[i]
            for j in range(i, len(aminoAcids)):
                b = aminoAcids[j]
                # Count the number of times each pair of amino acids is aligned
                fab = 0
                for column in range(columns):
                    # Count number of each amino acid in each column
                    col = [seq[column] for seq in self.seqs]
                    if a == b:
                        # Number of ways of pairing up n occurrences of amino
                        # acid a is n*(n-1)/2
                        cnt = col.count(a)
                        fab += cnt * (cnt-1)/2
                    else:
                        # Number of ways of pairing up n & m occurrences of
                        # amino acids a & b is n*m
                        fab += col.count(a)*col.count(b)
                # Calculate proportion of all aligned pairs of amino acids
                q[a+b] = q[b+a] = float(fab) / aaPairs
                if q[a+b] == 0:   # This is so we don't end up doing log(0)
                    q[a+b] = q[b+a] = 0.001
        # Background frequency calculation if required
        p = background or self.calcBackground()
        # Calculate log-odds ratio for each pair of amino acids
        s = SubstMatrix(self.alphabet)
        for a in aminoAcids:
            for b in aminoAcids:
                # Calculate random chance probabilitity (eab)
                if a == b:
                    eab = p[a]**2
                else:
                    eab = 2*p[a]*p[b]
                if eab == 0:
                    eab = 0.001
                # Calculate final score to be set in the substitution matrix
                odds = q[a+b] / eab
                sab = math.log(odds, 2) # log_2 transform
                sab = sab * 2 # units in half bits
                s.set(a, b, int(round(sab)))
        return s

    def calcDistances(self, measure, a=1.0):
        """ Calculate the evolutionary distance between all pairs of sequences
        in this alignment, using the given measure. Measure can be one of
        'fractional', 'poisson', 'gamma', 'jc' or 'k2p'. If 'gamma' or 'k2p' is
        given, then the parameter a must also be specified (or else it will use
        the default value of 1.0).
        Definitions of each distance metric are found in Zvelebil and Baum p268-276.
        These are mostly intended for DNA, but adapted for protein (as below).
        Note however that there are alternative distance matrices for proteins (p276).
        """
        measure = measure.lower()
        if not measure in ['fractional', 'poisson', 'gamma', 'jc', 'k2p']:
            raise RuntimeError('Unsupported evolutionary distance measure: %s' % measure)
        a = float(a)
        if len(self.alphabet) == 4:
            oneless = 3
            alphalen = 4
        elif len(self.alphabet) == 20:
            oneless = 19
            alphalen = 20
        else:
            raise RuntimeError('Invalid sequence alphabet: %s' % str(self.alphabet))
        distmat = numpy.zeros((len(self.seqs), len(self.seqs)))
        # Loop through each pair of sequences
        for i in range(len(self.seqs)):
            for j in range(i + 1, len(self.seqs)):
                seqA = self.seqs[i]
                seqB = self.seqs[j]
                # Calculate the fractional distance (p) first
                # The two sequences of interest are in seqA and seqB
                L = 0
                D = 0
                for k in range(self.alignlen):
                    # For every non-gapped column, put to L
                    # For every non-gapped column where the sequences are
                    # different, put to D
                    if seqA[k] != '-' and seqB[k] != '-':
                        L += 1
                        if seqA[k] != seqB[k]:
                            D += 1
                p = float(D)/L
                # Now calculate the specified measure based on p
                if measure == 'fractional':
                    dist = p
                elif measure == 'poisson':
                    dist = -math.log(1-p)
                elif measure == 'jc':
                    dist = -(float(oneless)/alphalen)*math.log(1 - (float(alphalen)/oneless)*p)
                elif measure == 'k2p':
                    dist = (float(oneless)/alphalen)*a*((1 - (float(alphalen)/oneless)*p)**(-1/a) - 1)
                else: # measure == 'gamma'
                    dist = a*((1-p)**(-1/a) - 1)
                distmat[i, j] = distmat[j, i] = dist
        return distmat

    def writeHTML(self, filename = None):
        """ Generate HTML that displays the alignment in color.
            Requires that the alphabet is annotated with the label 'html-color' (see Sequence.annotateSym)
            and that each symbol maps to a text string naming the color, e.g. 'blue'
        """
        html = '''<html><head><meta content="text/html; charset=ISO-8859-1" http-equiv="Content-Type">\n<title>Sequence Alignment</title>\n</head><body><pre>\n'''
        maxNameLength =  self.getnamelen()
        html += ''.ljust(maxNameLength) + ' '
        for i in range(self.alignlen - 1):
            if (i+1) % 10 == 0:
                html += str(i/10+1)[0]
            else:
                html += ' '
        html += '%s\n' % (self.alignlen)

        if self.alignlen > 10:
            html += ''.ljust(maxNameLength) + ' '
            for i in range(self.alignlen - 1):
                if (i+1) % 10 == 0:
                    index = len(str(i/10 + 1).split('.')[0])
                    html += str(i / 10 + 1).split('.')[0][(index * -1) + 1 ] if (len(str(i / 10 + 1).split('.')[0]) > 1) else '0'
                else:
                    html += ' '
            html += '\n'

        if self.alignlen > 100:
            html += ''.ljust(maxNameLength) + ' '
            for i in range(self.alignlen - 1):
                if (i+1) % 10 == 0 and i >= 99:
                    index = len(str(i/10 + 1).split('.')[0])
                    html += str(i / 10 + 1).split('.')[0][-1] if (len(str(i / 10 + 1).split('.')[0]) >2) else '0'

                else:
                    html += ' '
            html += '\n'

        if self.alignlen > 1000:
            html += ''.ljust(maxNameLength) + ' '
            for i in range(self.alignlen - 1):
                if (i+1) % 10 == 0:
                    html += '0' if (len(str(i / 10 + 1).split('.')[0]) > 2) else ' '

                else:
                    html += ' '
            html += '\n'
        for seq in self.seqs:
            html += seq.name.ljust(maxNameLength) + ' '
            for sym in seq:
                color = self.alphabet.getAnnotation('html-color', sym)
                if not color:
                    color = 'white'
                html += '<font style="BACKGROUND-COLOR: %s">%s</font>' % (color, sym)
            html += '\n'
        html += '</pre></body></html>'
        if filename:
            fh = open(filename, 'w')
            fh.write(html)
            fh.close()
        return html

def saveConsensus(aln, theta1 = 0.99, theta2 = 0.01, countgaps = False, consensus = True, filename = None):
    """ Display a table with rows for each alignment column, showing
        column index, entropy, number of gaps, and symbols in order of decreasing probability.
        theta1 is the percent threshold for consensus (when achieved, all other symbols are ignored)
        theta2 is the percent threshold for inclusion (symbols below are ignored).
        countgaps, if true, count gaps (default false).
        consensus, if true, always print the consensus symbol.
        filename is name of file to save the output to (default stdout)."""
    if filename == None:
        f = sys.stdout
    else:
        filename = ''.join(e for e in filename if e.isalnum() or e == '_' or e == '.')
        f = open(filename, 'w')
    if consensus:
        f.write("Alignment of %d sequences, with %d columns\n" % (len(aln.seqs), aln.alignlen))
        f.write("Consensus>=%.2f;Inclusion>=%.2f)\n" % (theta1, theta2))
    for col in range(aln.alignlen):
        # collect probabilities for column, with or without gap
        myalpha = aln.alphabet
        if countgaps:
            alist = list(aln.alphabet)
            alist.append('-')
            myalpha = Alphabet(alist)
        d = Distrib(myalpha)
        for seq in aln.seqs:
            if seq[col] in myalpha:
                d.observe(seq[col])
        symprobs = d.getProbsort() # the symbols sorted by probability
        ninclusions = 0
        for (s, p) in symprobs:
            if p >= theta2:
                ninclusions += 1
            else:
                break
        if consensus or ninclusions > 1:
            f.write("%d " % (col + 1))
        (maxs, maxp) = symprobs[0]
#        if maxp >= theta1 or consensus:
#            f.write("%c" % maxs)
        for (s, p) in symprobs[1:]:
            if p >= theta2:
                f.write("%c" % s)
        f.write("; ")
    f.write('\n')
    f.close()

def alignGlobal(seqA, seqB, substMatrix, gap = -1):
    """ Align seqA with seqB using the Needleman-Wunsch
    (global) algorithm. subsMatrix is the substitution matrix to use and
    gap is the linear gap penalty to use. """
    lenA, lenB = len(seqA), len(seqB)
    # Create the scoring matrix (S)
    S = numpy.zeros((lenA + 1, lenB + 1))
    # Fill the first row and column of S with multiples of the gap penalty
    for i in range(lenA + 1):
        S[i, 0] = i * gap
    for j in range(lenB + 1):
        S[0, j] = j * gap
    # Calculate the optimum score at each location in the matrix S
    # (where the score represents the best possible score for an alignment
    #  that ends at sequence indices i and j, for A and B, resp.)
    for i in range(1, lenA + 1):
        for j in range(1, lenB + 1):
            match  = S[i-1, j-1] + substMatrix.get(seqA[i-1], seqB[j-1])
            delete = S[i-1, j  ] + gap
            insert = S[i  , j-1] + gap
            S[i, j] = max([match, delete, insert])
    # Traceback the optimal alignment
    alignA = '' # a string for sequence A when aligned (e.g. 'THIS-LI-NE-', initially empty).
    alignB = '' # a string for sequence B when aligned (e.g. '--ISALIGNED', initially empty).
    # Start at the end (bottom-right corner of S)
    i = lenA
    j = lenB
    # Stop when we hit the beginning of at least one sequence
    while i > 0 and j > 0:
        if S[i, j] == S[i-1, j] + gap:
            # Got here by a gap in sequence B (go up)
            alignA = seqA[i-1] + alignA
            alignB = '-' + alignB
            i -= 1
        elif S[i, j] == S[i, j-1] + gap:
            # Got here by a gap in sequence A (go left)
            alignA = '-' + alignA
            alignB = seqB[j-1] + alignB
            j -= 1
        else:
            # Got here by aligning the bases (go diagonally)
            alignA = seqA[i-1] + alignA
            alignB = seqB[j-1] + alignB
            i -= 1
            j -= 1
    # Fill in the rest of the alignment if it begins with gaps
    # (i.e., traceback all the way to S[0, 0])
    while i > 0:
        # Go up
        alignA = seqA[i-1] + alignA
        alignB = '-' + alignB
        i -= 1
    while j > 0:
        # Go left
        alignA = '-' + alignA
        alignB = seqB[j-1] + alignB
        j -= 1
    return Alignment([Sequence(alignA, seqA.alphabet, seqA.name, gappy = True), Sequence(alignB, seqB.alphabet, seqB.name, gappy = True)])

def alignLocal(seqA, seqB, substMatrix, gap = -1):
    """ Align seqA with seqB using the Smith-Waterman
    (local) algorithm. subsMatrix is the substitution matrix to use and
    gap is the linear gap penalty to use. """
    lenA, lenB = len(seqA), len(seqB)
    # Create the scoring matrix (S)
    S = numpy.zeros((lenA + 1, lenB + 1))
    # Fill the first row and column of S with multiples of the gap penalty
    for i in range(lenA + 1):
        S[i, 0] = 0  # Local: init 0
    for j in range(lenB + 1):
        S[0, j] = 0  # Local: init 0
    # Calculate the optimum score at each location in the matrix S
    # (where the score represents the best possible score for an alignment
    #  that ends at sequence indices i and j, for A and B, resp.)
    for i in range(1, lenA + 1):
        for j in range(1, lenB + 1):
            match  = S[i-1, j-1] + substMatrix.get(seqA[i-1], seqB[j-1])
            delete = S[i-1, j  ] + gap
            insert = S[i  , j-1] + gap
            S[i, j] = max([match, delete, insert, 0])  # Local: add option that we re-start alignment from "0"
    # Trace back the optimal alignment
    alignA = ''
    alignB = ''
    # Local: start at the cell which has the highest score; find it
    i = 0
    j = 0
    for ii in range(1, lenA + 1):
        for jj in range(1, lenB + 1):
            if S[ii, jj] > S[i, j]:
                i = ii
                j = jj

    # Stop when we hit the end of a sequence
    # Local: also stop when we hit a score 0
    while i > 0 and j > 0 and S[i, j] > 0:
        if S[i, j] == S[i-1, j] + gap:
            # Got here by a gap in sequence B (go up)
            alignA = seqA[i-1] + alignA
            alignB = '-' + alignB
            i -= 1
        elif S[i, j] == S[i, j-1] + gap:
            # Got here by a gap in sequence A (go left)
            alignA = "-" + alignA
            alignB = seqB[j-1] + alignB
            j -= 1
        else:
            # Got here by aligning the bases (go diagonally)
            alignA = seqA[i-1] + alignA
            alignB = seqB[j-1] + alignB
            i -= 1
            j -= 1
    return Alignment([Sequence(alignA, seqA.alphabet, seqA.name, gappy = True), Sequence(alignB, seqB.alphabet, seqB.name, gappy = True)])

def tripletAlignGlobal(seqA, seqB, seqC, subsMatrix, gap = -1):
    """ Triplet-wise align this sequence with sequences seqB and seqC,
    using the Needleman-Wunsch (global) algorithm. subsMatrix is the
    substitution matrix to use and gap is the linear gap penalty to use. """

    lenA, lenB, lenC = [s.length for s in [seqA, seqB, seqC]]

    # Create the 3D scoring matrix
    traceback = numpy.zeros((lenA+1, lenB+1, lenC+1))
    # Fill the first row (in each dimension) with multiples of the gap penalty
    S = numpy.zeros((lenA+1, lenB+1, lenC+1))
    for i in range(lenA+1):
        S[i,0,0] = i * gap
    for j in range(lenB+1):
        S[0,j,0] = j * gap
    for k in range(lenC+1):
        S[0,0,k] = k * gap
    # Calculate the optimum __getitem__ at each location in the matrix
    for i in range(1, lenA+1):
        for j in range(1, lenB+1):
            for k in range(1, lenC+1):
                # Scored using sum-of-pairs
                matchABC = S[i-1, j-1, k-1] + subsMatrix.get(seqA[i-1], seqB[j-1]) \
                           + subsMatrix.get(seqA[i-1], seqC[k-1]) \
                           + subsMatrix.get(seqB[j-1], seqC[k-1])
                matchAB = S[i-1, j-1, k] + 2*gap + subsMatrix.get(seqA[i-1], seqB[j-1])
                matchBC = S[i, j-1, k-1] + 2*gap + subsMatrix.get(seqB[j-1], seqC[k-1])
                matchAC = S[i-1, j, k-1] + 2*gap + subsMatrix.get(seqA[i-1], seqC[k-1])
                gapAB = S[i, j, k-1] + 3*gap
                gapBC = S[i-1, j, k] + 3*gap
                gapAC = S[i, j-1, k] + 3*gap
                # Use maximum of the 7 options for this location
                S[i, j, k] = max([matchABC, matchAB, matchBC, matchAC, gapAB, gapBC, gapAC])
                # Remember which one was max., for the traceback
                if S[i, j, k] == matchABC:
                    traceback[i, j, k] = 0 #"matchABC"
                elif S[i, j, k] == matchBC:
                    traceback[i, j, k] = 1 #"matchBC"
                elif S[i, j, k] == matchAC:
                    traceback[i, j, k] = 2 #"matchAC"
                elif S[i, j, k] == matchAB:
                    traceback[i, j, k] = 3 #"matchAB"
                elif S[i, j, k] == gapAB:
                    traceback[i, j, k] = 4 #"gapAB"
                elif S[i, j, k] == gapBC:
                    traceback[i, j, k] = 5 #"gapBC"
                elif S[i, j, k] == gapAC:
                    traceback[i, j, k] = 6 #"gapAC"

    # Traceback the optimal alignment
    alignA = ""
    alignB = ""
    alignC = ""
    # Start at the end
    i = lenA
    j = lenB
    k = lenC
    # Stop when we hit the end of all but one sequence
    while (i>0 and j>0) or (j>0 and k>0) or (i>0 and k>0):
        if traceback[i, j, k] == 0: #"matchABC":
            alignA = seqA[i-1] + alignA
            alignB = seqB[j-1] + alignB
            alignC = seqC[k-1] + alignC
            i -= 1
            j -= 1
            k -= 1
        elif traceback[i, j, k] == 3: #"matchAB":
            alignA = seqA[i-1] + alignA
            alignB = seqB[j-1] + alignB
            alignC = "-" + alignC
            i -= 1
            j -= 1
        elif traceback[i, j, k] == 2: #"matchAC":
            alignA = seqA[i-1] + alignA
            alignB = "-" + alignB
            alignC = seqC[k-1] + alignC
            i -= 1
            k -= 1
        elif traceback[i, j, k] == 1: #"matchBC":
            alignA = "-" + alignA
            alignB = seqB[j-1] + alignB
            alignC = seqC[k-1] + alignC
            j -= 1
            k -= 1
        elif traceback[i, j, k] == 4: #"gapAB":
            alignA = "-" + alignA
            alignB = "-" + alignB
            alignC = seqC[k-1] + alignC
            k -= 1
        elif traceback[i, j, k] == 6: #"gapAC":
            alignA = "-" + alignA
            alignB = seqB[j-1] + alignB
            alignC = "-" + alignC
            j -= 1
        elif traceback[i, j, k] == 5: #"gapBC":
            alignA = seqA[i-1] + alignA
            alignB = "-" + alignB
            alignC = "-" + alignC
            i -= 1
    # Fill in the rest of the alignment if it begins with gaps
    # (i.e., traceback all the way to S[0, 0, 0])
    while i > 0:
        alignA = seqA[i-1] + alignA
        alignB = "-" + alignB
        alignC = "-" + alignC
        i -= 1
    while j > 0:
        alignA = "-" + alignA
        alignB = seqB[j-1] + alignB
        alignC = "-" + alignC
        j -= 1
    while k > 0:
        alignA = "-" + alignA
        alignB = "-" + alignB
        alignC = seqC[k-1] + alignC
        k -= 1

    return Alignment([Sequence(alignA, seqA.alphabet, seqA.name, gappy = True),
                      Sequence(alignB, seqB.alphabet, seqB.name, gappy = True),
                      Sequence(alignC, seqC.alphabet, seqC.name, gappy = True)])

def readClustal(string, alphabet):
    """ Read a ClustalW2 alignment in the given string and return as an
    Alignment object. """
    seqs = {} # sequence data
    for line in string.splitlines():
        if line.startswith('CLUSTAL') or line.startswith('STOCKHOLM') \
           or line.startswith('#'):
            continue
        if len(line.strip()) == 0:
            continue
        if line[0] == ' ' or '*' in line or ':' in line:
            continue
        sections = line.split()
        name, seqstr = sections[0:2]
        index = name.find('/')
        if index >= 0:
            name = name[0:index]
        if name in seqs:
            seqs[name] += seqstr
        else:
            seqs[name] = seqstr
    sequences = []
    for name, seqstr in list(seqs.items()):
        sequences.append(Sequence(seqstr, alphabet, name, gappy = True))
    return Alignment(sequences)

def readClustalFile(filename, alphabet):
    """ Read a ClustalW2 alignment file and return an Alignment object
    containing the alignment. """
    fh = open(filename)
    data = fh.read()
    fh.close()
    aln = readClustal(data, alphabet)
    return aln

# Substitution Matrix ------------------

class SubstMatrix():

    scoremat = None
    alphabet = None

    def __init__(self, alphabet):
        self.alphabet = alphabet
        self.scoremat = {}

    def setScores(self, scoremat):
        """ Set all scores in one go.
            scoremat is a (sym1, sym2)-keyed dictionary of scores. """
        self.scoremat = scoremat

    def _getkey(self, sym1, sym2):
        """ Construct canonical (unordered) key for two symbols """
        if sym1 <= sym2:
            return tuple([sym1, sym2])
        else:
            return tuple([sym2, sym1])

    def set(self, sym1, sym2, score):
        """ Add a score to the substitution matrix """
        self.scoremat[self._getkey(sym1, sym2)] = score

    def get(self, sym1, sym2):
        return self.scoremat[self._getkey(sym1, sym2)]

    def __str__(self):
        symbols = self.alphabet.symbols # what symbols are in the alphabet
        i = len(symbols)
        string = ''
        for a in symbols:
            string += a + ' '
            for b in symbols[:len(symbols)-i+1]:
                score = self.scoremat[self._getkey(a, b)]
                if score != None:
                    string += str(score).rjust(3) + ' '
                else:
                    string += "?".rjust(3) + ' '
            string += '\n'
            i -= 1
        string += '    ' + '   '.join(self.alphabet.symbols)
        return string

    def writeFile(self, filename):
        """ Write this substitution matrix to the given file. """
        fh = open(filename, 'w')
        file = ''
        for key in self.scoremat:
            file += ''.join(key) + ': ' + str(self.scoremat[key]) + '\n'
        fh.write(file)
        fh.close()


def readSubstMatrix(filename, alphabet):
    """ Read in the substitution matrix stored in the given file. """
    mat = SubstMatrix(alphabet)
    fh = open(filename, 'r')
    data = fh.read()
    fh.close()
    lines = data.splitlines()
    for line in lines:
        if len(line.strip()) == 0:
            continue
        symbols, score = line.split(':')
        score = int(score)
        mat.set(symbols[0], symbols[1], score)
    return mat

#import os
#os.chdir('/Users/mikael/workspace/binf/data/')  # set to the directory where you keep your files
#BLOSUM62 = readSubstMatrix('blosum62.matrix', Protein_Alphabet)

# Motifs -------------------

class Regexp(object):

    """ A class that defines a sequence pattern in terms of a
    given regular expression, with . indicating any symbol and square brackets
    indicating a selection. See standard regexp definitions for more. """

    def __init__(self, pattern):
        """ Create a new consensus sequence with the given pattern. """
        try:
            self.pattern = pattern
            self.regex = re.compile(pattern)
        except:
            raise RuntimeError('invalid consensus sequence given: %s' % pattern)

    def __str__(self):
        return self.pattern

    def search(self, sequence):
        """ Find matches to the motif in the specified sequence. Returns a list
        of triples, of the form (position, matched string, score). Note that
        the score is always 1.0 because a consensus sequence either matches
        or doesn't. """
        if not type(sequence) is Sequence:
            sequence = Sequence(sequence)
        sequenceString = sequence[:]

        results = []
        for match in self.regex.finditer(sequenceString):
            results.append((match.start(), match.group(), 1.0))
        return results


class PWM(object):

    """ A position weight matrix. """

    def __init__(self, foreground, background = None, start = 0, end = None, pseudo = 0.0):
        """ Create a new PWM from the given probability matrix/ces.
        foreground: can be either an Alignment, a list of Distrib's or an instance of IndepJoint.
        background: must be a Distrib instance or None (in which case a uniform background will be used)
        Specify only a section of the matrix to use with start and end. """
        if isinstance(foreground, Alignment):
            foreground = foreground.getProfile(pseudo = pseudo)
        if isinstance(foreground, IndepJoint):
            foreground = foreground.store
        self.start = start
        self.end = end or len(foreground)
        self.length = self.end - self.start
        self.alphabet = foreground[self.start].alpha
        if False in [ col.alpha == self.alphabet for col in foreground[self.start + 1 : self.end] ]:
            raise RuntimeError("All positions need to be based on the same alphabet")
        self.symbols = self.alphabet.symbols
        # Set foreground probabilities from given alignment
        self.m = numpy.zeros((len(self.symbols), self.length))
        self.fg = foreground[self.start:self.end]
        self.bg = background or Distrib(self.alphabet, 1.0) # specified background or uniform
        if not self.alphabet == self.bg.alpha:
            raise RuntimeError("Background needs to use the same alphabet as the foreground")
        p = self.bg.prob()
        for i in range(self.length):
            q = self.fg[i].prob()
            for j in range(len(self.alphabet)):
                self.m[j][i] = self.logme(q[j], p[j])

    def __len__(self):
        return self.length

    def getRC(self, swap = [('A', 'T'), ('C', 'G')] ):
        """ Get the reverse complement of the current PWM.
            Use for DNA sequences with default params.
        """
        new_fg = self.fg[::-1]  # backwards
        for s in swap:
            new_fg = [d.swapxcopy(s[0], s[1]) for d in new_fg]
        return PWM(new_fg, self.bg)

    MIN_VALUE = 0.00000000001

    def logme(self, fg, bg):
        if fg > self.MIN_VALUE and bg > self.MIN_VALUE:
            ratio = fg / bg
            return math.log(ratio)
        # if not, one of fg and bg is practically zero
        if fg > self.MIN_VALUE: # bg is zero
            return math.log(fg / self.MIN_VALUE)
        else: # fg is zero
            return math.log(self.MIN_VALUE)

    def getMatrix(self):
        return self.m

    def __str__(self):
        str = ''
        for j in range(len(self.alphabet)):
            str += "%s\t%s\n" % (self.alphabet[j], ' '.join("%+6.2f" % (y) for y in self.m[j]))
        return str

    def display(self, format = 'COLUMN'):
        if format == 'COLUMN':
            print((" \t%s" % (' '.join(" %5d" % (i + 1) for i in range(self.length)))))
            for j in range(len(self.alphabet)):
                print(("%s\t%s" % (self.alphabet[j], ' '.join("%+6.2f" % (y) for y in self.m[j]))))
        elif format == 'JASPAR':
            for j in range(len(self.alphabet)):
                print(("%s\t[%s]" % (self.alphabet[j], ' '.join("%+6.2f" % (y) for y in self.m[j]))))

    def search(self, sequence, lowerBound=0):
        """ Find matches to the motif in a specified sequence. Returns a list
        of  results as triples: (position, matched string, score).
        The optional argument lowerBound specifies a lower bound on reported
        scores. """
        results = []
        for i in range(len(sequence)-self.length+1):
            subseq = sequence[i:i + self.length]
            ndxseq = [ self.alphabet.index(sym) for sym in subseq ]
            score = 0.0
            for w in range(len(ndxseq)):
                score += self.m[ ndxseq[w] ][ w ]
            if score > lowerBound:
                results.append((i, subseq, score))
        return results

    def maxscore(self, sequence):
        """ Find matches to the motif in a specified sequence.
            Returns the maximum score found in the sequence and its index as a tuple:
            (maxscore, maxindex) """
        maxscore = None
        maxindex = None
        for i in range(len(sequence)-self.length+1):
            subseq = sequence[i:i + self.length]
            ndxseq = [ self.alphabet.index(sym) for sym in subseq ]
            score = 0.0
            for w in range(len(ndxseq)):
                score += self.m[ ndxseq[w] ][ w ]
            if maxscore == None:
                maxscore = score
                maxindex = i
            elif maxscore < score:
                maxscore = score
                maxindex = i
        return (maxscore, maxindex)

# Web Service Functions -------------------

def getSequence(id, database = 'uniprotkb', start=None, end=None):
    """ Get the sequence identified by the given ID from the given database
    (e.g. 'uniprotkb', 'refseqn' or 'refseqp'), and return it as a Sequence
    object. An error is caused if the sequence ID is not found. If start and
    end are given, then only that section of the sequence is returned.
    Note: more flexible search options are supported by using webservice.fetch
    directly."""

    MAX_TRY = 5

    for i in range(MAX_TRY):
        try:
            fastaData = fetch(id, database)
            seq = readFasta(fastaData)[0]
            break
        except:
            from time import sleep
            print(('Failed on {i}th try for id {id}'.format(i=i, id=id)))
            sleep(0.1)
    try:
        return Sequence(seq[start:end], seq.alphabet, seq.name, seq.info)
    except:
        raise RuntimeError('An error occurred while retrieving the specified sequence: %s (maybe the ID doesn\'t exist)' % id)

def searchSequences(query, database='uniprot'):
    """ Search for sequences matching the given query in the given database
    (must be 'uniprot'), and return a list of sequence IDs. """
    ids = search(query, limit = None)
    return ids

def runClustal(sequences, method='slow'):
    """ Run a ClustalOmega alignment of the given list of Sequence objects.
    Return an Alignment object. Method should be one of 'fast' or 'slow'. """
    alpha = None
    for seq in sequences:
        if alpha == None:
            alpha = seq.alphabet
        elif alpha != seq.alphabet:
            raise RuntimeError("Invalid alphabet: " + str(seq.alphabet) + ". Not compatible with " + str(alpha))
    serviceName = 'clustalo'
    resultType = 'aln-clustal'
    fastaSeqs = ''.join([seq.writeFasta() for seq in sequences])
    params = {'alignment': method.lower(), 'sequence': fastaSeqs}
    service = EBI(serviceName)
    result = service.submit(params, resultType)
    alignment = readClustal(result, alpha)
    return alignment

def createTree(alignment, type):
    """ Run a ClustalW 2 phylogeny tree creation of either a 'Neighbour-joining'
    or 'UPGMA' type tree from the given multiple sequence Alignment object. """
    if not type in ['Neighbour-joining', 'UPGMA']:
        raise RuntimeError('type must be either \'Neighbour-joining\' or \'UPGMA\'.')
    serviceName = 'clustalw2_phylogeny'
    resultType = 'tree'
    output = 'dist'
    clustalAln = alignment.writeClustal()
    params = {'tree': output, 'sequence': clustalAln, 'clustering': type, 'tossgaps': 'true'}
    service = EBI(serviceName)
    tree = service.submit(params, resultType)
    return tree

def runBLAST(sequence, program='blastp', database='uniprotkb', exp='1e-1'):
    """ Run a BLAST search of nucleotide mouse databases using the given
    sequence as a query. Return a list of matched sequence IDs, in descending
    order of similarity to query sequence.
    program: either blastn (nucleotide) or blastp (protein)
    database: many available, e.g. uniprotkb, pdb (protein); em_rel, nrnl1 (EMBL nucleotide, non-redundant resp)
        (for protein see http://www.ebi.ac.uk/Tools/sss/ncbiblast/help/index-protein.html#database)
        (for nucleotide see http://www.ebi.ac.uk/Tools/sss/ncbiblast/help/index-nucleotide.html#database)
    exp: E-value threshold (select only hits that have a better E-value than this)
    """
    if sequence.alphabet == predefAlphabets['DNA']:
        stype = 'dna'
    elif sequence.alphabet == predefAlphabets['RNA']:
        stype = 'rna'
    else:
        stype = 'protein'
    serviceName = 'ncbiblast'
    resultTypes = ['ids', 'out'] # request
    fastaSeq = sequence.writeFasta()
    databases = [database]
    params = {'program': program, 'database': databases, 'sequence': fastaSeq,
              'stype': stype, 'exp': exp}
    service = EBI(serviceName)
    idsData, output = service.submit(params, resultTypes)
    ids=[]
    for id in idsData.splitlines():
        if len(id) > 0:
            ids.append(id.split(':')[1])
    return ids

if __name__ == '__main__':
    seqs = readFastaFile('/Users/mikael/ASR/CYP11/CYP11_aln_full.fa', Protein_wX, gappy=True)
    print(('Read', len(seqs), 'sequences'))
