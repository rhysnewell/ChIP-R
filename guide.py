###################################################
# This module is a supplement to the Python guide #
# Version 2017.3  (10/03/2017)                    #
###################################################
'''
This module contains code that can help solve bioinformatics problems.
See the accompanying Python guide document for more explanations and examples.

Alphabet is a class that defines valid symbols that we then use to make up valid
biological sequences. Note that we also define variables corresponding to
DNA, RNA and Protein sequences that can be used directly.

Sequence is a class that defines basic parts and operations on biological sequences.

Alignment is a class that defines an alignment of sequences (how symbols in different sequences line
up when placed on-top of one another). Alignment methods should generate instances of this class.

SubstMatrix is a class that defines a substitution matrix, i.e. a scoring system for performing
alignments. You can read these from files or construct them manually.

GeneProfile is a class that defines parts and operations for gene expression profiles. Essentially,
the class will help to index expression data by gene name (rows) and by sample name (columns).

There are several methods not tied to a particular class because they construct new instances,
e.g. reading from file, retrieving from the internet, creating an alignment from sequences etc.

You need to have numpy installed (see http://www.numpy.org/).
Should work with Python v3.5 (see http://www.python.org/).
The code may contain bugs--please report to m.boden@uq.edu.au
'''

import math, numpy, urllib.request, urllib.parse, urllib.error, urllib.request, urllib.error, urllib.parse

###############################################################################
# Alphabet                                                                    #
###############################################################################

class Alphabet():
    """ A minimal class for alphabets
        Alphabets include DNA, RNA and Protein """
    def __init__(self, symbolString):
        self.symbols = symbolString
    def __len__(self):              # implements the "len" operator, e.g. "len(Alphabet('XYZ'))" results in 3
        return len(self.symbols)    # will tell you the length of the symbols in an Alphabet instance
    def __contains__(self, sym):    # implements the "in" operator, e.g. "'A' in Alphabet('ACGT')" results in True
        return sym in self.symbols  # will tell you if 'A' is in the symbols in an Alphabet instance
    def __iter__(self):             # method that allows us to iterate over all symbols, e.g. "for sym in Alphabet('ACGT'): print sym" prints A, C, G and T on separate lines
        tsyms = tuple(self.symbols)
        return tsyms.__iter__()
    def __getitem__(self, ndx):
        """ Retrieve the symbol(s) at the specified index (or slice of indices) """
        return self.symbols[ndx]
    def index(self, sym):
        """ Retrieve the index of the given symbol in the alphabet. """
        return self.symbols.index(sym)
    def __str__(self):
        return self.symbols

""" Below we declare alphabet variables that are going to be available when
this module (this .py file) is imported """
DNA_Alphabet = Alphabet('ACGT')
RNA_Alphabet = Alphabet('ACGU')
Protein_Alphabet = Alphabet('ACDEFGHIKLMNPQRSTVWY')
Protein_wX = Alphabet('ACDEFGHIKLMNPQRSTVWYX')
Protein_wGAP = Alphabet('ACDEFGHIKLMNPQRSTVWY-')

###############################################################################
# Sequence                                                                    #
###############################################################################

class Sequence():
    """ A biological sequence class. Stores the sequence itself,
        the alphabet and a name.
        Usage:
        Create an instance of Sequence - have to pass in sequence, alphabet and name - gappy is an optional argument
        >>> seq1 = Sequence('ACGGGAGAGG', DNA_Alphabet, 'ABC')
        >>> print(seq1)
        ABC: ACGGGAGAGG
        >>> 'C' in seq1
        True
        >>> for sym in seq1:
        ...     print(sym)
        """
    def __init__(self, sequence, alphabet, name = '', gappy = False, annot = ''):
        """ Construct a sequence from a string, an alphabet (gappy or not) and a name.
            The parameter gappy is for sequences when used in alignments. """
        for sym in sequence:
            if not sym in alphabet and (sym != '-' or not gappy):  # error check: bail out
                raise RuntimeError('Invalid symbol: ' + sym)
        self.sequence = sequence # Store sequence
        self.alphabet = alphabet # Store alphabet
        self.name = name         # Store name
        self.gappy = gappy
        self.annot = annot # some annotation, e.g. species

    def __len__(self):      # the "len" operator
        return len(self.sequence)
    def __iter__(self):     # method that allows us to iterate over a sequence
        tsyms = tuple(self.sequence)
        return tsyms.__iter__()
    def __contains__(self, item):   # test for membership (the "in" operator)
        for sym in self.sequence:
            if sym == item:
                return True
        return False
    def __getitem__(self, ndx):     # [ndx] operator (retrieve a specified index (or a "slice" of indices) of the sequence data.
        return self.sequence[ndx]
    def writeFasta(self):
        """ Write one sequence in FASTA format to a string and return it. """
        fasta = '>' + self.name + ' ' + self.annot + '\n'
        data = self.sequence
        nlines = (len(self.sequence) - 1) // 60 + 1
        for i in range(nlines):
            lineofseq = ''.join(data[i*60 : (i+1)*60]) + '\n'
            fasta += lineofseq
        return fasta
    def __str__(self):      # "pretty" print sequence
        str = self.name + ': '
        for sym in self:
            str += sym
        return str
    def count(self, findme):
        """ Get the number of occurrences of specified symbol """
        cnt = 0
        for sym in self.sequence:
            if findme == sym:
                cnt = cnt + 1
        return cnt
    def find(self, findme):
        """ Find the position of the specified symbol or sub-sequence """
        return self.sequence.find(findme)

###############################################################################
# Alignment                                                                   #
###############################################################################

class Alignment():
    """ A sequence alignment class. Stores two or more sequences of equal length where
    one symbol is gap '-'. The number of columns in the alignment is given by alignlen.
    Example usage:
    >>> seqs = [Sequence('THIS-LI-NE', Protein_Alphabet, gappy = True), Sequence('--ISALIGNED', Protein_Alphabet, gappy = True)]
    >>> print(Alignment(seqs))
     THIS-LI-NE-
     --ISALIGNED """
    def __init__(self, seqs):
        self.alphabet = None
        self.alignlen = -1
        self.seqs = seqs
        self.namelen = 0
        for s in seqs:
            if self.alphabet == None:
                self.alphabet = s.alphabet
            elif self.alphabet != s.alphabet:
                raise RuntimeError("Alignment invalid: contains a mix of alphabets")
            if self.alignlen == -1:
                self.alignlen = len(s)
            elif self.alignlen != len(s):
                raise RuntimeError("Alignment invalid: lengths vary")
            self.namelen = max(len(s.name), self.namelen)
    def __str__(self):
        string = u''
        for seq in self.seqs:
            string += seq.name.ljust(self.namelen+1)
            for sym in seq:
                string += sym
            string += '\n'
        return string
    def __len__(self):
        """ Defines what the "len" operator returns for an instance of Alignment: the number of sequences. """
        return len(self.seqs)
    def __getitem__(self, ndx):
        return self.seqs[ndx]
    def calcDistances(self, measure = 'fractional', a=1.0):
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
                p = D / L
                # Now calculate the specified measure based on p
                if measure == 'fractional':
                    dist = p
                else:
                    raise RuntimeError('Not implemented: %s' % measure)
                distmat[i, j] = distmat[j, i] = dist
        return distmat
    def writeClustal(self):
        """ Write the alignment to a string using the Clustal file format. """
        symbolsPerLine = 60
        maxNameLength =  self.namelen + 1
        mystring = u''
        wholeRows = self.alignlen // symbolsPerLine
        for i in range(wholeRows):
            for j in range(len(self.seqs)):
                mystring += self.seqs[j].name.ljust(maxNameLength) + u' '
                mystring += self.seqs[j][i*symbolsPerLine:(i+1)*symbolsPerLine] + u'\n'
            mystring += u'\n'
        # Possible last row
        lastRowLength = self.alignlen - wholeRows*symbolsPerLine
        if lastRowLength > 0:
            for j in range(len(self.seqs)):
                if maxNameLength > 0:
                    mystring += self.seqs[j].name.ljust(maxNameLength) + u' '
                mystring += self.seqs[j][-lastRowLength:] + '\n'
        return mystring
    def writeHTML(self, filename):
        """ Generate HTML that displays the alignment in colour. """
        fh = open(filename, 'wt')
        fh.write(u'<html><head><meta content="text/html; charset=ISO-8859-1" http-equiv="Content-Type">\n<title>Sequence Alignment</title>\n</head><body><pre>\n')
        html = u''.ljust(self.namelen) + u' '
        for i in range(self.alignlen - 1):
            if (i+1) % 10 == 0:
                html += str(i//10+1)[-1]
            else:
                html += ' '
        html += u'%s\n' % (self.alignlen)
        fh.write(html)
        if self.alignlen > 10:
            html = u''.ljust(self.namelen) + u' '
            for i in range(self.alignlen - 1):
                if (i+1) % 10 == 0:
                    html += u'0'
                else:
                    html += u' '
            html += u'\n'
            fh.write(html)
        if len(self.alphabet) <= 5: # DNA or RNA
            colours = {'A':u'green','C':u'orange','G':u'red','T':u'#66bbff','U':u'#66bbff'}
        else: # amino acids
            colours = {'G':u'orange','P':u'orange','S':u'orange','T':u'orange','H':u'red','K':u'red','R':u'red','F':u'#66bbff','Y':u'#66bbff','W':u'#66bbff','I':u'green','L':u'green','M':u'green','V':u'green'}
        for seq in self.seqs:
            html = seq.name.ljust(self.namelen) + u' '
            for sym in seq:
                try:
                    colour = colours[sym]
                except KeyError:
                    colour = u'white'
                html += u'<font style="BACKGROUND-COLOR: %s">%s</font>' % (colour, sym)
            html += u'\n'
            fh.write(html)
        fh.write(u'</pre></body></html>\n')
        fh.close()

def scoreAlignment(aln, substmat = None, gap = -1):
    """Score an alignment (aln) using a substitution matrix (substmat).
       If the alignment consists of more than two sequences, the minimum
       score of each column is used.
       If substmat is not specified (None), the count of matches is returned.
    """
    nseqs = len(aln.seqs)
    total = 0
    for pos in range(aln.alignlen):
        min = None
        for i in range(nseqs):
            for j in range(i+1, nseqs):
                gap_here = aln.seqs[i][pos] == '-' or aln.seqs[j][pos] == '-'
                score = 0
                if substmat == None:
                    if aln.seqs[i][pos] == aln.seqs[j][pos]:
                        score = 1
                else: # we have a substitution matrix
                    if gap_here:
                        score = gap
                    else:
                        score = substmat.get(aln.seqs[i][pos], aln.seqs[j][pos])
                if min == None:
                    min = score
                elif min > score:
                    min = score
        total += min
    return total

###############################################################################
# Methods to create instances of Alignment                                    #
###############################################################################

def align(seqA, seqB, substMatrix, gap=-1):
    """ Align seqA with seqB using the Needleman-Wunsch
    (global) algorithm. substMatrix is the substitution matrix to use and
    gap is the linear gap penalty to use. """
    stringA, stringB = seqA.sequence, seqB.sequence
    lenA, lenB = len(seqA), len(seqB)
    # Create the scoring matrix (S) and a matrix for traceback
    S = numpy.zeros((lenA + 1, lenB + 1))
    Traceback = numpy.zeros((lenA + 1, lenB + 1))
    # Fill the first row and column of S with multiples of the gap penalty
    for i in range(lenA + 1):
        S[i, 0] = i * gap
    for j in range(lenB + 1):
        S[0, j] = j * gap
    # Calculate the optimum score at each location in the matrix, note which option that was chosen for traceback
    for i in range(1, lenA + 1):
        for j in range(1, lenB + 1):
            match = S[i - 1, j - 1] + substMatrix.get(stringA[i - 1], stringB[j - 1])
            delete = S[i - 1, j] + gap
            insert = S[i, j - 1] + gap
            Traceback[i, j] = numpy.argmax([match, delete, insert])
            S[i, j] = max([match, delete, insert])
    # Trace back the optimal alignment
    alignA = ''
    alignB = ''
    # Start at the end
    i = lenA
    j = lenB
    # Stop when we hit the end of a sequence
    while i > 0 and j > 0:
        if Traceback[i, j] == 1:
            # Got here by a gap in sequence B (go up)
            alignA = stringA[i - 1] + alignA
            alignB = '-' + alignB
            i -= 1
        elif Traceback[i, j] == 2:
            # Got here by a gap in sequence A (go left)
            alignA = "-" + alignA
            alignB = stringB[j - 1] + alignB
            j -= 1
        else:
            # Got here by aligning the bases (go diagonally)
            alignA = stringA[i - 1] + alignA
            alignB = stringB[j - 1] + alignB
            i -= 1
            j -= 1
    # Fill in the rest of the alignment if it begins with gaps
    # (i.e., trace back all the way to S[0, 0])
    while i > 0:
        # Go up
        alignA = stringA[i - 1] + alignA
        alignB = '-' + alignB
        i -= 1
    while j > 0:
        # Go left
        alignA = '-' + alignA
        alignB = stringB[j - 1] + alignB
        j -= 1
    return Alignment([Sequence(alignA, seqA.alphabet, seqA.name, gappy=True),
                      Sequence(alignB, seqB.alphabet, seqB.name, gappy=True)])

###############################################################################
# SubstMatrix                                                                 #
###############################################################################

class SubstMatrix():
    """ Create a substitution matrix for an alphabet.
    Example usage:
    >>> sm = SubstMatrix(DNA_Alphabet)
    >>> for a in DNA_Alphabet:
    ...     for b in DNA_Alphabet:
    ...         if a > b:
    ...             sm.set(a, b, -1)
    ...         elif a == b:
    ...             sm.set(a, b, +1)
    ...
    >>> print(sm)
    A   1
    C  -1   1
    G  -1  -1   1
    T  -1  -1  -1   1
        A   C   G   T
    >>> sm.get('C', 'T')
    -1
    """
    def __init__(self, alphabet, scoremat = None):
        self.scoremat = scoremat or {}      # start with empty dictionary
        self.alphabet = alphabet
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
        string = u''
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
        file = u''
        for key in self.scoremat:
            file += ''.join(key) + ': ' + str(self.scoremat[key]) + '\n'
        fh.write(file)
        fh.close()

###############################################################################
# Below are some useful methods for loading data from strings and files.      #
# They recognize the FASTA and Clustal formats (nothing fancy).               #
###############################################################################

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

def readFastaString(string, alphabet, gappy = False):
    """ Read the given string as FASTA formatted data and return the list of
    sequences contained within it. """
    seqlist = []    # list of sequences contained in the string
    seqname = ''  # name of *current* sequence
    seqannot = '' # annotation of *current* sequence
    seqdata = []    # sequence data for *current* sequence
    for line in string.splitlines():    # read every line
        if len(line) == 0:              # ignore empty lines
            continue
        if line[0] == '>':  # start of new sequence
            if seqname:     # check if we've got one current
                current = Sequence(''.join(seqdata), alphabet, seqname, gappy, seqannot)
                seqlist.append(current)
            # now collect data about the new sequence
            parts = line[1:].split() # skip first char
            seqname = ''  # name of *current* sequence
            seqannot = '' # annotation of *current* sequence
            if len(parts) > 0: seqname = parts[0]
            if len(parts) > 1: seqannot = line[len(seqname) + 2:] # the rest of the line
            seqdata = []
        else:               # we assume this is (more) data for current
            cleanline = line.split()
            for thisline in cleanline:
                seqdata.extend(tuple(thisline.strip('*')))
    # we're done reading the file, but the last sequence remains
    if seqname:
        lastseq = Sequence(''.join(seqdata), alphabet, seqname, gappy, seqannot)
        seqlist.append(lastseq)
    return seqlist

def readFastaFile(filename, alphabet):
    """ Read the given FASTA formatted file and return the list of sequences
    contained within it. """
    fh = open(filename)
    data = fh.read()
    fh.close()
    seqlist = readFastaString(data, alphabet)
    return seqlist

def writeFastaFile(filename, seqs):
    """ Write the specified sequences to a FASTA file. """
    fh = open(filename, 'wt')
    for seq in seqs:
        fh.write(seq.writeFasta())
    fh.close()

def readClustalString(string, alphabet):
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
        name, seq = sections[0:2]
        if name in seqs:
            seqs[name] += seq
        else:
            seqs[name] = seq
    sequences = []
    for name, seq in seqs.items():
        sequences.append(Sequence(seq, alphabet, name, gappy = True))
    return Alignment(sequences)

def readClustalFile(filename, alphabet):
    """ Read a ClustalW2 alignment file and return an Alignment object
    containing the alignment. """
    fh = open(filename)
    data = fh.read()
    fh.close()
    aln = readClustalString(data, alphabet)
    return aln

def writeClustalFile(filename, aln):
    """ Write the specified alignment to a Clustal file. """
    fh = open(filename, 'wt')
    fh.write('CLUSTAL W (1.83) multiple sequence alignment\n\n\n') # fake header so that clustal believes it
    fh.write(aln.writeClustal())
    fh.close()

###############################################################################
# GeneProfile                                                                 #
###############################################################################

class GeneProfile():
    """ A class for gene expression data.
    Example usage:
    >>> gp = GeneProfile('MyMicroarray', ['Exp1', 'Exp2'])
    >>> gp['gene1'] = [0.1, 0.5]
    >>> gp['gene2'] = [2, 1]
    >>> gp.getSample('Exp2')
    {'gene1': [0.5], 'gene2': [1.0]}
    """
    def __init__(self, dataset_name='', sample_names=[], profiles = None):
        """ Create a gene profile set. """
        self.name = dataset_name
        self.samples = sample_names
        self.genes = profiles or {} # dictionary for storing all gene--measurement pairs
    def __setitem__(self, name, probevalues):
        if len(probevalues) == len(self.samples):
            self.genes[name] = [float(y) for y in probevalues]
        else:
            raise RuntimeError('Invalid number of measurements for probe ' + name)
    def __getitem__(self, name):
        return self.genes[name]
    def getSorted(self, index, descending=True):
        """Get a list of (gene, value) tuples in descending order by value"""
        key_fn = lambda v: v[1][index]
        return sorted(list(self.genes.items()), key=key_fn, reverse=descending)
    def addSample(self, sample_name, sample_dict):
        """Add a sample to the current data set.
           sample_dict is a dictionary with the same keys as the current gene set.
           Only values for genes in the current set will be added. """
        self.headers.extend(sample_name)
        if not self.genes:
            self.genes = sample_dict
        else:
            for gene in self.genes:
                values = sample_dict[gene]
                if values:
                    self.genes[gene].extend([float(y) for y in values])
                else:
                    self.genes[gene].extend([0.0 for _ in sample_name])
        return self.genes
    def getSample(self, sample_name):
        """Construct a gene dictionary including only named samples. """
        mygenes = {}
        if isinstance(sample_name, str):    # a single sample-name
            mysamples = [sample_name]
        else:                               # a list of sample-names
            mysamples = sample_name
        for gene in self.genes:
            mygenes[gene] = []
            for name in mysamples:
                mygenes[gene].append(self.genes[gene][self.samples.index(name)])
        return mygenes
    def getRatio(self, sample1, sample2):
        """Get the ratio of two samples in the data set. """
        mygenes = {}
        index1 = self.samples.index(sample1)
        index2 = self.samples.index(sample2)
        for gene in self.genes:
            mygenes[gene] = []
            mygenes[gene].append(self.genes[gene][index1] / self.genes[gene][index2])
        return mygenes
    def __str__(self):
        """ Display data as a truncated GEO SOFT file named filename. """
        line = '^DATASET = ' + self.dataset + '\n'
        line += '!dataset_table_begin\nID_REF\t'
        for header in self.headers:
            line += header + '\t'
        line += '\n'
        for gene in self.genes:
            line += gene + '\t'
            values = self.genes[gene]
            for value in values:
                line += format(value, '5.3f') + '\t'
            line += '\n'
        line += '!dataset_table_end\n'
    def writeGeoFile(self, filename):
        fh = open(filename, 'w')
        fh.write(str(self))
        fh.close()

def getLog(genedict, base=2):
    """Get the log-transformed value of a sample/column. """
    mygenes = {}
    for gene in genedict:
        mygenes[gene] = []
        for sample in genedict[gene]:
            mygenes[gene].append(math.log(sample, base))
    return mygenes

def readGeoFile(filename, id_column = 0):
    """ Read a Gene Expression Omnibus SOFT file. """
    dataset = None
    fh = open(filename, "rU")
    manylines = fh.read()
    fh.close()
    data_rows = False  # Indicates whether we're reading the data section or metadata
    name = u'Unknown'
    cnt_data = 0
    for line in manylines.splitlines():
        if line.startswith('^DATASET'):
            name = line.split('= ')[1]
            continue
        data_rows = line.startswith('!dataset_table_begin')
        data_rows = not line.startswith('!dataset_table_end')
        if len(line.strip()) == 0 or line.startswith('!') or line.startswith('#') or line.startswith('^'):
            continue
        if data_rows:
            cnt_data += 1
            if (cnt_data == 1):  # First line contains the headers
                headers = line.split('\t')
                dataset = GeneProfile(name, headers[2:])  # Create the data set
                continue
            ignore = (dataset == None)  # ignore the row if the dataset is not initialised
            id = line.split('\t')[id_column]
            values = []
            cnt_word = 0
            for word in line.split('\t'):
                cnt_word += 1
                if cnt_word <= (id_column + 1): # ignore the gene names
                    continue
                if word == 'null':
                    ignore = True # ignore gene if a value is null
                    continue
                try:
                    values.append(float(word))
                except:  # ignore values that are not "float"
                    continue
            if not ignore:
                dataset[id] = tuple(values)
    print('Data set %s contains %d genes' % (name, len(dataset.genes)))
    return dataset

###############################################################################
# Web service methods that find data in online databases.
# Our implementations are mainly serviced by EBI.
###############################################################################

def getSequence(entryId, dbName = 'uniprotkb', alphabet = Protein_Alphabet, format = 'fasta', debug: bool = True):
    """ Retrieve a single entry from a database
    entryId: ID for entry e.g. 'P63166' or 'SUMO1_MOUSE'
    dbName: name of database e.g. 'uniprotkb' or 'pdb' or 'refseqn'; see http://www.ebi.ac.uk/Tools/dbfetch/dbfetch/dbfetch.databases for available databases
    format: file format specific to database e.g. 'fasta' or 'uniprot' for uniprotkb (see http://www.ebi.ac.uk/Tools/dbfetch/dbfetch/dbfetch.databases)
    See http://www.ebi.ac.uk/Tools/dbfetch/syntax.jsp for more info re URL syntax
    """
    if not isinstance(entryId, str):
        entryId = entryId.decode("utf-8")
    url ='http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?style=raw&db=' + dbName + '&format=' + format + '&id=' + entryId
    try:
        if debug:
            print('DEBUG: Querying URL: {0}'.format(url))
        data = urllib.request.urlopen(url).read()
        if format == 'fasta':
            return readFastaString(data.decode("utf-8"), alphabet)[0]
        else:
            return data.decode("utf-8")
    except urllib.error.HTTPError as ex:
        raise RuntimeError(ex.read())

def searchSequences(query, dbName = 'uniprot'):
    """
    Retrieve multiple entries matching query from a database currently only via UniProtKB
    query: search term(s) e.g. 'organism:9606+AND+antigen'
    dbName: name of database e.g. 'uniprot', "refseq:protein", "refseq:pubmed"
    See http://www.uniprot.org/faq/28 for more info re UniprotKB's URL syntax
    See http://www.ncbi.nlm.nih.gov/books/NBK25499/ for more on NCBI's E-utils
    """
    if dbName.startswith('uniprot'):
        # Construct URL
        url = 'http://www.uniprot.org/uniprot/?format=list&query=' + query
        try:
            data = urllib.request.urlopen(url).read()
            return data.decode("utf-8").splitlines()
        except urllib.error.HTTPError as ex:
            raise RuntimeError(ex.read())
    elif dbName.startswith('refseq'):
        dbs = dbName.split(":")
        if len(dbs) > 1:
            dbName = dbs[1]
        base = u'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
        url = base + "esearch.fcgi?db=" + dbName + "&term=" + query
        # Get the entries
        try:
            data = urllib.request.urlopen(url).read()
            words = data.decode("utf-8").split("</Id>")
            words = [w[w.find("<Id>")+4:] for w in words[:-1]]
            return words
        except urllib.error.HTTPError as ex:
            raise RuntimeError(ex.read())
    return

def idmap(identifiers, frm='ACC', to='P_REFSEQ_AC'):
    """
    Map identifiers between databases (based on UniProtKB;
    see http://www.uniprot.org/help/programmatic_access)
    identifiers: a list of identifiers (list of strings)
    frm: the abbreviation for the identifier FROM which to idmap
    to: the abbreviation for the identifier TO which to idmap
    Returns a dictionary with key (from) -> value (to) """
    url = u'http://www.uniprot.org/uploadlists/'
    # construct query by concatenating the list of identifiers
    if isinstance(identifiers, str):
        query = identifiers.strip()
    else: # assume it is a list of strings
        query = u''
        for id in identifiers:
            query = query + id + ' '
        query = query.strip() # remove trailing spaces
    params = {
        'from' : frm,
        'to' : to,
        'format' : 'tab',
        'query' : query
    }
    if len(query) > 0:
        data = urllib.parse.urlencode(params).encode("utf-8")
        request = urllib.request.Request(url, data)
        response = urllib.request.urlopen(request).read()
        d = dict()
        for row in response.decode("utf-8").splitlines()[1:]:
            pair = row.split('\t')
            d[pair[0]] = pair[1]
        return d
    else:
        return dict()

###############################################################################
# Gene Ontology services.
# See http://www.ebi.ac.uk/QuickGO/WebServices.html for more info
###############################################################################

def getGODef(goterm):
    """
    Retrieve information about a GO term
    goterm: the identifier, e.g. 'GO:0002080'
    """
    # Construct URL
    url = u'http://www.ebi.ac.uk/QuickGO/GTerm?format=obo&id=' + goterm
    # Get the entry: fill in the fields specified below
    try:
        entry={'id': None, 'name': None, 'def': None}
        data = urllib.request.urlopen(url).read()
        for row in data.decode("utf-8").splitlines():
            index = row.find(u':')
            if index > 0 and len(row[index:]) > 1:
                field = row[0:index].strip()
                value = row[index+1:].strip(u' "') # remove spaces
                if field in list(entry.keys()):         # check if we need field
                    if entry[field] == None:      # check if assigned
                        entry[field] = value
        return entry
    except urllib.error.HTTPError as ex:
        raise RuntimeError(ex.read())

def getGOTerms(genes, db='UniProtKB'):
    """
    Retrieve all GO terms for a given set of genes (or single gene).
    db: use specified database, e.g. 'UniProtKB', 'UniGene',
    or 'Ensembl'.
    The result is given as a map (key=gene name, value=list of unique
    terms) OR in the case of a single gene as a list of unique terms.
    """
    if type(genes) != list and type(genes) != set and type(genes) != tuple:
        genes = [genes]  # if 'genes' is a single gene, we make a single item list
    map = dict()
    uri = u'http://www.ebi.ac.uk/QuickGO/GAnnotation?format=tsv&db=' + db + u'&protein='
    for gene in genes:
        terms = set()  # empty result set
        url = uri + gene # Construct URL
        try: # Get the entry: fill in the fields specified below
            data = urllib.request.urlopen(url).read()
            for row in data.decode("utf-8").splitlines()[1:]:  # we ignore header row
                values = row.split('\t')
                if len(values) >= 7:
                    terms.add(values[6]) # add term to result set
            map[gene] = list(terms)      # make a list of the set
        except urllib.error.HTTPError as ex:
            raise RuntimeError(ex.read())
    if len(genes) == 1:
        return map[genes[0]]
    else:
        return map

def getGenes(goterms, db='UniProtKB', taxo=None):
    """
    Retrieve all genes/proteins for a given set of GO terms
    (or single GO term).
    db: use specified database, e.g. 'UniProtKB', 'UniGene',
    or 'Ensembl'
    taxo: use specific taxonomic identifier, e.g. 9606 (human)
    The result is given as a map (key=gene name, value=list of unique
    terms) OR in the case of a single gene as a list of unique terms.
    """
    if type(goterms) != list and type(goterms) != set and type(goterms) != tuple:
        goterms = [goterms]
    map = dict()
    if taxo == None:
        uri = u'http://www.ebi.ac.uk/QuickGO/GAnnotation?format=tsv&db=' + db + u'&term='
    else:
        uri = u'http://www.ebi.ac.uk/QuickGO/GAnnotation?format=tsv&db=' + db + u'&tax=' + str(taxo) + u'&term='
    for goterm in goterms:
        genes = set()   # start with empty result set
        url = uri + goterm # Construct URL
        try: # Get the entry: fill in the fields specified below
            data = urllib.request.urlopen(url).read()
            for row in data.decode("utf-8").splitlines()[1:]:  # we ignore first (header) row
                values = row.split('\t')
                if len(values) >= 7:
                    genes.add(values[1])  # add gene name to result set
            map[goterm] = list(genes)
        except urllib.error.HTTPError as ex:
            raise RuntimeError(ex.read())
    if len(goterms) == 1:
        return map[goterms[0]]
    else:
        return map

###############################################################################
# PhyloTree                                                                   #
###############################################################################

class PhyloTree:
    """ Rooted, binary (bifurcating) tree for representing phylogenetic relationships.
        Functionality includes labelling and traversing nodes; reading and writing to Newick format;
        association with sequence alignment; maximum parsimony inference of ancestral sequence;
        generation of single, bifurcating rooted tree by UPGMA.
        Known issues: Binary only; Parsimony does not handle gaps in alignment.
        Programmers should note that almost all functionality is implemented through recursion. """

    def __init__(self, root):
        """ Create a tree from a node that is "root" in the tree."""
        self.root = root

    def putAlignment(self, aln):
        """ Associate the tree with a set of sequences/alignment.
            Involves assigning the sequence to the leaf nodes. """
        self.aln = aln
        self.root._assignAlignment(aln)

    def __str__(self):
        """ Produce a printable representation of the tree, specifically the root of the tree. """
        return str(self.root)

    def strSequences(self, start=None, end=None):
        """ Produce a sequence representation of the tree, specifically the root of the tree.
            Specify the start and end positions in the alignment for the sequence to be printed
            (if None the min and max positions will be used). """
        if self.aln != None:
            my_start = start or 0
            my_end = end or self.aln.alignlen
            return self.root._printSequences(my_start, my_end)

    def findLabel(self, label):
        """ Retrieve/return the node with the specified label.
            Returns None if not found."""
        return self.root._findLabel(label)

    def getDescendantsOf(self, node, transitive=False):
        """ Retrieve and return the (list of) descendants (children) of a specified node.
            Node can be the label or the instance.
            transitive indicates if only the direct descendants (False) or if all descendants
            should be returned.
            If node does not exist, None is returned.
            If node has no descendants, an empty list will be returned."""
        if not isinstance(node, PhyloNode):
            node = self.root.findLabel(node)
        if node:
            return node.getDescendants(transitive)
        return None

    def getAncestorsOf(self, node, transitive=False):
        """ Retrieve and return the ancestor (transitive=False) or
            ancestors (transitive=True) of a specified node.
            Node can be the label or the instance.
            If node does not exist, None is returned.
            If node is the root of the tree, None is returned."""
        if not isinstance(node, PhyloNode):
            node = self.root.findLabel(node)
        if node:
            myroot = self.root
            found = False
            branching = []
            while not found and myroot != None:
                branching.append(myroot)
                if myroot.left == node or myroot.right == node:
                    found = True
                    break
                if myroot.left:
                    if myroot.left.isAncestorOf(node, transitive=True):
                        myroot = myroot.left
                    else:  # must be right branch then...
                        myroot = myroot.right
                else:  # must be right branch then...
                    myroot = myroot.right
            if found and transitive:
                return branching
            elif found and len(branching) > 0:
                return branching[len(branching) - 1]
            return None

    def parsimony(self):
        """ Solve the "small parsimony problem",
            i.e. find the sequences on each of the internal nodes.
            See Jones and Pevzner, p. 368 and onwards, for details. """
        self.root._forwardParsimony(self.aln)  # setup and compute scores for all nodes
        self.root._backwardParsimony(self.aln)  # use scores to determine sequences
        return self.root.getSequence()  # return the sequence found at the root

###############################################################################
# PhyloNode                                                                   #
###############################################################################

class PhyloNode:
    """ A class for a node in a rooted, binary (bifurcating) tree.
        Contains pointers to descendants/daughters (left and right),
        optional fields include data, label, sequence and dist.
        If parsimony is used scores and traceback pointers are available.
        A number of methods are named with a _ prefix. These can be, but
        are not intended to be used from outside the class. """

    def __init__(self, label=''):
        """ Initialise an initially unlinked node.
            Populate fields left and right to link it with other nodes.
            Set label to name it.
            Use field data for any type of information associated with node.
            Use dist to indicate the distance to its parent (if any).
            Other fields are used internally, including sequence for associated alignment,
            seqscores, backleft and backright for maximum parsimony. """
        self.left = None
        self.right = None
        self.data = None
        self.label = label
        self.dist = None
        self.sequence = None  # The sequence after an alignment have been mapped (leaf) or the most parsimonous sequence (ancestral)
        self.seqscores = None  # The scores propagated from leaves via children
        self.backleft = None  # Pointers back to left child: what symbol rendered current/parent symbols
        self.backright = None  # Pointers back to right child: what symbol rendered current/parent symbols

    def __str__(self):
        """ Returns string with node (incl descendants) in a Newick style. """
        left = right = label = dist = ''
        if self.left:
            left = str(self.left)
        if self.right:
            right = str(self.right)
        if self.dist or self.dist == 0.0:
            dist = ':' + str(self.dist)
        if self.label != None:
            label = str(self.label)
            if not self.left and not self.right:
                return label + dist
            else:
                return '(' + left + ',' + right + ')' + label + dist
        else:  # there is no label
            if not self.left and self.right:
                return ',' + right
            elif self.left and not self.right:
                return left + ','
            elif self.left and self.right:
                return '(' + left + ',' + right + ')' + dist

    def _printSequences(self, start, end):
        """ Returns string with node (incl descendants) in a Newick style. """
        left = right = label = dist = ''
        if self.left:
            left = self.left._printSequences(start, end)
        if self.right:
            right = self.right._printSequences(start, end)
        if self.dist:
            dist = ':' + str(self.dist)
        if self.sequence != None:
            label = "".join(self.sequence[start:end]) + ""
            if not self.left and not self.right:
                return label + dist
            else:
                return '(' + left + ',' + right + ')' + label + dist
        else:  # there is no label
            if not self.left and self.right:
                return ',' + right
            elif self.left and not self.right:
                return left + ','
            elif self.left and self.right:
                return '(' + left + ',' + right + ')' + dist

    def _findLabel(self, label):
        """ Find a node by label at this node or in any descendants (recursively). """
        if self.label == label:
            return self
        else:
            if self.left:
                foundLeft = self.left._findLabel(label)
                if foundLeft:
                    return foundLeft
            if self.right:
                return self.right._findLabel(label)
            return None

    def _propagateDistance(self, parent_dist):
        """ Convert absolute distances to relative.
            The only parameter is the absolute distance to the parent of this node. """
        travelled = self.dist  # absolute distance to this node
        self.dist = parent_dist - self.dist  # relative distance to this node
        if self.left != None:  # if there is a child node...
            self.left._propagateDistance(travelled)  # pass absolute distance to this node
        if self.right != None:
            self.right._propagateDistance(travelled)

    def _assignAlignment(self, aln):
        """ Assign an alignment to the node, which implies assigning a sequence to it if one is
            available in the alignment. """
        self.sequence = None
        if self.left != None:
            self.left._assignAlignment(aln)
        if self.right != None:
            self.right._assignAlignment(aln)
        for seq in aln.seqs:
            if seq.name == self.label:
                self.sequence = seq
                break

    def _forwardParsimony(self, aln):
        """ Internal function that operates recursively to first initialise each node (forward),
            stopping only once a sequence has been assigned to the node,
            then to propagate scores from sequence assigned nodes to root (backward). """
        if self.sequence == None:  # no sequence has been assigned
            if self.left == None and self.right == None:  # no children, so terminal, cannot propagate scores
                raise RuntimeError("No sequence assigned to leaf node:", self.label)
            scoresleft = scoresright = None
            if self.left != None:
                scoresleft = self.left._forwardParsimony(aln)
            if self.right != None:
                scoresright = self.right._forwardParsimony(aln)
            # for each position in the alignment,
            # introduce (initially zero) score for each symbol in alphabet
            self.seqscores = [[0 for _ in aln.alphabet] for col in range(aln.alignlen)]
            # for each position in the alignment,
            # allocate a position to put the left child symbol from which each current node symbol score was determined
            self.backleft = [[None for _ in aln.alphabet] for _ in range(aln.alignlen)]
            # allocate a position to put the right child symbol from which each current node symbol score was determined
            self.backright = [[None for _ in aln.alphabet] for _ in range(aln.alignlen)]
            for col in range(aln.alignlen):
                for a_parent in range(len(aln.alphabet)):
                    best_score_left = +9999999
                    best_score_right = +9999999
                    best_symb_left = 0
                    best_symb_right = 0
                    for a_left in range(len(aln.alphabet)):
                        score = (scoresleft[col][a_left] + (
                        1 if a_left != a_parent else 0))  # if we want to weight scores, this would need to change
                        if score < best_score_left:
                            best_symb_left = a_left
                            best_score_left = score
                    for a_right in range(len(aln.alphabet)):
                        score = (scoresright[col][a_right] + (
                        1 if a_right != a_parent else 0))  # if we want to weight scores, this would need to change
                        if score < best_score_right:
                            best_symb_right = a_right
                            best_score_right = score
                    self.seqscores[col][a_parent] = best_score_left + best_score_right
                    self.backleft[col][a_parent] = best_symb_left
                    self.backright[col][a_parent] = best_symb_right
        else:
            self.seqscores = [[0 if a == sym else 999999 for a in aln.alphabet] for sym in
                              self.sequence]  # if we want to weight scores, this would need to change
        return self.seqscores

    def _backwardParsimony(self, aln, seq=None):
        """ Internal function that operates recursively to inspect scores to determine
            most parsimonious sequence, from root to leaves. """
        if self.sequence == None:  # no sequence has been assigned
            leftbuf = []
            rightbuf = []
            if self.left == None and self.right == None:  # no children, so terminal, cannot propagate scores
                raise RuntimeError("No sequence assigned to leaf node:", self.label)
            if seq == None:  # Only root can do this, no parents to consider, so we pick the lowest scoring symbol
                currbuf = []
                for col in range(aln.alignlen):
                    min_score = 999999
                    min_symb = None
                    left_symb = None
                    right_symb = None
                    for a_parent in range(len(aln.alphabet)):
                        if self.seqscores[col][a_parent] < min_score:
                            min_score = self.seqscores[col][a_parent]
                            min_symb = a_parent
                            left_symb = self.backleft[col][a_parent]
                            right_symb = self.backright[col][a_parent]
                    currbuf.append(aln.alphabet[min_symb])
                    leftbuf.append(aln.alphabet[left_symb])
                    rightbuf.append(aln.alphabet[right_symb])
                self.sequence = Sequence(currbuf, aln.alphabet, self.label, gappy=True)
            else:  # Non-root, but not leaf
                self.sequence = seq
                col = 0
                for sym_parent in self.sequence:
                    a_parent = aln.alphabet.index(sym_parent)
                    left_symb = self.backleft[col][a_parent]
                    right_symb = self.backright[col][a_parent]
                    leftbuf.append(aln.alphabet[left_symb])
                    rightbuf.append(aln.alphabet[right_symb])
                    col += 1
            self.left._backwardParsimony(aln, Sequence(leftbuf, aln.alphabet, self.label, gappy=True))
            self.right._backwardParsimony(aln, Sequence(rightbuf, aln.alphabet, self.label, gappy=True))
        return self.sequence

    def getSequence(self):
        """ Get the sequence for the node. Return None if no sequence is assigned.
            Requires that an alignment is associated with the tree, and that sequence names match node labels.
            If the explored node is not a leaf, the sequence can be determined by parsimony. """
        if self.sequence != None:  # a sequence has been assigned
            return self.sequence
        elif self.seqscores != None:  # inferred by parsimony but not yet assigned
            return None  # determine most parsimonous sequence, not yet implemented

    def isAncestorOf(self, node, transitive=True):
        """ Decide if this node is the ancestor of specified node.
            If transitive is True (default), all descendants are included.
            If transitive is False, only direct descendants are included. """
        if node == self.left or node == self.right:
            return True
        elif transitive:
            if self.left:
                statusLeft = self.left.isAncestorOf(node, transitive)
                if statusLeft: return True
            if self.right:
                return self.right.isAncestorOf(node, transitive)
        else:
            return False

    def getDescendants(self, transitive=False):
        """ Retrieve and return (list of) nodes descendant of this.
            If transitive is False (default), only direct descendants are included.
            If transitive is True, all descendants are (recursively) included. """
        children = []
        if self.left:
            children.append(self.left)
        if self.right:
            children.append(self.right)
        if not transitive:
            return children
        else:
            grandchildren = []
            for c in children:
                d = c.getDescendants(transitive)
                if d:
                    grandchildren.extend(d)
            children.extend(grandchildren)
            return children


###############################################################################
# Methods for generating a single tree by clustering, here UPGMA Zvelebil and Baum p. 278
# Methods for processing files of trees on the Newick format
###############################################################################

def runUPGMA(aln, measure, absoluteDistances=False):
    """ Generate an ultra-metric, bifurcating, rooted tree from an alignment based on pairwise distances.
        Use specified distance metric (see sequence.calcDistances).
        If absoluteDistances is True, the tree will be assigned the total distance from provided species.
        Otherwise, the relative addition at each path will be assigned."""
    D = {}
    N = {}  # The number of sequences in each node
    M = aln.calcDistances(measure)  # determine all pairwise distances
    nodes = [PhyloNode(seq.name) for seq in aln.seqs]  # construct all leaf nodes
    """ For each node-pair, assign the distance between them. """
    for i in range(len(nodes)):
        nodes[i].sequence = aln.seqs[i]
        nodes[i].dist = 0.0
        N[nodes[i]] = 1  # each cluster contains a single sequence
        for j in range(0, i):
            D[frozenset([nodes[i], nodes[j]])] = M[i, j]
    """ Now: treat each node as a cluster,
        until there is only one cluster left,
        find the *closest* pair of clusters, and
        merge that pair into a new cluster (to replace the two that merged).
        In each case, the new cluster is represented by the (phylo)node that is formed. """
    while len(N) > 1:  # N will contain all "live" clusters, to be reduced to a signle below
        closest_pair = (None, None)  # The two nodes that are closest to one another according to supplied metric
        closest_dist = None  # The distance between them
        for pair in D:  # check all pairs which should be merged
            dist = D[pair]
            if closest_dist == None or dist < closest_dist:
                closest_dist = dist
                closest_pair = list(pair)
        # So we know the closest, now we need to merge...
        x = closest_pair[0]  # See Zvelebil and Baum p. 278 for notation
        y = closest_pair[1]
        z = PhyloNode()  # create a new node for the cluster z
        z.dist = D.pop(frozenset([x, y])) / 2.0  # assign the absolute distance, travelled so far, note: this will change to relative distance later
        Nx = N.pop(x)  # find number of sequences in x, remove the cluster from list N
        Ny = N.pop(y)  # find number of sequences in y, remove the cluster from list N
        dz = {}  # new distances to cluster z
        for w in N:  # for each node w ...
            # we will merge x and y into a new cluster z, so need to consider w (which is not x or y)
            dxw = D.pop(frozenset([x, w]))  # retrieve and remove distance from D: x to w
            dyw = D.pop(frozenset([y, w]))  # retrieve and remove distance from D: y to w
            dz[w] = (Nx * dxw + Ny * dyw) / (Nx + Ny)  # distance: z to w
        N[z] = Nx + Ny  # total number of sequences in new cluster, insert new cluster in list N
        for w in dz:  # we have to run through the nodes again, now not including the removed x and y
            D[frozenset([z, w])] = dz[w]  # for each "other" cluster, update distance per EQ8.16 (Z&B p. 278)
        z.left = x  # link the phylogenetic tree
        z.right = y
        nodes.append(z)
    if not absoluteDistances:
        x._propagateDistance(z.dist)  # convert absolute distances to relative by recursing down left path
        y._propagateDistance(z.dist)  # convert absolute distances to relative by recursing down right path
        z.dist = 0.0  # root z is at distance 0 from merged x and y
    return PhyloTree(z)  # make it to tree, return

def _findComma(string, level=0):
    """ Find first comma at specified level of embedding """
    mylevel = 0
    for i in range(len(string)):
        if string[i] == '(':
            mylevel += 1
        elif string[i] == ')':
            mylevel -= 1
        elif string[i] == ',' and mylevel == level:
            return i
    return -1


def parseNewickNode(string):
    """ Utility function that recursively parses embedded string using Newick format. """
    first = string.find('(')
    last = string[::-1].find(')')  # look from the back
    if first == -1 and last == -1:  # we are at leaf
        y = string.split(':')
        node = PhyloNode(y[0])
        if len(y) >= 2:
            node.dist = float(y[1])
        return node
    elif first >= 0 and last >= 0:
        # remove parentheses
        last = len(string) - last - 1  # correct index to refer from start instead of end of string
        embed = string[first + 1:last]
        tail = string[last + 1:]
        # find where corresp comma is
        comma = _findComma(embed)
        if comma == -1:
            raise RuntimeError('Invalid format: invalid placement of "," in sub-string "' + embed + '"')
        left = embed[0:comma].strip()
        right = embed[comma + 1:].strip()
        y = tail.split(':')
        node = PhyloNode(y[0])
        if len(y) >= 2:
            node.dist = float(y[1])
        node.left = parseNewickNode(left)
        node.right = parseNewickNode(right)
        return node
    else:
        raise RuntimeError('Invalid format: unbalanced parentheses in sub-string "' + string + '"')


def parseNewick(string):
    """ Main method for parsing a Newick string into a (phylogenetic) tree.
        Handles labels (on both leaves and internal nodes), and includes distances (if provided).
        Returns an instance of a PhyloTree. """
    if string.find(';') != -1:
        string = string[:string.find(';')]
    return PhyloTree(parseNewickNode(string))


def readNewickFile(filename):
    """ Read file on Newick format.
        Returns an instance of a PhyloTree."""
    f = open(filename)
    string = ''.join(f)
    return parseNewick(string)


def writeNewickFile(filename, tree):
    """ Write the specified tree to a Newick file. """
    fh = open(filename, 'w')
    fh.write(tree.__str__())
    fh.close()
