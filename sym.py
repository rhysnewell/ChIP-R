"""
Module symbol is for defining alphabets (of symbols), and
for storing and operating on symbols and tuples (ordered or
unordered).
"""
import os

# ------------------ Alphabet ------------------

class Alphabet(object):
    """ Defines an immutable biological alphabet (e.g. the alphabet for DNA is AGCT)
    that can be used to create sequences (see sequence.py).
    We use alphabets to define "tuple" tables, where entries are keyed by combinations
    of symbols of an alphabet (see class TupleStore below).
    Alphabets are used to define probability distributions for stochastic events
    (see prob.py). """

    def __init__(self, symbolString):
        """ Construct an alphabet from a string of symbols. Lower case characters
        will be converted to upper case, repeated characters are ignored.
        Example of constructing the DNA alphabet:
        >>> alpha = Alphabet('ACGTttga')
        >>> alpha.symbols
        ('A', 'C', 'G', 'T') """

        # Add each symbol to the symbols list, one at a time, and ignore doubles (could use "set" here...)
        _symbols = [] # create a temporary list
        for s in symbolString:
            if not str(s).upper()[0] in _symbols:
                _symbols.append(str(s).upper()[0])
        _symbols.sort() # we put them in alphabetical (one canonical) order
        # OK done extracting, put them in place
        self.symbols = tuple(_symbols); # create the immutable tuple from the extracted list
        self.length = len(self.symbols)
        self.annotations = {}

    def __str__(self):
        return str(self.symbols)

    def __len__(self):
        return len(self.symbols)

    def __iter__(self):
        return self.symbols.__iter__()

    def __getitem__(self, ndx):
        """ Retrieve the symbol(s) at the specified index (or slice of indices) """
        return self.symbols[ndx]

    def __contains__(self, sym):
        """ Check if the given symbol is a member of the alphabet. """
        return sym in self.symbols

    def index(self, sym):
        """ Retrieve the index of the given symbol in the alphabet. """
        # If the symbol is valid, use the tuple's index function
        if sym in self.symbols:
            syms = self.symbols
            return syms.index(sym)
        else:
            raise RuntimeError('Symbol %s is not indexed by alphabet %s' % (sym, str(self.symbols)))

    def __eq__(self, rhs):
        """ Test if the rhs alphabet is equal to ours. """
        if rhs == None:
            return False
        if len(rhs) != len(self):
            return False
        # OK we know they're same size...
        for sym in self.symbols:
            if not sym in rhs:
                return False
        return True

    def isSubsetOf(self, alpha2):
        """ Test if this alphabet is a subset of alpha2. """
        for sym in self.symbols:
            if not alpha2.isValidSymbol(sym):
                return False
        return True

    def isSupersetOf(self, alpha2):
        """ Test if this alphabet is a superset of alpha2. """
        return alpha2.isSubsetOf(self)

    def annotateSym(self, label, sym, value):
        try:
            lookup = self.annotations[label]
        except KeyError:
            lookup = self.annotations[label] = {}
        lookup[sym] = value

    def annotateAll(self, label, symdictOrFilename):
        if isinstance(symdictOrFilename, str): # we assume it is a filename
            fh = open(symdictOrFilename)
            string = fh.read()
            d = {}
            for line in string.splitlines():
                if len(line.strip()) == 0:
                    continue
                sections = line.split()
                symstr, value = sections[0:2]
                for sym in symstr:
                    d[sym] = value
            fh.close()
        else: # we assume it is a dictionary
            d = symdictOrFilename
        for sym in d:
            self.annotateSym(label, sym, d[sym])

    def getAnnotation(self, label, sym):
        try:
            lookup = self.annotations[label]
            return lookup[sym]
        except KeyError:
            return None


""" Below we declare alphabets that are going to be available when
this module is imported """
Bool_Alphabet = Alphabet('TF')
DNA_Alphabet = Alphabet('ACGT')
DNA_Alphabet_wN = Alphabet('ACGTN')
RNA_Alphabet_wN = Alphabet('ACGUN')
RNA_Alphabet = Alphabet('ACGU')
Protein_Alphabet = Alphabet('ACDEFGHIKLMNPQRSTVWY')
Protein_Alphabet_wX = Protein_wX = Alphabet('ACDEFGHIKLMNPQRSTVWYX')
Protein_Alphabet_wSTOP = Protein_wSTOP = Alphabet('ACDEFGHIKLMNPQRSTVWY*')
DSSP_Alphabet = Alphabet('GHITEBSC')
DSSP3_Alphabet = Alphabet('HEC')

predefAlphabets = {'Bool_Alphabet': Bool_Alphabet,
                   'DNA': DNA_Alphabet,
                   'RNA': RNA_Alphabet,
                   'DNAwN': RNA_Alphabet_wN,
                   'RNAwN': DNA_Alphabet_wN,
                   'Protein': Protein_Alphabet,
                   'ProteinwX': Protein_wX,
                   'ProteinwSTOP' : Protein_wSTOP,
                   'DSSP_Alphabet' : DSSP_Alphabet,
                   'DSSP3_Alphabet' : DSSP3_Alphabet}
# The preferred order in which a predefined alphabet is assigned to a sequence
# (e.g., we'd want to assign DNA to 'AGCT', even though Protein is also valid)
preferredOrder = ['Bool_Alphabet', 'DNA', 'RNA', 'DNAwN', 'RNAwN', 'Protein', 'ProteinwX', 'ProteinwSTOP', 'DSSP_Alphabet', 'DSSP3_Alphabet']
# Useful annotations
DNA_Alphabet.annotateAll('html-color', {'A':'green','C':'orange','G':'red','T':'#66bbff'})
RNA_Alphabet.annotateAll('html-color', {'A':'green','C':'orange','G':'red','U':'#66bbff'})
Protein_Alphabet.annotateAll('html-color', {'G':'orange','P':'orange','S':'orange','T':'orange','H':'red','K':'red','R':'red','F':'#66bbff','Y':'#66bbff','W':'#66bbff','I':'green','L':'green','M':'green','V':'green'})

# ------------------ Substitution Matrix ------------------

class TupleStore(dict):
    """ Internal utility class that can be used for associating
    a value with ordered n-tuples (n=1..N).
    Read/write functions are defined for instances of this class.
    """

    def __init__(self, alphas=None, entries=None, sparse=True):
        """
        Manage entries keyed by symbol-tuples with values of arbitrary type.
        If alphas is None, the alphabet(s) are inferred from the provided entries.
        If entries is None, all entries are defined by possible combinations of symbols from specified alphabets,
        and are assumed to be None until specified. Either alphas or entries must be supplied.
        If sparse is True, a sparse memory-saving encoding is used, if false, a time-saving, more flexible encoding is used.
        >>> matrix = TupleStore({'AA': 2, 'AW': -3, 'WW': 4, 'AR': -1})
        >>> matrix[('A', 'W')]
        -3
        >>> matrix['AR']
        -1
        """
        assert sparse, "Currently only sparse encoding is implemented."
        assert alphas or entries, "Either alphabets or entries (from which alphabets can be inferred) must be supplied."
        self.sparse = sparse         # sparse encoding if true
        if alphas == None:
            self.alphas = None       # need to figure out alphabet from supplied entries
            self.keylen = None       # tuple length not known yet
        elif type(alphas) is Alphabet:
            self.alphas = tuple ([ alphas ]) # make it into a tuple
            self.keylen = 1          # tuple length 1
        else:
            self.alphas = alphas     # alphabets are supplied
            self.keylen = len(alphas)# length of tuples is the same as the number alphabets

        # Check if entries are supplied to the constructor
        if entries == None:
            self.entries = entries = {}
        elif type(entries) is dict:
            raise RuntimeError("When specified, entries must be a dictionary")
        # Check length of tuples, must be the same for all
        for entry in entries:
            if self.keylen == None:
                self.keylen = len(entry)
            elif self.keylen != len(entry):
                raise RuntimeError("All entries must have the same number of symbols")

        # go through each position in tuples, to check what alphabet is right
        myalphas = []                   # my suggestions from entries (need to be subsets of specified)
        for idx in range(self.keylen):
            symset = set()              # we collect all symbols in position idx here
            for key in entries:
                symset.add(key[idx])
            myalpha = Alphabet(symset)
            myalphas.append(myalpha)
            if self.alphas != None:     # if specified it needs to be a superset of that we constructed
                if not self.alphas[idx].isSupersetOf(myalpha):
                    raise RuntimeError("Specified alphabet is not compatible with specified entries")

        if self.alphas == None:     # if not specified to constructor use those we found
            self.alphas = tuple(myalphas)

        for key in entries:
            self[key] = entries[key]

    def _isValid(self, symkey):
        for idx in range(self.keylen):
            if not symkey[idx] in self.alphas[idx]:
                return False
        return True

    def __setitem__(self, symkey, value):
        assert self.keylen == len(symkey), "All entries in dictionary must be equally long"
        assert self._isValid(symkey), "Invalid symbol in entry"
        self.entries[symkey] = value

    def __getitem__(self, symkey):
        """ Return the score matching the given symbols together."""
        assert self.keylen == len(symkey), "Entries must be of the same length"
        try:
            return self.entries[symkey]
        except KeyError:
            return None

    def __iadd__(self, symkey, ivalue):
        assert self.keylen == len(symkey), "All entries in dictionary must be equally long"
        assert self._isValid(symkey), "Invalid symbol in entry"
        try:
            self.entries[symkey] += ivalue
        except KeyError:
            self.entries[symkey] = ivalue

    def __isub__(self, symkey, ivalue):
        assert self.keylen == len(symkey), "All entries in dictionary must be equally long"
        assert self._isValid(symkey), "Invalid symbol in entry"
        try:
            self.entries[symkey] -= ivalue
        except KeyError:
            self.entries[symkey] = -ivalue

    def getAll(self, symkey=None):
        """ Return the values matching the given symbols together.
        symkey: tuple (or list) of symbols or None (symcount symbol); if tuple is None, all entries are iterated over.
        """
        if symkey == None:
            symkey = []
            for idx in range(self.keylen):
                symkey.append(None)
        else:
            assert self.keylen == len(symkey), "Entries must be of the same length"
        for idx in range(self.keylen):
            if symkey[idx] != None:
                if not symkey[idx] in self.alphas[idx]:
                    raise RuntimeError("Invalid entry: must be symbols from specified alphabet or None")
        return TupleEntries(self, symkey)

    def __iter__(self):
        return TupleEntries(self, tuple([None for _ in range(self.keylen)]))

    def items(self, sort = False):
        """ In a dictionary-like way return all entries as a list of 2-tuples (key, prob).
        If sort is True, entries are sorted in descending order of value.
        Note that this function should NOT be used for big (>5 variables) tables."""
        ret = []
        for s in self.entries:
            if self[s] != None:
                ret.append((s, self[s]))
        if sort:
            return sorted(ret, key=lambda v: v[1], reverse=True)
        return ret

class TupleEntries(object):
    """ Iterator class for multiple entries in a tuple store.
    """
    def __init__(self, tuplestore, symkey):
        self.tuplestore = tuplestore
        self.symkey = symkey
        self.symcount = []
        self.indices = []
        for ndx in range(tuplestore.keylen):
            if symkey[ndx] == None:
                self.indices.append(ndx)
                self.symcount.append(0)        # start at this index to alter symbol
            else:
                self.symcount.append(None)     # do not alter this symbol
        self.nextIsLast = False

    def __iter__(self):
        return self

    def __next__(self):
        """ Step through sequence of entries, either
        (if not sparse) with a step-size based on alphabet-sizes and what symbols are specified or
        (if sparse) with calls to tuple store based on all possible symbol combinations."""

        if self.nextIsLast:
            raise StopIteration

        mykey = [] # construct current combination from known and unspecified symbols
        for ndx in range(self.tuplestore.keylen):
            if (self.symkey[ndx] == None):
                sym = self.tuplestore.alphas[ndx][self.symcount[ndx]]
                mykey.append(sym)
            else:
                mykey.append(self.symkey[ndx])

        # decide which ndx that should be increased (only one)
        self.nextIsLast = True # assume this is the last round (all counters are re-set)
        for ndx in self.indices:
            if self.symcount[ndx] == len(self.tuplestore.alphas[ndx]) - 1: # if we just entered the last symbol of this alphabet
                self.symcount[ndx] = 0                  # reset count here
            else:
                self.symcount[ndx] = self.symcount[ndx] + 1
                self.nextIsLast = False
                break

        return tuple(mykey)

