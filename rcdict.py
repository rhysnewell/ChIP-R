import sym

class RCDict(dict):
    """ Class that extends a standard dictionary to accept only fixed-length DNA symbol strings as keys.
        Additionally, it maps the reverse complement to the same value. """

    def __init__(self, alpha = sym.DNA_Alphabet):
        """ Initialise a reverse-complement dictionary to accept strings of a given alphabet (DNA by default) """
        self.alpha = alpha
        self.length = None

    def __setitem__(self, key, value):
        """ Set the value for a key.
            Checks to see that if
            (a) previous keys have been used that the length is the same, and
            (b) the key consists only of valid symbols. """
        if self.length == None:
            self.length = len(key)
        elif len(key) != self.length:
            raise RuntimeError("Invalid key: " + str(key))
        for i in range(len(key)):
            if not key[i] in sym.DNA_Alphabet:
                raise RuntimeError("Invalid symbol in key: " + str(key[i]))
        super(RCDict, self).__setitem__(self.canonical(key), value)

    def canonical(self, key):
        """ Figures out the canonical key (original or its reverse complement).
            Note that is you use other than DNA you may need to modify this code. """
        if self.length == None:
            return key
        alpha = self.alpha
        rcindx = [0 for _ in range(self.length)]
        fwindx = [alpha.index(sym) for sym in key]
        undecided = True
        for forw in range(self.length):
            backw = self.length - forw - 1
            rcindx[forw] = 3 - fwindx[backw] # formula for converting A <-> T, C <-> G
            if undecided and rcindx[forw] > fwindx[forw]: # RC is "higher" than original
                return key
            undecided = rcindx[forw] == fwindx[forw]
        return ''.join([alpha.symbols[indx] for indx in rcindx])

    def __getitem__(self, key):
        """ Retrieve the value associated with a specified key. """
        return super(RCDict, self).__getitem__(self.canonical(key))

    def getSum(self, IUPAC_key):
        """ Retrieve the sum of all the entries that match the specified IUPAC key. """
        # TODO
        pass
