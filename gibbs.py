"""
Motif discovery using Gibb's sampling
@author: mikael
"""

import math
import random

import sym
import prob
import sequence

class GibbsMotif():
    """
    A class for discovering linear motifs in sequence data.
    Uses Gibb's sampling (Lawrence et al., Science 262:208-214 1993).

    Also see http://bayesweb.wadsworth.org/gibbs/content.html which has info
    on "site sampling", "motif sampling", "recursive sampling" and "centroid
    sampling". The first is implemented (roughly) below.
    """
    def __init__(self, seqs, length, alignment = None):
        """ Construct a "discovery" session by providing the sequences that will be used.
            seqs: sequences in which the motif is sought
            length: length of sought pattern (W)
            alignment: positions in each sequence for the initial alignment (use only if the alignment
            has been determined from a previous run).
            """
        self.seqs = seqs
        self.length = length # length of motif 1..W
        seqs = self.seqs
        self.alphabet = None
        k = 0
        for s in seqs:
            if self.alphabet != None and self.alphabet != s.alphabet:
                raise RuntimeError("Sequences invalid: different alphabets")
            self.alphabet = s.alphabet
            if alignment:
                if alignment[k] < 0 or alignment[k] >= len(s):
                    raise RuntimeError("Initial alignment invalid: does not match sequence " + s.name)
            k += 1
        """ Initialise parameters that are part of the setup (below) """
        self.alignment = alignment or [ random.randint(0, len(s) - length) for s in seqs ] # starting positions defining alignment

    def discover(self, pseudocount = None, niter = None):
        """ Find the most probable common pattern represented by a
            position weight matrix (PWM), based on W+1 distributions
            pseudocount: the distribution used for pseudo-counts (default is uniform)
            niter: number of iterations (if None, 100*N is used; where N is number of seqs).
        """
        """ Initialise parameters necessary for the discovery run (below) """
        N = len(self.seqs) # number of sequences 1..N
        seqs = self.seqs
        W = self.length    # motif width
        """ background that will be used as pseudo-counts """
        pseudocount = pseudocount or prob.Distrib(self.alphabet, 1.0)
        """ q: the foreground distribution (specifying the W distributions in aligned columns)
            p: the background distribution (for non-aligned positions in all sequences) """
        q = [ prob.Distrib(self.alphabet, pseudocount) for _ in range(W) ]
        p = prob.Distrib(self.alphabet, pseudocount)
        a = self.alignment

        new_z = random.randint(0, N-1) # pick a random sequence to withhold
        for k in range(N):
            if k != new_z:
                k_len = len(seqs[k]) # length of current seq
                offset = 0
                for i in range(k_len):
                    if i >= a[k] and i < a[k] + W: # within pattern
                        q[offset].observe(seqs[k][i])
                        offset += 1
                    else: # outside pattern
                        p.observe(seqs[k][i])

        """ Main loop: predictive update step THEN sampling step, repeat... """
        niter = niter or 100 * N # use specified number of iterations or default
        for round in range(niter):

            """ Predictive update step:
                One of the N sequences are chosen at random: z.
                We will not use it in the profile, nor background so we
                exclude it from our counts. """
            prev_z = new_z
            new_z = random.randint(0, N - 1)
            # q's and p's are updated from current a's and all sequences except z,
            # which is the same as use old q's and p's and subtract z's contribs...
            offset = 0
            for i in range(len(seqs[new_z])):
                if i >= a[new_z] and i < a[new_z] + W: # within pattern
                    q[offset].observe(seqs[new_z][i], -1) # subtract the count
                    offset += 1
                else: # outside pattern
                    p.observe(seqs[new_z][i], -1) # subtract the count
            # ... and add back the previous and now updated z
            offset = 0
            for i in range(len(seqs[prev_z])):
                if i >= a[prev_z] and i < a[prev_z] + W: # within pattern
                    q[offset].observe(seqs[prev_z][i], +1) # add the count
                    offset += 1
                else: # outside pattern
                    p.observe(seqs[prev_z][i], +1) # add the count

            """ Sampling step:
                Consider each position x in z as a match: find a weight Ax """
            z_len = len(seqs[new_z]) # length of seq z
            A = [ 0.0 for _ in range(z_len) ]
            Asum = 0.0
            for x in range(z_len - W + 1): # look at all starts for a W-wide pattern
                Px = 1.0; Qx = 1.0
                for w in range(W):
                    Px *= p[seqs[new_z][x+w]]
                    Qx *= q[w][seqs[new_z][x+w]]
                try:
                    A[x] = Qx / Px
                except ZeroDivisionError:
                    pass
                Asum += A[x]
            for x in range(z_len - W + 1): # score all starts for a W-wide pattern
                A[x] /= Asum               # normalise so that all Ax's sum to 1.0
            # Pick the next a[z], with a probability proportional to Ax
            pick = random.random()         # any value between 0 and 1
            cumul = 0.0                    # cumulative probability
            for x in range(z_len - W + 1): # check starts for a W-wide pattern
                cumul += A[x]
                if pick <= cumul:          # check if our random pick is smaller than the cumulative prob
                    a[new_z] = x
                    break

            """ Evaluate data log-likelihood """
            if round % 100 == 0: # but only every 100th round
                LL = 0.0
                for k in range(N):
                    Pk = 1.0; Qk = 1.0
                    for w in range(W):
                        Pk *= p[seqs[k][a[k]+w]]
                        Qk *= q[w][seqs[k][a[k]+w]]
                    try:
                        LL += math.log(Qk / Pk)
                    except ZeroDivisionError:
                        pass
                print("LL @ %5d=\t%5.2f" % (round, LL))

        # end main for-loop
        self.q = q
        self.p = p
        self.alignment = a
        return q

    def getForeground(self):
        """ Return the probability distributions for columns in the discovered alignment. """
        return self.q

    def getBackground(self):
        """ Return the probability distributions for the background used in the discovery. """
        return self.p

def getAlignment(seqs, motif, background):
    """ Retrieve the best alignment (positions) in provided sequences defined by the specified
        motif params.
        seqs: sequence data
        motif: the foreground distribution (specifying the W distributions in aligned columns)
        background: the background distribution (for non-aligned positions in all sequences)
        Note that this is similar but not the same as the stochastically selected alignment that
        is kept while training. It can be implemented using a PWM constructed from a previous session.
        Note also that this alignment can be used as input to continue an earlier discovery session
        when motif distributions had been saved.  """
    N = len(seqs)
    q = motif
    p = background
    W = len(q)
    a = [0 for _ in range(N)] # start positions unknown
    for k in range(N):
        k_len = len(seqs[k]) # length of seq k
        Amax = None
        xmax = 0
        for x in range(k_len - W + 1):
            Px = 1.0; Qx = 1.0
            for w in range(W):
                Px *= p[seqs[k][x+w]]
                Qx *= q[w][seqs[k][x+w]]
            try:
                Atmp = math.log(Qx / Px)
            except ZeroDivisionError:
                pass
            if Amax == None or Amax < Atmp:
                Amax = Atmp
                xmax = x
        a[k] = xmax
    return a

class GibbsAlign():
    """ A class for performing ungapped sequence alignment.
        Uses Gibb's sampling (Lawrence et al., Science 262:208-214 1993).
    """

    def __init__(self, seqs, length, alignment = None):
        """ Construct a "discover" session by providing the sequences that will be aligned.
            seqs: sequences that will be aligned
            length: maximum length of alignment (must be equal or greater than max sequence length)
            alignment: positions in each sequence for the initial alignment (use only if the alignment
            has been determined from a previous run).
            """
        self.seqs = seqs
        self.length = length # length of motif 1..W
        seqs = self.seqs
        self.alphabet = None
        k = 0
        for s in seqs:
            if self.alphabet != None and self.alphabet != s.alphabet:
                raise RuntimeError("Sequences invalid: different alphabets")
            self.alphabet = s.alphabet
            if alignment:
                if alignment[k] < 0 or alignment[k] >= len(s):
                    raise RuntimeError("Initial alignment invalid: does not match sequence " + s.name)
            k += 1
        """ Initialise parameters that are part of the setup (below) """
        self.alignment = alignment or [ random.randint(0, length - len(s)) for s in seqs ] # starting offsets defining alignment

    def discover(self, pseudocount = None, niter = None):
        """ Find the most probable common pattern represented by a
            position weight matrix (PWM), based on W+1 distributions
            pseudocount: the distribution used for pseudo-counts (default is uniform)
            niter: number of iterations (if None, 100*N is used; where N is number of seqs).
        """
        """ Initialise parameters necessary for the discovery run (below) """
        N = len(self.seqs) # number of sequences 1..N
        seqs = self.seqs
        W = self.length    # alignment width
        """ background that will be used as pseudo-counts """
        pseudocount = pseudocount or prob.Distrib(self.alphabet, 1.0)
        """ q: the foreground distribution (specifying the W distributions in aligned columns)
            p: the background distribution (for non-aligned positions in all sequences) """
        q = [ prob.Distrib(self.alphabet, pseudocount) for _ in range(W) ]
        p = prob.Distrib(self.alphabet, pseudocount)
        a = self.alignment

        new_z = random.randint(0, N-1) # pick a random sequence to withhold
        for k in range(N):
            if k != new_z:
                k_len = len(seqs[k]) # length of current seq
                offset = 0
                for i in range(k_len):
                    if i >= a[k] and i < a[k] + W: # within pattern
                        q[offset].observe(seqs[k][i])
                        offset += 1
                    else: # outside pattern
                        p.observe(seqs[k][i])

        """ Main loop: predictive update step THEN sampling step, repeat... """
        niter = niter or 100 * N # use specified number of iterations or default
        for round in range(niter):

            """ Predictive update step:
                One of the N sequences are chosen at random: z.
                We will not use it in the profile, nor background so we
                exclude it from our counts. """
            prev_z = new_z
            new_z = random.randint(0, N - 1)
            # q's and p's are updated from current a's and all sequences except z,
            # which is the same as use old q's and p's and subtract z's contribs...
            offset = 0
            for i in range(len(seqs[new_z])):
                if i >= a[new_z] and i < a[new_z] + W: # within pattern
                    q[offset].observe(seqs[new_z][i], -1) # subtract the count
                    offset += 1
                else: # outside pattern
                    p.observe(seqs[new_z][i], -1) # subtract the count
            # ... and add back the previous and now updated z
            offset = 0
            for i in range(len(seqs[prev_z])):
                if i >= a[prev_z] and i < a[prev_z] + W: # within pattern
                    q[offset].observe(seqs[prev_z][i], +1) # add the count
                    offset += 1
                else: # outside pattern
                    p.observe(seqs[prev_z][i], +1) # add the count

            """ Sampling step:
                Consider each position x in z as a match: find a weight Ax """
            z_len = len(seqs[new_z]) # length of seq z
            A = [ 0.0 for _ in range(z_len) ]
            Asum = 0.0
            for x in range(z_len - W + 1): # look at all starts for a W-wide pattern
                Px = 1.0; Qx = 1.0
                for w in range(W):
                    Px *= p[seqs[new_z][x+w]]
                    Qx *= q[w][seqs[new_z][x+w]]
                try:
                    A[x] = Qx / Px
                except ZeroDivisionError:
                    pass
                Asum += A[x]
            for x in range(z_len - W + 1): # score all starts for a W-wide pattern
                A[x] /= Asum               # normalise so that all Ax's sum to 1.0
            # Pick the next a[z], with a probability proportional to Ax
            pick = random.random()         # any value between 0 and 1
            cumul = 0.0                    # cumulative probability
            for x in range(z_len - W + 1): # check starts for a W-wide pattern
                cumul += A[x]
                if pick <= cumul:          # check if our random pick is smaller than the cumulative prob
                    a[new_z] = x
                    break

            """ Evaluate data log-likelihood """
            if round % 100 == 0: # but only every 100th round
                LL = 0.0
                for k in range(N):
                    Pk = 1.0; Qk = 1.0
                    for w in range(W):
                        Pk *= p[seqs[k][a[k]+w]]
                        Qk *= q[w][seqs[k][a[k]+w]]
                    try:
                        LL += math.log(Qk / Pk)
                    except ZeroDivisionError:
                        pass
                print("LL @ %5d=\t%5.2f" % (round, LL))

        # end main for-loop
        self.q = q
        self.p = p
        self.alignment = a
        return q

    def getForeground(self):
        """ Return the probability distributions for columns in the discovered alignment. """
        return self.q

    def getBackground(self):
        """ Return the probability distributions for the background used in the discovery. """
        return self.p
