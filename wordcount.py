#!/usr/bin/python

import sys, math, random, getopt
import numpy as np
import matplotlib.pyplot as plt
import prob as prb
import sequence
import stats
from rcdict import *
import operator        # for use with key= in max() function
import binomial

def slidewin(seq, winsize):
    """ Produce a list of sub-sequences of a given length from a complete sequence """
    subseqs = []
    for i in range(len(seq) - winsize + 1):
        subseqs.append(seq[i : i + winsize])
    return subseqs

def countWordsReport(seqs, WordWidth = 8, PeakWidth = 100, PeakMargin = 100):
    """ Produce a report of enriched words of specified length.
        seqs: DNA sequence data
        WordWidth: length of sought words
        PeakWidth: width of window around centre of sequence
        PeakMargin: the width of the margin on each side of the centre window
        (which delineates the positives around peak from negatives away from peak). """
    pos = RCDict() # reverse complement-aware dictionary for DNA
    neg = RCDict() # reverse complement-aware dictionary for DNA
    for seq in seqs:
        centre = len(seq)/2 # find peak
        """ Construct all words around peak (positives) and count their presence """
        words = set(slidewin(seq[centre-PeakWidth/2:centre+PeakWidth/2], WordWidth))
        for word in words:
            try:
                pos[word] += 1
            except KeyError:
                pos[word] = 1
        """ Construct all words away from peak (negatives) and count """
        words = set(slidewin(seq[:centre-PeakWidth/2-PeakMargin], WordWidth))
        words.union(slidewin(seq[centre+PeakWidth/2+PeakMargin:], WordWidth))
        for word in words:
            try:
                neg[word] += 1
            except KeyError:
                neg[word] = 1

    logratio = RCDict() # DNA dictionary for storing the log-ration between pos and neg
    for (word, cnt_pos) in list(pos.items()):
        cnt_neg = 0.0001
        try:
            cnt_neg = neg[word]
        except KeyError:
            pass
        logratio[word] = math.log(float(cnt_pos) / float(cnt_neg))

    allpos = list(logratio.items()) # extract all pairs of words:log-ratio
    sortpos = sorted(allpos, key=lambda v: v[1], reverse=True) # sort them
    print("Enriched words (sorted by ln pos/neg)")
    print("Word    \tln pos/neg\tE-value")
    for (word, lgr) in sortpos[0:100]: # Look at the top-entries according to log-ratio, compute e-values
        cnt_pos = int(pos[word])
        try: cnt_neg = int(neg[word])
        except KeyError: cnt_neg = 0
        # Compute p-value using Fisher's Exact test
        pval = stats.getFETpval(cnt_pos, cnt_neg, len(seqs) * (PeakWidth - WordWidth + 1) - cnt_pos, len(seqs) * (len(seq) - (PeakMargin * 2 + PeakWidth) - (WordWidth - 1) * 2) - cnt_neg, False)
        # Correct for multiple testing (very conservatively)
        eval = pval * len(allpos)
        print("%s\t%6.3f  \t%e" % (word, lgr, eval))

def getReverse(distribs):
    """ Construct a new list of probability distributions of DNA, by
        1. swapping their order, and
        2. swapping A's and T's, and C's and G's """
    return [d.swapxcopy('A','T').swapxcopy('C','G') for d in distribs[::-1]] # backwards


def scanMotifReport(seqs, motif, threshold=0, jaspar = 'JASPAR_matrices.txt'):
    """ Produce a plot for a scan of the specified motif.
        The plot has as its x-axis position of sequence, and
        the y-axis the cumulative, non-negative PWM score over all sequences. """
    # check that all sequences are the same length and set sequence length
    seq_len = len(seqs[0])
    for seq in seqs:
        if len(seq) != seq_len:
            usage(sys.argv[0], "All sequences must have same length")
            return

    # create the motif and its reverse complemennt
    bg = prb.Distrib(sym.DNA_Alphabet, sequence.getCount(seqs))
    d = prb.readMultiCounts(jaspar)
    try:
        fg1 = d[motif]
        fg2 = getReverse(d[motif])
    except KeyError:
        usage(sys.argv[0], "Unknown motif %s" % motif)
        return
    print("Motif %s:" % motif)
    pwm1 = sequence.PWM(fg1, bg)
    pwm1.display(format='JASPAR')
    print("Motif %s (reverse complement):" % motif)
    pwm2 = sequence.PWM(fg2, bg)
    pwm2.display(format='JASPAR')

    # initialize things to zero
    avg_motif_score = np.zeros(seq_len)

    # compute average score at each position (on both strands) in sequences
    i_seq = 0
    motif_width = pwm1.length
    for seq in seqs:
        i_seq += 1
        # print >> sys.stderr, "Scoring seq: %4d\r" % (i_seq),

        # positive strand
        hits = pwm1.search(seq, threshold)
        pos_scores = seq_len * [0]
        for hit in hits:
            # mark hit at *center* of site (hence motif_width/2)
            pos_scores[hit[0]+(motif_width/2)] = hit[2]

        # negative strand
        hits = pwm2.search(seq, threshold)
        neg_scores = seq_len * [0]
        for hit in hits:
            neg_scores[hit[0]+(motif_width/2)] = hit[2]

        # use maximum score on two strands
        for i in range(seq_len):
            score = max(pos_scores[i], neg_scores[i])
            if (score > threshold):
                avg_motif_score[i] += score

    # compute average score
    for i in range(seq_len):
        avg_motif_score[i] /= len(seqs)

    # hw = 5 # window width is 2*hw + 1
    # smoothed_avg_motif_score = np.zeros(seq_len)
    # for i in range(hw, seq_len-motif_width+1-hw):
    #    smoothed_avg_motif_score[i]=sum(avg_motif_score[i-hw:i+hw+1])/(2*hw+1)

    # plot the average score curve
    # print >> sys.stderr, ""
    x = list(range(-(seq_len/2), (seq_len/2)))    # call center of sequence X=0
    lbl = "%s" % (motif)
    plt.plot(x, avg_motif_score, label=lbl)
    #plt.plot(x, smoothed_avg_motif_score, label=lbl)
    plt.axhline(color='black', linestyle='dotted')
    plt.legend(loc='lower center')
    plt.xlabel('position')
    plt.ylabel('average motif score')
    plt.title(motif)
    plt.show()


def scanMotifReport_new(seqs, motif, threshold=3.4567, jaspar = 'JASPAR_matrices.txt', seed=0):
    """ Produce a plot for a scan of the specified motif.
        The plot has as its x-axis position of sequence, and
        the y-axis the number of sequences with a best hit at position x.
        Sequences with no hit above 'threshold' are ignored.
        Ties for best hit are broken randomly.
        The p-value of the central region that is most "centrally enriched"
        and the width of the best central region is printed in the label
        of the plot.
    """

    # set the random seed for repeatability
    random.seed(seed)

    # Copy the code from your "improved" version of scanMotifReport()
    # to here, and follow the instructions in the Prac to develop this
    # new function.

    # check that all sequences are the same length and set sequence length
    seq_len = len(seqs[0])
    for seq in seqs:
        if len(seq) != seq_len:
            usage(sys.argv[0], "All sequences must have same length")
            return

    # create the motif and its reverse complemennt
    bg = prb.Distrib(sym.DNA_Alphabet, sequence.getCount(seqs))
    d = prb.readMultiCounts(jaspar)
    try:
        fg1 = d[motif]
        fg2 = getReverse(d[motif])
    except KeyError:
        usage(sys.argv[0], "Unknown motif %s" % motif)
        return
    print("Motif %s:" % motif)
    pwm1 = sequence.PWM(fg1, bg)
    pwm1.display(format='JASPAR')
    print("Motif %s (reverse complement):" % motif)
    pwm2 = sequence.PWM(fg2, bg)
    pwm2.display(format='JASPAR')

    # initialize things to zero
    hit_count = np.zeros(seq_len)
    n_seqs_with_hits = 0.0

    # Scan each sequence for all hits on both strands and record
    # the number of "best hits" at each sequence position.
    #
    motif_width = pwm1.length
    i_seq = 0
    for seq in seqs:
        i_seq += 1
        # print >> sys.stderr, "Scoring seq: %4d\r" % (i_seq),
        # scan with both motifs
        hits = pwm1.search(seq, threshold) + pwm2.search(seq, threshold)
        # Record position of best hit
        if (hits):
                n_seqs_with_hits += 1
                # find best hit score
                best_score = max(hits, key=operator.itemgetter(1))[2]
                # find ties
                best_hits = [ hit for hit in hits if hit[2] == best_score ]
                # break ties at random
                best_hit = random.choice(best_hits)
                # mark hit at *center* of site (hence pwm1.length/2)
                hit_count[best_hit[0] + pwm1.length/2] += 1
    # divide number of sequences with hit by total number of hits
    site_probability = [ (cnt/n_seqs_with_hits) for cnt in hit_count ]

    print("Number of sequences with hit (score >= %f): %d" % (threshold, n_seqs_with_hits), file=sys.stderr)

    # STATISTICS
    # Get the cumulative hit counts in concentric windows
    # and perform the Binomial Test.  Report best region and its p-value.
    #
    best_r = 0
    best_log_pvalue = 1
    center = seq_len/2                  # center of sequence
    cum_hit_count = np.zeros(seq_len)   # total hits in window of width i
    for i in range(1, (seq_len - pwm1.length/2 + 1)/2):
        cum_hit_count[i] = cum_hit_count[i-1] + hit_count[center-i] + hit_count[center+i]
        # Compute probability of observed or more best hits in central window
        # assuming uniform probability distribution in each sequence.
    #   successes = cum_hit_count[i]
    #   trials = n_seqs_with_hits
    #    p_success = ?
    #    log_pvalue = ?
    #    if (log_pvalue < best_log_pvalue):
    #        best_log_pvalue = log_pvalue
    #        best_r = 2*i
    # End STATISTICS

    hw = 5
    smoothed_site_probability = np.zeros(seq_len)
    for i in range(hw, seq_len-motif_width+1-hw):
        smoothed_site_probability[i]=sum(site_probability[i-hw:i+hw+1])/(2*hw+1)

    x = list(range(-(seq_len/2), (seq_len/2)))        # call center of sequence X=0
    lbl = "%s, t=%.2f" % (motif, threshold)
    #lbl = "%s, t=%.2f, w=%d, p=%.2e" % (motif, threshold, best_r, math.exp(best_log_pvalue))
    plt.plot(x, smoothed_site_probability, label=lbl)
    plt.axhline(color='black', linestyle='dotted')
    plt.legend(loc='lower center')
    plt.xlabel('Position of best site')
    plt.ylabel('Smoothed probability')
    plt.title(motif)
    plt.show()

def usage(name, errmsg = None):
    if errmsg != None:
        print("Error: %s" % errmsg)
    print("""Usage: %s [options]
                -f <fasta-filename> (required)
                -d discover enriched words
                -w <word width, default 8>
                -p <peak width, default 100>
                -m <peak margin, default 100>
                -s <JASPAR-ID> scan for JASPAR motif
                -h print this help""" % name)

if __name__ == '__main__':
    try:
        optlst, args = getopt.getopt(sys.argv[1:], 'f:hds:j:w:p:m:')
    except getopt.GetoptError as err:
        usage(sys.argv[0], str(err))
        sys.exit(2)
    FILENAME =      None
    DISCOVER_MODE = False
    SCAN_MODE =     False
    WORD_WIDTH =    8
    PEAK_WIDTH =    100
    PEAK_MARGIN =   100
    MOTIF_ID =      'MA0112.2'
    JASPAR_FILE =   'JASPAR_matrices.txt'
    for o, a in optlst:
        if   o == '-h': usage(sys.argv[0])
        elif o == '-f': FILENAME = a
        elif o == '-d': DISCOVER_MODE = True
        elif o == '-w': WORD_WIDTH = int(a)
        elif o == '-p': PEAK_WIDTH = int(a)
        elif o == '-m': PEAK_MARGIN = int(a)
        elif o == '-s': SCAN_MODE = True; MOTIF_ID = a
        elif o == '-j': JASPAR_FILE = a
    if FILENAME == None:
        usage(sys.argv[0], "Filename not specified")
        sys.exit(3)
    seqs = sequence.readFastaFile(FILENAME, sym.DNA_Alphabet_wN)
    if DISCOVER_MODE:
        print("Discover (f=%s; w=%d; p=%d; m=%d)" % (FILENAME, WORD_WIDTH, PEAK_WIDTH, PEAK_MARGIN))
        countWordsReport(seqs, WORD_WIDTH, PEAK_WIDTH, PEAK_MARGIN)
    elif SCAN_MODE:
        scanMotifReport(seqs, MOTIF_ID)
    else:
        usage(sys.argv[0], "No run mode selected")

