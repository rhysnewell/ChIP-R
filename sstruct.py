'''
Module sstruct -- methods for protein secondary structure
'''

import sequence
import sym

cf_dict = {  # Chou-Fasman table
#     P(a), P(b), P(t),    f(i), f(i+1), f(i+2), f(i+3)
'A': ( 142,   83,   66,   0.060,  0.076,  0.035,  0.058 ),    # Alanine
'R': (  98,   93,   95,   0.070,  0.106,  0.099,  0.085 ),    # Arginine
'N': ( 101,   54,  146,   0.147,  0.110,  0.179,  0.081 ),    # Aspartic Acid
'D': (  67,   89,  156,   0.161,  0.083,  0.191,  0.091 ),    # Asparagine
'C': (  70,  119,  119,   0.149,  0.050,  0.117,  0.128 ),    # Cysteine
'E': ( 151,   37,   74,   0.056,  0.060,  0.077,  0.064 ),    # Glutamic Acid
'Q': ( 111,  110,   98,   0.074,  0.098,  0.037,  0.098 ),    # Glutamine
'G': (  57,   75,  156,   0.102,  0.085,  0.190,  0.152 ),    # Glycine
'H': ( 100,   87,   95,   0.140,  0.047,  0.093,  0.054 ),    # Histidine
'I': ( 108,  160,   47,   0.043,  0.034,  0.013,  0.056 ),    # Isoleucine
'L': ( 121,  130,   59,   0.061,  0.025,  0.036,  0.070 ),    # Leucine
'K': ( 114,   74,  101,   0.055,  0.115,  0.072,  0.095 ),    # Lysine
'M': ( 145,  105,   60,   0.068,  0.082,  0.014,  0.055 ),    # Methionine
'F': ( 113,  138,   60,   0.059,  0.041,  0.065,  0.065 ),    # Phenylalanine
'P': (  57,   55,  152,   0.102,  0.301,  0.034,  0.068 ),    # Proline
'S': (  77,   75,  143,   0.120,  0.139,  0.125,  0.106 ),    # Serine
'T': (  83,  119,   96,   0.086,  0.108,  0.065,  0.079 ),    # Threonine
'W': ( 108,  137,   96,   0.077,  0.013,  0.064,  0.167 ),    # Tryptophan
'Y': (  69,  147,  114,   0.082,  0.065,  0.114,  0.125 ),    # Tyrosine
'V': ( 106,  170,   50,   0.062,  0.048,  0.028,  0.053 ),    # Valine
'Y': (  69,  147,  114,   0.082,  0.065,  0.114,  0.125 ),    # Tyrosine
'V': ( 106,  170,   50,   0.062,  0.048,  0.028,  0.053 ),}   # Valine

prot_alpha = sym.Protein_Alphabet
sstr_alpha = sym.DSSP3_Alphabet

def makesstr(seq, sym = '*', gap = '-'):
    """ Create a string from a list of booleans (seq) that indicate with sym what elements are true.
        gap is used for elements that are false.
    """
    sstr = ''
    for yes in seq:
        if yes:
            sstr += sym
        else:
            sstr += gap
    return sstr

def markCountAbove(scores, width = 6, call_cnt = 4):
    """ Create a list of booleans that mark all positions within a window
        of specified width that have scores above 100.
        scores: a list of scores (one for each position in sequence)
        width: width of window
        call_cnt: required number of positions with score 100 or more
        return: list of "calls" (positions in windows with at least call_cnt)
    """
    above = [False for _ in range(len(scores))]
    cnt = 0 # keep track of how many in the current window that are > 100
    for i in range(len(scores)):
        if scores[i] > 100: cnt += 1
        if i >= width:
            if scores[i - width] > 100: cnt -= 1
        if cnt >= call_cnt:
            for j in range(max(0, i - width + 1), i + 1):
                above[j] = True
    return above

def markAvgAbove(scores, width = 4, call_avg = 100.0):
    """ Create a list of booleans that mark all positions within a window of specified width
        that have an average score above specified call_avg.
    """
    above = [False for _ in range(len(scores))]
    sum = 0.0 #
    for i in range(len(scores)):
        sum += scores[i]
        if i >= width: #
            sum -= scores[i - width]
        if sum >= call_avg * width:
            for j in range(max(0, i - width + 1), i + 1):
                above[j] = True
    return above

def extendDownstream(scores, calls, width = 4):
    """ Create a list of booleans that mark all positions that are contained
        in supplied calls list AND extend this list downstream containing a
        specified width average of 100.
    """
    sum = 0.0
    order = list(range(0, len(calls) - 1, +1))  # we are extending calls downstream
    cnt = 0
    for i in order:  # extend to the right
        if calls[i]: # to extend a call is required in the first place
            cnt += 1
            sum += scores[i] # keep a sum to be able to average
            if cnt >= width:   # only average over a width
                sum -= scores[i - width + 1]
            if not calls[i + 1] and sum + scores[i + 1] > width * 100: # check
                calls[i + 1] = True
        else: # no call, reset sum
            cnt = 0
            sum = 0.0
    return calls

def extendUpstream(scores, calls, width = 4):
    """ Create a list of booleans that mark all positions that are contained in supplied calls list
        AND extend this list upstream containing a specified width average of 100.
    """
    sum = 0.0
    order = list(range(len(calls) - 1, 0, -1))  # we are extending calls upstream/to-the-left
    cnt = 0
    for i in order:  # extend to the right
        if calls[i]: # a requirement to extend is to have a call in the first place
            cnt += 1
            sum += scores[i] # keep a sum to be able to average
            if cnt >= width:   # only average over a width
                sum -= scores[i + width - 1]
            if not calls[i - 1] and sum + scores[i - 1] > width * 100: # check average
                calls[i - 1] = True
        else: # no call, reset sum
            cnt = 0
            sum = 0.0
    return calls

def calcRegionAverage(scores, calls):
    """ Determine for each position in a calls list the average score over the region
        in which it is contained.
    """
    region_avg = []
    sum = 0.0
    cnt = 0
    # First determine the average for each region
    for i in range(len(scores)):  # go through each position
        if calls[i]:              # position is part of a "called" region
            sum += scores[i]      # add the score of that position to the average
            cnt += 1              # keep track of the number of positions in the region
        else:                     # we are outside a "called" region
            if cnt > 0:           # if it is the first AFTER a called region
                region_avg.append(sum/cnt)   # save the average
            sum = 0.0             # reset average
            cnt = 0
    if cnt > 0:           # if it is the first AFTER a called region
        region_avg.append(sum/cnt)   # save the average
    # with all averages known, we'll populate the sequence of "averages"
    region = 0
    pos_avg = []
    cnt = 0
    for i in range(len(scores)):
        if calls[i]:
            pos_avg.append(region_avg[region])
            cnt += 1
        else:
            pos_avg.append(0)
            if cnt > 0:
                region += 1
            cnt = 0
    return pos_avg

def checkSupport(calls, diff):
    """ Create a list of booleans indicating if each true position is supported
        by a positive score """
    supported = []
    for i in range(len(calls)):   # go through each position
        supported.append(calls[i] and diff[i] > 0)
    return supported

def getScores(seq, index = 0):
    """ Create a score list for a sequence by referencing the Chou-Fasman table.
    """
    return [cf_dict[s.upper()][index] for s in seq]
