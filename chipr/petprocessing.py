from chipr import bed
import scipy
import numpy as np
from chipr import rankprod
from chipr import multipletesting

"""
Dead code soon to be deleted
"""
def verifyChIAPeaks(chiaData, chipData, filename, alpha):
    listinters = []
    PETs = []
    probs = []
    inters = {}
    PETgroups = {}
    prs = []
    fdr = []
    RP = []
    if getattr(chiaData, 'format') == 'BED12':
        chiaData = bed.BED12toBEDPE(chiaData)
    for link in chiaData:
        if link.PETs <= 2:
            continue
        o1 = chipData.getOverlap(link.partner1)
        o2 = chipData.getOverlap(link.partner2)
        try:
            PETgroups[link.PETs].append(link)
        except KeyError:
            PETgroups[link.PETs] = [link]

        try:
            if len(o1) > 0 and len(o2) > 0:
                for peaks1, peaks2 in zip(o1, o2):
                    # Create entries in dictionary with peaks as support
                    try:

                        inters[link].append([peaks1, peaks2])
                    except KeyError:

                        inters[link] = [peaks1, peaks2]
        except TypeError:
            continue
    dists = []
    for k, v in inters.items():
        pvals = []
        rs = []
        # if k.PETs <= 1:
        #     continue
        # kvals = []
        for peak in v:
            try:
                pvals.append(peak.pValue)
                rs.append(peak.signalValue)

                # combined = scipy.stats.combine_pvalues([peak1.pValue, peak2.pValue], 'fisher')
                # pvals.append(combined[1])
            except AttributeError:
                if type(peak) == 'list':
                    for p in peak:
                        try:
                            pvals.append(p.pValue)
                            rs.append(peak.signalValue)
                        except AttributeError:
                            continue
                else:
                    continue
        if len(pvals) != 0:
            combined = scipy.stats.combine_pvalues(pvals, 'stouffer', rs)
            probs.append(combined[1])
            # rs = iterflatten(rs)
            # RP.append(rs)
            listinters.append(k)
            PETs.append(k.PETs)
            dists.append(k.getDistance())

    # probcor = multipletesting.fdrcorrection(probs, alpha)
    PETs, listinters, probcor, dists = (list(x) for x in zip(*sorted(zip(PETs, listinters, probs, dists),
                                                              key=lambda pair: pair[0], reverse=True)))

    bins = np.arange(min(dists), max(dists), max(dists)/50, dtype=int)
    brobs = []
    ints = []
    pets2 = []
    ainters = []
    tags = 0

    n = len(RP)
        # k = sum(pets)
    t = sum(PETs)
        # brobs = []
    tags = 0

    for i, p in enumerate(zip(probcor, PETs)):
        tags += p[1]
        b = scipy.stats.binom.cdf(tags, t, 1-p[0])
        print(p[0], p[1])
        brobs.append(b)
        pets2.append(p[1])
        ints.append(listinters[i])

    corrected = multipletesting.fdrcorrection(brobs, alpha)
    thresh = rankprod.thresholdCalc(brobs)
    try:
        binomAlpha = round(min(thresh[2]), 3)
        if binomAlpha==0:
            binomAlpha=0.05
    except ValueError:
        binomAlpha = 0.05
    # brobs.append(rpb)
    print(binomAlpha)
    for i, p in zip(ints, corrected[1]):
        i.addOption(pValue=p)
        prs.append(p)
        ainters.append(i)
        # print(p, binomAlpha)
        if p <= 0.05:
            if p == 0:
                i.addOption(pValue=0.00000000000001)
            fdr.append(i)

    fdr = bed.BedFile(fdr, "BEDPE")
    allinters = bed.BedFile(ainters, "BEDPE")
    print(len(fdr), 'Interactions pass FDR threshold')
    bed.writeBedFile(allinters, filename + '.bedpe', format='BEDPE')
    bed.writeBedFile(fdr, filename + '.bedpe', format="BEDPE")
    intersBED12 = bed.BEDPEtoBED12(allinters)
    bed.writeBedFile(intersBED12, filename + '.bed', format='BED12')

    return ainters, fdr, allinters, PETs, dists, brobs, corrected[1]


'''
joinChIA:
Helper function that puts all the interactions from multiple bedPE files in to one file
'''


def joinChIA(chia):
    interactions = []
    for rep in chia:
        for intr in rep:
            interactions.append(intr)

    pe = bed.BedFile(interactions, 'bedpe')
    return pe


def ChIAreproducibility(chia, chip, minentries=None, rank='signalValue', threshold='all',
                        alpha=0.05, filename='NA'):
    if minentries is None:
        minentries = len(chip)
    if len(chia) == 1 and len(chip) == 1:
        return verifyChIAPeaks(chia[0], chip[0], filename=filename, alpha=alpha)
    else:
        try:  # Check to see if multiple chiapet files are entered
            chia[0]
            inters = joinChIA(chia)
        except TypeError:
            inters = chia
        rankprod.performrankprod(chip, minentries=minentries, rankmethod=rank, alpha=alpha, filename=filename+'.bed')
        if threshold == 'binom':
            RP = bed.BedFile('T2_' + filename + '.bed', 'Peaks')

            return verifyChIAPeaks(inters, RP, filename=filename, alpha=alpha)
        elif threshold == 'alpha':
            RP = bed.BedFile('T1_' + filename + '.bed', 'Peaks')
            return verifyChIAPeaks(inters, RP, filename=filename, alpha=alpha)

        elif threshold == 'all':
            RP = bed.BedFile('ALL_' + filename + '.bed', 'Peaks')
            return verifyChIAPeaks(inters, RP, filename=filename, alpha=alpha)

