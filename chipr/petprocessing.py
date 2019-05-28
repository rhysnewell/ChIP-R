from chipr import bed
import scipy
import numpy as np
from chipr import rankprod
from chipr import multipletesting

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
    # for bin in range(1, 50):
    #     pets = []
    #     lints = []
    #     ps = []
    #     ds = []
    #     # rp = []
    #     for idx, j in enumerate(np.digitize(dists, bins)):
    #         if j == bin:
    #             pets.append(PETs[idx])
    #             lints.append(listinters[idx])
    #             ps.append(probcor[idx])
    #             ds.append(dists[idx])
                # rp.append(RP[idx])
        # try:
        #     pets, lints, ps, ds, rp = (list(x) for x in zip(*sorted(zip(pets, lints, ps, ds, rp),
        #                                                                  key=lambda pair: pair[0], reverse=True)))
        # except ValueError:
        #     continue
        #
    n = len(RP)
        # k = sum(pets)
    t = sum(PETs)
        # brobs = []
    tags = 0
        # RPranks = scipy.stats.rankdata(RP)
        # mRPrank = max(RPranks)
        # RPranks1 = [mRPrank+1 - x for x in RPranks]
        # probranks = scipy.stats.rankdata(probcor)
        # mprobs = max(probranks)
        # probranks1 = [mprobs+1 -x for x in probranks]
        # # distranks = scipy.stats.rankdata(dists)
        # # mdist = max(distranks)
        # # distranks1 = [mdist+1 - x for x in distranks]
        # # PETranks1 = scipy.stats.rankdata(PETs)
        # # mpet = max(PETranks1)
        # rankprods = []
        # # PETranks = [mpet + 1 - x for x in PETranks1]
        # PETdists = [x/y for x,y in zip(dists, PETs)]
        # PDranks = scipy.stats.rankdata(PETdists)
        # mPD = max(PDranks)
        # PDranks1 = [mPD+1 - x for x in PDranks]
        # for rank in range(len(RP)):
        #     print(PDranks[rank], PETdists[rank])
        #     try:
        #         rankprods.append(PDranks1[rank]*RPranks1[rank])
        #         pets2.append(PETs[rank])
        #     except IndexError:
        #         continue
        # print(rankprods)
        # rpb = rankprodbounds(rankprods, len(rankprods), 2, 'geometric')
        # for p in rpb:
        #     brobs.append(p)
        # print(rpb)
        # mRP = max(rankprods)
        # mPETs = max(PETs)
        # mdists = max(dists)
    for i, p in enumerate(zip(probcor, PETs)):
        tags += p[1]
        b = scipy.stats.binom.cdf(tags, t, 1-p[0])
        # b2 = scipy.stats.combine_pvalues([b, [i]], 'stouffer', [(mPETs-PETs[i])/mPETs, (mRP-rankprods[i])/mRP])
        print(p[0], p[1])
        brobs.append(b)
        pets2.append(p[1])
        ints.append(listinters[i])
    #     # secinters = []
    #     # secPETs = []
    #     # secprobcor = []
    #     print(rpb)
    # ints, probs = (list(x) for x in zip(*sorted(zip(ainters, probs), key=lambda pair: pair[1], reverse=False)))
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


                # else:
                #     secinters.append(i)
                #     secPETs.append(pet)
                #     secprobcor.append(pc)

    # n = sum(secPETs)
    # brobs = []
    # tags = 0
    # for i, p in zip(secPETs, secprobcor):
    #     tags += i
    #     b = scipy.stats.binom.cdf(tags, n, p)
    #     brobs.append(b)
    #
    # corrected = multipletesting.fdrcorrection(brobs, alpha)
    # for i, p in zip(listinters, brobs):
    #     i.addOption(pValue=p)
    #     if p <= alpha:
    #         if p == 0:
    #             i.addOption(pValue=0.00000000000001)
    #         fdr.append(i)


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

