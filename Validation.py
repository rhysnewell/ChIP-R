
# from statsmodels.stats import multitest
# from scipy.interpolate import UnivariateSpline
# import statsmodels.api as sm
# import importlib
print("Loading modules...")
import bed
import prob
import sequence
import numpy as np
import pandas as pd
import os
import math
#import matplotlib.pyplot as plt
import RankProd
import twobit
import sym
import multipletesting


# import BedPairFile
# import scipy


print("Loading dorksouls...")
def dorksouls(directory, id=''):
    cnt = 1
    for p in os.listdir(directory):
        p = os.fsdecode(p)
        # if p.endswith('.bed'):
            # newdir = directory+'/'+p
            # for file in os.listdir(newdir):
            #     filename = os.fsdecode(file)
        if p.endswith('.fa'):
            os.system(
                '/home/rhys/meme/libexec/meme-5.0.1/fasta-get-markov -m 1 ' + directory + '/' + p + ' > temp.txt')
            os.system(str('/home/rhys/meme/bin/fimo --bfile temp.txt --max-stored-scores 1000000 --thresh 1 --oc '+
                          directory+p+' --motif '+id+' '+directory+'/MEME_JASPAR.txt '+directory+'/'+p))


            cnt += 1

print("Sorter")
def thegreatsorter(bedf, use_score=False, use_signal=False):
    peaks = []
    signals = []
    sv = False
    for i in bedf:
        peaks.append(i)
        if use_score is True:
            sig = i.score
        elif use_signal is True:

            sig = i.signalValue
        else:
            try:
                sig = i.pValue
            except AttributeError:
                sig = 0.000000001
        if sig == -1:
            sv = True
            signals.append(i.signalValue)
        if sig > 1:
            sv = True
            signals.append(sig)
        else:
            signals.append(sig)
    print(sv)
    if sv is True:
        sortedpeaks = [x for _, x in sorted(zip(signals, peaks), reverse=True, key=lambda pair: pair[0])]
    else:
        sortedpeaks = [x for _, x in
                       sorted(zip(signals, peaks), reverse=False, key=lambda pair: pair[0])]
    # print(sortedpeaks[100].pValue, sortedpeaks[200].pValue, min(signals), max(signals))
    return sortedpeaks

print("Loading motif_discovery...")
def motif_discovery(peak, refgenome=None, JASPAR_ID=None, multicounts=None, pwm1=None, pwm2=None):
    # if JASPAR_ID is None:
    # return print("Please provide a JASPAR_ID for your motif")
    # hg = twobit.TwoBitFile(refgenome)
    # d = prob.readMultiCounts('JASPAR_matrices.txt')
    # pwm1 = sequence.PWM(d[JASPAR_ID])
    # pwm2 = pwm1.getRC()

    hg = refgenome
    # d = multicounts
    pwm1 = pwm1
    pwm2 = pwm2
    # seq = input
    seqs = []
    cnt=1
    # for b in entries:
    try:
        seqstr = hg[peak.chrom][peak.chromStart:peak.chromEnd]
        newseq = sequence.Sequence(seqstr.upper(), sym.DNA_Alphabet, peak.chrom+':'+str(peak.chromStart)+'-'+str(peak.chromEnd))
        # seq_len.append(len(newseq))
        cnt += 1
        seqs.append(newseq)
    except:
        print('Failed to map', peak.chrom, peak.chromStart, peak.chromEnd)
    # seqs200 = []
    # for seq in seqs:
    #     centre = len(seq) // 2
    #     seqs200.append(sequence.Sequence(seq[centre - 100:centre + 100], sym.DNA_Alphabet, str(seq.name)))
    # sequence.writeFastaFile('Max200.fa', seqs200)

    # seq_len = 200
    # for seq in seqs200:
    #     if len(seq) != seq_len:
    #         print('Error: all sequences must be of the same length, and', seq.name, 'is', len(seq), 'and not', seq_len)
    #         break

    """
        Produce a plot for the "average scan" of the specified motif.
        The plot has as its x-axis position of sequence, and
        the y-axis the average PWM score over all sequences.
        Make sure the following variables are set correctly before running:
            seqs200 - all sequences to be scanned, must be same lengths
            pwm1 - PWM
            pwm2 - PWM reverse strand
    """
    threshold = 0
    # seq_len = len(seqs200[0])
    avg_motif_score1 = []
    i_seq = 0
    motif_width = pwm1.length

    for seq in seqs:
        i_seq += 1
        seq_len = len(seq)
        hits = pwm1.search(seq, threshold)
        pos_scores = [0] * seq_len
        for hit in hits:
            # mark hit at *center* of site (hence motif_width//2)
            pos_scores[hit[0] + (motif_width // 2)] = hit[2]
        # negative strand
        # print(pos_scores[0])
        hits = pwm2.search(seq, threshold)
        neg_scores = [0] * seq_len
        for hit in hits:
            neg_scores[hit[0] + (motif_width // 2)] = hit[2]
        # for each position use the maximum score of the two strands
        # print(neg_scores[0])
        seq_scores = []
        for i in range(seq_len):
            p = pos_scores[i]
            n = neg_scores[i]
            score = max(p, n)
            if score > 0:
                seq_scores.append(score)
            else:
                seq_scores.append(0.0)
            # if score > threshold:
            #     avg_motif_score1[i] += score
        # print(seq_scores)
        if len(seq_scores) > 0:
            avg_motif_score1.append(max(seq_scores))
    # compute average score
    # print(avg_motif_score1)
    if len(avg_motif_score1) > 0:
        n = len(avg_motif_score1)
        # avg_max_score = np.mean(avg_motif_score1)
        # std_max_score = np.std(avg_motif_score1)
        return avg_motif_score1, avg_motif_score1, n
    else:
        avg_max_score = 0
        std_max_score = 0
        return 0.0, 0.0, 1
    # for i in range(seq_len):
    #     avg_motif_score1[i] /= len(seqs200)



print("Electricboogaloo...")
def electricboogaloo(directory, JASPAR_ID=None, refgenome='hg19.2bit', output='scores', idr_only=False):
    hg = twobit.TwoBitFile(refgenome)
    d = prob.readMultiCounts('JASPAR_matrices.txt')
    pwm1 = sequence.PWM(d[JASPAR_ID])
    pwm2 = pwm1.getRC()
    cnt = 0
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith('_T1.bed'):
            continue
        elif filename.endswith('_T2.bed'):
            continue
        elif filename.endswith('_idrValues.txt'):
            # continue
            peaks = bed.BedFile(directory + "/" + filename, format='2idr')
            sortedpeaks = thegreatsorter(peaks)
            # os.system('mkdir '+directory+'/fastas/'+filename)
            # tops = range(0, len(sortedpeaks), len(sortedpeaks))
            # for i in range(len(tops)):
            #     try:
            scores = []
            # scores.append(motif_discovery(sortedpeaks, refgenome=hg,
            #                             JASPAR_ID=JASPAR_ID, multicounts=d, pwm1=pwm1, pwm2=pwm2))
            for peaks in sortedpeaks:
                # cnt += 1
                # # if peaks > 50000:
                # #     break
                # if cnt == 100:
                #     print(peaks, ':', peaks + 100)
                #     cnt = 0
                try:
                    scores.append(motif_discovery(peaks, refgenome=hg,
                                                  JASPAR_ID=JASPAR_ID, multicounts=d, pwm1=pwm1, pwm2=pwm2))
                except IndexError:
                    scores.append(motif_discovery(peaks, refgenome=hg,
                                                  JASPAR_ID=JASPAR_ID, multicounts=d, pwm1=pwm1, pwm2=pwm2))
            df = pd.DataFrame(scores, columns=['avg_max', 'std', 'n'])
            df.to_csv(directory + '/' + output + '/' + filename.split('.')[0] + '.csv')
        elif filename.endswith('_ALL.bed') and idr_only is False:
            # continue
            input_peaks=bed.BedFile(directory+"/"+filename, format='idr')
            sortedpeaks = thegreatsorter(input_peaks)
            # # os.system('mkdir '+directory+'/fastas/'+filename)
            # # tops = range(0, len(sortedpeaks), len(sortedpeaks))
            # # for i in range(len(tops)):
            # #     try:
            scores = []
            # scores.append(motif_discovery(sortedpeaks, refgenome=hg,
            #                               JASPAR_ID=JASPAR_ID, multicounts=d, pwm1=pwm1, pwm2=pwm2))
            for peaks in sortedpeaks:
                # cnt += 1
                # # if peaks > 50000:
                # #     break
                # if cnt == 100:
                #     print(peaks, ':', peaks + 100)
                #     cnt = 0
                try:
                    scores.append(motif_discovery(peaks, refgenome=hg,
                                                  JASPAR_ID=JASPAR_ID, multicounts=d, pwm1=pwm1, pwm2=pwm2))
                except IndexError:
                    scores.append(motif_discovery(peaks, refgenome=hg,
                                                  JASPAR_ID=JASPAR_ID, multicounts=d, pwm1=pwm1, pwm2=pwm2))
            # print(scores)
            df = pd.DataFrame(scores, columns=['avg_max', 'std', 'n'])
            df.to_csv(directory+'/'+output+'/'+filename.split('.')[0]+'.csv')
        elif filename.endswith('.narrowPeak'):
            # continue
            input_peaks=bed.BedFile(directory+"/"+filename, format='Peaks')
            sortedpeaks = thegreatsorter(input_peaks)
            # # os.system('mkdir '+directory+'/fastas/'+filename)
            # # tops = range(0, len(sortedpeaks), len(sortedpeaks))
            # # for i in range(len(tops)):
            # #     try:
            scores = []
            # scores.append(motif_discovery(sortedpeaks, refgenome=hg,
            #                               JASPAR_ID=JASPAR_ID, multicounts=d, pwm1=pwm1, pwm2=pwm2))
            for peaks in sortedpeaks:
                # cnt += 1
                # # if peaks > 50000:
                # #     break
                # if cnt == 100:
                #     print(peaks, ':', peaks + 100)
                #     cnt = 0
                try:
                    scores.append(motif_discovery(peaks, refgenome=hg,
                                                  JASPAR_ID=JASPAR_ID, multicounts=d, pwm1=pwm1, pwm2=pwm2))
                except IndexError:
                    scores.append(motif_discovery(peaks, refgenome=hg,
                                                  JASPAR_ID=JASPAR_ID, multicounts=d, pwm1=pwm1, pwm2=pwm2))
            df = pd.DataFrame(scores, columns=['avg_max', 'std', 'n'])
            df.to_csv(directory+'/'+output+'/'+filename.split('.')[0]+'.csv')

        elif filename.endswith('_bf.bed'):
            # continue
            input_peaks=bed.BedFile(directory+"/"+filename, format='bed6')
            sortedpeaks = thegreatsorter(input_peaks, use_score=True)
            # # os.system('mkdir '+directory+'/fastas/'+filename)
            # # tops = range(0, len(sortedpeaks), len(sortedpeaks))
            # # for i in range(len(tops)):
            # #     try:
            scores = []
            # scores.append(motif_discovery(sortedpeaks, refgenome=hg,
            #                               JASPAR_ID=JASPAR_ID, multicounts=d, pwm1=pwm1, pwm2=pwm2))
            for peaks in sortedpeaks:
                # cnt += 1
                # # if peaks > 50000:
                # #     break
                # if cnt == 100:
                #     pr
                #     cnt = 0
                try:
                    scores.append(motif_discovery(peaks, refgenome=hg,
                                                  JASPAR_ID=JASPAR_ID, multicounts=d, pwm1=pwm1, pwm2=pwm2))
                except IndexError:
                    scores.append(motif_discovery(peaks, refgenome=hg,
                                                  JASPAR_ID=JASPAR_ID, multicounts=d, pwm1=pwm1, pwm2=pwm2))
            df = pd.DataFrame(scores, columns=['avg_max', 'std', 'n'])
            df.to_csv(directory+'/'+output+'/'+filename.split('.')[0]+'.csv')

        elif filename.endswith('_mspc.bed'):
            # continue
            input_peaks=bed.BedFile(directory+"/"+filename, format='mspc')
            sortedpeaks = thegreatsorter(input_peaks, use_signal=True)
            # # os.system('mkdir '+directory+'/fastas/'+filename)
            # # tops = range(0, len(sortedpeaks), len(sortedpeaks))
            # # for i in range(len(tops)):
            # #     try:
            scores = []
            # scores.append(motif_discovery(sortedpeaks, refgenome=hg,
            #                               JASPAR_ID=JASPAR_ID, multicounts=d, pwm1=pwm1, pwm2=pwm2))
            for peaks in sortedpeaks:
                # cnt += 1
                # # if peaks > 50000:
                # #     break
                # if cnt == 100:
                #     print(peaks, ':', peaks + 100)
                #     cnt = 0
                try:
                    scores.append(motif_discovery(peaks, refgenome=hg,
                                                  JASPAR_ID=JASPAR_ID, multicounts=d, pwm1=pwm1, pwm2=pwm2))
                except IndexError:
                    scores.append(motif_discovery(peaks, refgenome=hg,
                                                  JASPAR_ID=JASPAR_ID, multicounts=d, pwm1=pwm1, pwm2=pwm2))
            df = pd.DataFrame(scores, columns=['avg_max', 'std', 'n'])
            df.to_csv(directory+'/'+output+'/'+filename.split('.')[0]+'.csv')

# def mammamia(directory):
#     tss = bed.BedFile('NAR_Runs/refGene_hg19_TSS.bed', 'TSS')
#     for file in os.listdir(directory):
#         filename = os.fsdecode(file)
#         # if filename.endswith('_ALL.bed'):
#         #     pass
#         if filename.endswith('.bed'):
#             peaks=bed.BedFile(directory+"/"+filename, format='idr')
#             sortedpeaks = thegreatsorter(peaks)
#             o = overlap(peaks, tss)
#             print('Overlap for '+filename+':')
#             print(o[1], o[2], len(peaks))

def getMatrixMax(JASPAR_ID='MA0058.1'):
    d = prob.readMultiCounts('JASPAR_matrices.txt')
    pwm1 = sequence.PWM(d[JASPAR_ID])
    pwm2 = pwm1.getRC()
    mat1 = pwm1.getMatrix()
    mat2 = pwm2.getMatrix()
    maxsum = 0
    maxsum_rc = 0
    for i in range(len(mat1[0])):
        maxsum += max(mat1[:, i])
        maxsum_rc += max(mat2[:, i])
    print("max sum:", maxsum)
    return maxsum


def singleremix(filename, JASPAR_ID=None, refgenome='hg19.2bit', output='scores'):
    hg = twobit.TwoBitFile(refgenome)
    d = prob.readMultiCounts('JASPAR_matrices.txt')
    pwm1 = sequence.PWM(d[JASPAR_ID])
    pwm2 = pwm1.getRC()
    print(pwm1.getMatrix())
    print(pwm2.getMatrix())
    input_peaks = bed.BedFile(filename, format='idr')
    sortedpeaks = thegreatsorter(input_peaks)
    scores = []
    for peaks in range(0, len(sortedpeaks), 100):
        if peaks >= len(sortedpeaks):
            break
        # print(peaks, ':', peaks + 100)
        try:
            scores.append(motif_discovery(sortedpeaks[peaks:peaks + 100], refgenome=hg,
                                          JASPAR_ID=JASPAR_ID, multicounts=d, pwm1=pwm1, pwm2=pwm2))
        except IndexError:
            try:
                scores.append(motif_discovery(sortedpeaks[peaks:len(sortedpeaks)], refgenome=hg,
                                              JASPAR_ID=JASPAR_ID, multicounts=d, pwm1=pwm1, pwm2=pwm2))
            except IndexError:
                print('frick')
                break
    # print(scores)
    df = pd.DataFrame(scores, columns=['avg_max', 'std', 'n'])
    df.to_csv(filename.split('.')[0] + '.csv')


def runthatcode(bedfs, thename='default', rankmethod='pvalue'):


    RankProd.performrankprod([bedfs[0], bedfs[1]], minentries=1, rankmethod=rankmethod, filename=thename + '_P1_min1')
    RankProd.performrankprod([bedfs[0], bedfs[2]], minentries=1, rankmethod=rankmethod, filename=thename + '_P2_min1')
    RankProd.performrankprod([bedfs[0], bedfs[3]], minentries=1, rankmethod=rankmethod, filename=thename + '_P3_min1')
    RankProd.performrankprod([bedfs[1], bedfs[2]], minentries=1, rankmethod=rankmethod, filename=thename + '_P4_min1')
    RankProd.performrankprod([bedfs[1], bedfs[3]], minentries=1, rankmethod=rankmethod, filename=thename + '_P5_min1')
    RankProd.performrankprod([bedfs[2], bedfs[3]], minentries=1, rankmethod=rankmethod, filename=thename + '_P6_min1')

    RankProd.performrankprod([bedfs[0], bedfs[1]], minentries=2, rankmethod=rankmethod, filename=thename + '_P1_min2')
    RankProd.performrankprod([bedfs[0], bedfs[2]], minentries=2, rankmethod=rankmethod, filename=thename + '_P2_min2')
    RankProd.performrankprod([bedfs[0], bedfs[3]], minentries=2, rankmethod=rankmethod, filename=thename + '_P3_min2')
    RankProd.performrankprod([bedfs[1], bedfs[2]], minentries=2, rankmethod=rankmethod, filename=thename + '_P4_min2')
    RankProd.performrankprod([bedfs[1], bedfs[3]], minentries=2, rankmethod=rankmethod, filename=thename + '_P5_min2')
    RankProd.performrankprod([bedfs[2], bedfs[3]], minentries=2, rankmethod=rankmethod, filename=thename + '_P6_min2')

    RankProd.performrankprod([bedfs[0], bedfs[1], bedfs[2]], minentries=1, rankmethod=rankmethod, filename=thename + '_T1_min1')
    RankProd.performrankprod([bedfs[0], bedfs[1], bedfs[3]], minentries=1, rankmethod=rankmethod, filename=thename + '_T2_min1')
    RankProd.performrankprod([bedfs[0], bedfs[2], bedfs[3]], minentries=1, rankmethod=rankmethod, filename=thename + '_T3_min1')
    RankProd.performrankprod([bedfs[1], bedfs[2], bedfs[3]], minentries=1, rankmethod=rankmethod, filename=thename + '_T4_min1')

    RankProd.performrankprod([bedfs[0], bedfs[1], bedfs[2]], minentries=2, rankmethod=rankmethod,
                             filename=thename + '_T1_min2')
    RankProd.performrankprod([bedfs[0], bedfs[1], bedfs[3]], minentries=2, rankmethod=rankmethod,
                             filename=thename + '_T2_min2')
    RankProd.performrankprod([bedfs[0], bedfs[2], bedfs[3]], minentries=2, rankmethod=rankmethod,
                             filename=thename + '_T3_min2')
    RankProd.performrankprod([bedfs[1], bedfs[2], bedfs[3]], minentries=2, rankmethod=rankmethod,
                             filename=thename + '_T4_min2')

    RankProd.performrankprod([bedfs[0], bedfs[1], bedfs[2]], minentries=3, rankmethod=rankmethod,
                             filename=thename + '_T1_min3')
    RankProd.performrankprod([bedfs[0], bedfs[1], bedfs[3]], minentries=3, rankmethod=rankmethod,
                             filename=thename + '_T2_min3')
    RankProd.performrankprod([bedfs[0], bedfs[2], bedfs[3]], minentries=3, rankmethod=rankmethod,
                             filename=thename + '_T3_min3')
    RankProd.performrankprod([bedfs[1], bedfs[2], bedfs[3]], minentries=3, rankmethod=rankmethod,
                             filename=thename + '_T4_min3')

    RankProd.performrankprod(bedfs, minentries=1, rankmethod=rankmethod, filename=thename + '_min1')
    RankProd.performrankprod(bedfs, minentries=2, rankmethod=rankmethod, filename=thename + '_min2')
    RankProd.performrankprod(bedfs, minentries=3, rankmethod=rankmethod, filename=thename + '_min3')
    RankProd.performrankprod(bedfs, minentries=4, rankmethod=rankmethod, filename=thename + '_min4')

def runthreerep(thelist, thename='default'):
    cnt1 = 0
    cnt2 = 2
    cnt3 = 1
    for file in thelist:
        if cnt3 > 2:
            break
        else:
            RankProd.performrankprod([thelist[0],thelist[1], thelist[cnt2]], minentries=3, rankmethod='pvalue',
                                     filename=thename + '1_3repmin3_' + str(cnt3))
            RankProd.performrankprod([thelist[2], thelist[3], thelist[cnt1]], minentries=3, rankmethod='pvalue',
                                     filename=thename + '2_3repmin3_' + str(cnt3))
            # RankProd.performrankprod(newlist, minentries=1, filename=thename + '_2repmin1_' + str(cnt3))
            cnt1 += 1
            cnt2 += 1
            cnt3 += 1



def dorksouls2(directory):
    #Collects scores from FIMO output
    cnt = 1
    scores = {}
    for p in os.listdir(directory):
        p = os.fsdecode(p)
        # if p.endswith('.bed'):
            # newdir = directory+'/'+p
            # for file in os.listdir(newdir):
            #     filename = os.fsdecode(file)
        if p.endswith('.fa'):
            f = open(directory+'/'+p+"/fimo.tsv")
            scores[p.split('.')[0]] = []
            cnt = 0
            for line in f:
                l = line.split('\t')
                try:
                    if l[6] != 'score':
                        scores[p.split('.')[0]].append(float(l[6]))
                        cnt += 1
                        if cnt == 1000:
                            break
                except IndexError:

                    if cnt != 999:
                        print(cnt)
                        print(p.split('.')[0], 'not sufficient scores')
                    break
    return scores

def getIDRPeaks(directory):
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith('_idrValues.txt'):
            peaks = bed.BedFile(directory + "/" + filename, format='2idr')
            qvalue = 0
            pvalue = 0
            for peak in peaks:
                if peak.pValue <= 0.05:
                    pvalue += 1
                if peak.qValue <= 0.05:
                    qvalue += 1
            print(filename + ' pvalues under 0.05: ' + str(pvalue))
            print(filename + ' qvalues under 0.05: ' + str(qvalue))

# print("Loading MAX files...")
# max1 = bed.BedFile('NAR_Runs/BAMS/MAX/max1_snyder1_peaks.narrowPeak', 'Peaks')
# max2 = bed.BedFile('NAR_Runs/BAMS/MAX/max1_snyder2_peaks.narrowPeak', 'Peaks')
# max3 = bed.BedFile('NAR_Runs/BAMS/MAX/max2_myers1_peaks.narrowPeak', 'Peaks')
# max4 = bed.BedFile('NAR_Runs/BAMS/MAX/max2_myers2_peaks.narrowPeak', 'Peaks')
# maxlist = [max1, max2, max3, max4]
#
# print("Loading REST files...")
# rest1 = bed.BedFile('NAR_Runs/BAMS/REST/rest1_snyder1_peaks.narrowPeak', 'Peaks')
# rest2 = bed.BedFile('NAR_Runs/BAMS/REST/rest1_snyder2_peaks.narrowPeak', 'Peaks')
# rest3 = bed.BedFile('NAR_Runs/BAMS/REST/rest2_myers1_peaks.narrowPeak', 'Peaks')
# rest4 = bed.BedFile('NAR_Runs/BAMS/REST/rest2_myers2_peaks.narrowPeak', 'Peaks')
# rest = [rest1, rest2, rest3, rest4]
#
# print("Loading SRF files...")
# srf3 = bed.BedFile('NAR_Runs/BAMS/SRF/srf2_myers1_peaks.narrowPeak', 'Peaks')
# srf4 = bed.BedFile('NAR_Runs/BAMS/SRF/srf2_myers2_peaks.narrowPeak', 'Peaks')
# srf1 = bed.BedFile('NAR_Runs/BAMS/SRF/srf1_myers1_peaks.narrowPeak', 'Peaks')
# srf2 = bed.BedFile('NAR_Runs/BAMS/SRF/srf1_myers2_peaks.narrowPeak', 'Peaks')
# srf = [srf1, srf2, srf3, srf4]
#
# nfi1 = bed.BedFile('NFI_chip/CD1_P7_NfiA_peaks.narrowPeak', 'Peaks')
# nfi2 = bed.BedFile('NFI_chip/CD1_P7_NfiB_peaks.narrowPeak', 'Peaks')
# nfi3 = bed.BedFile('NFI_chip/CD1_P7_NfiX_peaks.narrowPeak', 'Peaks')
# nfi4 = bed.BedFile('NFI_chip/Mixed_P7_NfiA_peaks.narrowPeak', 'Peaks')
# nfi5 = bed.BedFile('NFI_chip/Mixed_P7_NfiB_peaks.narrowPeak', 'Peaks')
# nfi6 = bed.BedFile('NFI_chip/Mixed_P7_NfiX_peaks.narrowPeak', 'Peaks')
# nfi = [nfi1, nfi2, nfi3, nfi4, nfi5, nfi6]
# #
# nf = RankProd.performrankprod(nfi, minentries=2, filename='nfi_test')
#
#

# print("Performing RP analysis...")
# runthatcode(maxlist, 'NAR_Runs/MAX/maxRP')
# runthatcode(rest, 'NAR_Runs/REST/restRP')
# runthatcode(srf, 'NAR_Runs/SRF/srfRP')

# runthreerep(maxlist, 'NAR_Runs/MAX/maxRP')
# runthreerep(rest, 'NAR_Runs/REST/restRP')
# runthreerep(srf, 'NAR_Runs/SRF/srfRP')

print("Motif Analysis...")
print("MAX")
electricboogaloo('NAR_Runs/MAX', 'MA0058.2')
print("REST")
electricboogaloo('NAR_Runs/REST', 'MA0138.2')
print("SRF")
electricboogaloo('NAR_Runs/SRF', 'MA0083.2')



# getIDRPeaks('NAR_Runs/MAX')
# getIDRPeaks('NAR_Runs/REST')
# getIDRPeaks('NAR_Runs/SRF')
#
# dorksouls('/home/rhys/binfpy/src/NAR_Runs/MAX/fastas', 'MA0058.1')
# dorksouls('/home/rhys/binfpy/src/NAR_Runs/REST/fastas', 'MA0138.1')
# dorksouls('/home/rhys/binfpy/src/NAR_Runs/SRF/fastas', 'MA0083.1')
#
# mammamia('NAR_Runs/MAX')
# mammamia('NAR_Runs/REST')
# mammamia('NAR_Runs/SRF')
#
# maxscores = dorksouls2('NAR_Runs/MAX')
# restscores = dorksouls2('NAR_Runs/REST')
# srfscores = dorksouls2('NAR_Runs/SRF')
#
# m = pd.DataFrame(maxscores)
# m.to_csv('~/Documents/Honours/NAR_dfs/maxscores.csv')
# r = pd.DataFrame(restscores)
# r.to_csv('~/Documents/Honours/NAR_dfs/restscores.csv')
# s = pd.DataFrame(srfscores)
# s.to_csv('~/Documents/Honours/NAR_dfs/srfscores.csv')
#
#
# def simRank(x, ncols):
#     allranks = []
#     for i in range(ncols):
#         ranks = scipy.stats.rankdata(x[:, i])
#         ranks = [len(ranks) - r for r in ranks]
#         allranks.append(ranks)
#     return allranks
#
#
# def simProduct(x):
#     zipd = zip(*x)
#     products = [RankProd.prod(i) for i in zipd]
#     return products
#
#
# def collapse(bedf, idr=False):
#     allOl = []
#     pvals = []
#     if idr is True:
#         for chrom in bedf.chroms:
#             chromarr = bedf.generate(chrom)
#             ol1 = []
#             ol2 = []
#             ol3 = []
#             ol4 = []
#
#             for entry in chromarr:
#
#                 if len(ol3) == 0:
#                     cnt = 1
#                     ol1.append(entry.chromStart)
#                     ol2.append(entry.chromEnd)
#                     ol3.append(entry)
#                     try:
#                         ol4.append(entry.score)
#                     except AttributeError:
#                         ol4.append(0)
#
#                 else:
#                     for o in ol3:
#                         if overlaps(o, entry):
#                             ol1.append(entry.chromStart)
#                             ol2.append(entry.chromEnd)
#                             ol3.append(entry)
#                             try:
#                                 ol4.append(entry.score)
#                             except AttributeError:
#                                 ol4.append(0)
#
#                             break
#                         else:
#                             cnt += 1
#                 if cnt > len(ol3):
#                     if len(ol3) > 1:
#                         mer = bed.BedEntry(chrom, min(ol1), max(ol2))
#                         newScore = max(ol4)  # Calculate combined p-value
#
#                         mer.addOption(name='TBD',
#                                       score=newScore,
#                                       strand='.',
#                                       signalValue=-1,
#                                       pValue=-1,
#                                       qValue=-1)
#                         allOl.append(mer)
#
#                     else:
#                         ol3[0].addOption(name='TBD',
#                                          score=newScore,
#                                          strand='.',
#                                          signalValue=-1,
#                                          pValue=-1,
#                                          qValue=-1)
#                         allOl.append(ol3[0])
#
#                     ol1 = []
#                     ol2 = []
#                     ol3 = []
#                     ol4 = []
#
#                     cnt = 1
#                     ol1.append(entry.chromStart)
#                     ol2.append(entry.chromEnd)
#                     ol3.append(entry)
#                     try:
#                         ol4.append(entry.score)
#                     except AttributeError:
#                         ol4.append(0)
#
#         return allOl
#     else:
#         for chrom in bedf.chroms:
#             chromarr = bedf.generate(chrom)
#             ol1 = []
#             ol2 = []
#             ol3 = []
#             ol4 = []
#             ol5 = []
#
#             for entry in chromarr:
#
#                 if len(ol3) == 0:
#                     cnt = 1
#                     ol1.append(entry.chromStart)
#                     ol2.append(entry.chromEnd)
#                     ol3.append(entry)
#                     try:
#                         ol4.append(entry.pValue)
#                     except AttributeError:
#                         ol4.append(0)
#                     try:
#                         ol5.append(entry.qValue)
#                     except AttributeError:
#                         ol5.append(0)
#
#                 else:
#                     for o in ol3:
#                         if overlaps(o, entry):
#                             ol1.append(entry.chromStart)
#                             ol2.append(entry.chromEnd)
#                             ol3.append(entry)
#                             try:
#                                 ol4.append(entry.pValue)
#                             except AttributeError:
#                                 ol4.append(0)
#                             try:
#                                 ol5.append(entry.qValue)
#                             except AttributeError:
#                                 ol5.append(0)
#                             break
#                         else:
#                             cnt += 1
#                 if cnt > len(ol3):
#                     if len(ol3) > 1:
#                         mer = bed.BedEntry(chrom, min(ol1), max(ol2))
#                         newp = scipy.stats.combine_pvalues(ol4)[1]  # Calculate combined p-value
#                         pvals.append(newp)
#
#                         mer.addOption(name='TBD',
#                                       score=min([abs(int(125 * math.log2(newp))), 1000]),
#                                       strand='.',
#                                       signalValue=-1,
#                                       pValue=newp,
#                                       qValue=-1)
#                         allOl.append(mer)
#
#                     else:
#                         ol3[0].addOption(name='TBD',
#                                          score=min([abs(int(125 * math.log2(ol4[0]))), 1000]),
#                                          strand='.',
#                                          signalValue=-1,
#                                          pValue=ol4[0],
#                                          qValue=-1)
#                         allOl.append(ol3[0])
#                         pvals.append(ol4[0])
#                     ol1 = []
#                     ol2 = []
#                     ol3 = []
#                     ol4 = []
#                     ol5 = []
#                     cnt = 1
#                     ol1.append(entry.chromStart)
#                     ol2.append(entry.chromEnd)
#                     ol3.append(entry)
#                     try:
#                         ol4.append(entry.pValue)
#                     except AttributeError:
#                         ol4.append(0)
#                     try:
#                         ol5.append(entry.qValue)
#                     except AttributeError:
#                         ol5.append(0)
#
#         return allOl, pvals
#
#
# """
# Function for calculating a suitable cutoff threshold for a list of P-values.
# Uses the binomial distribution to search for a minimum Pk, where Pk is the sum of binomial PMFs of 1 - p-values.
#
# params:
# @p - A list of P-values from the rank product distribution
# @niter - The number of times the same Pk value has to be seen before the function decides to stop
# @startingIndex - The starting position of to start looking through the p-value list. Lower numbers will provide more
#                  accurate results at the cost of speed.
#
# """
#
#
# def thresholdCalc(p, rev=False):
#     n = len(p)
#     ps = sorted(p, reverse=rev)
#     bs = []
#     for i, v in enumerate(ps):
#         b = scipy.stats.binom.cdf(i, n, v)
#         bs.append(b)
#
#     its = []
#     for p, b in zip(ps, bs):
#         if round(p, 1) == round(b, 1):
#             print(p, b)
#             its.append(ps.index(p))
#
#     return bs, ps, its
#
#
# def diversityCalc(ps, q=1):
#     # ps = [1-p for p in ps]
#     sumps = sum(ps)
#     props = [x / sumps for x in ps]
#     if q == 1:
#         lnps = []
#         for s in props:
#             lnps.append(s * math.log(s))
#         D = math.exp(-sum(lnps))
#     return D
#
#
# def matrixRanks(matrix):
#     ranks1 = [[]] * matrix.shape[1]
#     ranks2 = [[]] * matrix.shape[0]
#
#     for i in range(matrix.shape[1]):
#         r = (len(matrix[:, i]) + 1) - scipy.stats.rankdata(matrix[:, i]).astype(int)
#         ranks1[i] = r
#
#     for k in range(len(ranks1[0])):
#         for j in range(len(ranks1)):
#             ranks2[k].append(ranks1[j][k])
#     # prod = reduce(mul, ranks)
#     prod = [reduce(mul, i) for i in ranks2]
#     # print(prod)
#     # print(matrix.shape[0], matrix.shape[1])
#     pvals = RankProd.rankprodbounds(prod, matrix.shape[0], matrix.shape[1], 'geometric')
#     return pvals
#
#
# def retrieveRanks(ranks):
#     r = []
#     for i in range(len(ranks[1])):
#         rj = []
#         for j in ranks[1][i]:
#             # print(j)
#             rj.append(j[1])
#         r.append(rj)
#
#     return r
#
#
# def wilcoxonRankSign(ranks1, ranks2):
#     mad = []
#     sgn = []
#     for i in range(len(ranks1)):
#         diff = ranks2[i] - ranks1[i]
#         if diff != 0:
#             mad.append(abs(diff))
#             if diff > 0:
#                 sgn.append(1)
#             elif diff < 0:
#                 sgn.append(-1)
#         else:
#             mad.append(0)
#             sgn.append(0)
#             print(i)
#     Nr = len(mad)
#     R = scipy.stats.rankdata(mad, 'average').tolist()
#     W = []
#     for j in range(Nr):
#         W.append(sgn[j] * R[j])
#     return W, mad
#
#
# def MAD(zippedranks, n=4):
#     mad = []
#     sgn = []
#     sad = []
#     print(n)
#     for i in zippedranks:
#         abdiff = []
#         diff = []
#         m = np.mean(i)
#         n = len(i)
#         for j in i:
#             diff.append(j - m)
#             abdiff.append(abs(j - m))
#
#         sad.append((sum(diff)) / n)
#         sumdiff = sum(diff)
#         med = np.median(abdiff)
#         mad.append((sum(abdiff)) / n)
#         if sumdiff < 0:
#             sgn.append(-1)
#
#         elif sumdiff > 0:
#             sgn.append(1)
#         else:
#             sgn.append(0)
#
#     Nr = len(mad)
#     R = scipy.stats.rankdata(mad, 'average').tolist()
#     W = []
#     for j in range(Nr):
#         W.append(sgn[j] * R[j])
#
#     return sad, mad, W
#

#
# def ROC(classified_peaks, reps, IDR=None, minentries=2, alpha=np.arange(0.001, 0.01, 0.001), maxentries='maximum',
#         fn="RP.bed"):
#     if IDR is None:
#         # First create union and rank the entries in each replicate
#         ranks = RankProd.rankreps(reps, minentries, rankmethod='signalValue', duphandling='random', random_seed='0.5',
#                                   specifyMax=maxentries)
#
#         # Calculate rank product for each entry that contributes to a union entry
#         rp = RankProd.rankprod(ranks, method='all')
#
#         # Calculate the pvalues of the rank product values
#         rpb_up = RankProd.rankprodbounds(rp[1], len(rp[1]), len(reps), 'geometric')
#         ent = RankProd.rankEntropy(ranks)
#         # sdiversity = math.exp(ent)
#
#         TPR = []
#         FPR = []
#
#         for i in alpha:
#
#             idxs_t1 = []
#
#             # Calculate Maximum possible entropy value
#
#             for idx, val in enumerate(rpb_up):
#                 if val <= i:
#                     idxs_t1.append(idx)
#
#             t1_unions = []
#             for i in range(len(ranks[0])):
#                 if i in idxs_t1:
#                     # t1cnt += 1
#                     t1_unions.append(ranks[0][i])
#             runfile = bed.BedFile(t1_unions, 'Limited')
#             pos = 0
#             neg = 0
#             if len(t1_unions) != 0:
#                 for ind, cls in enumerate(classified_peaks):
#                     val = overlap(cls, runfile)
#                     if ind == 0:
#                         tpr = int(val[1]) / len(cls)
#                         TPR.append(tpr)
#                     else:
#                         fpr = int(val[1]) / len(cls)
#                         FPR.append(fpr)
#             else:
#                 TPR.append(0)
#                 FPR.append(0)
#
#     else:
#         alpha = np.arange(1000, 0, -10)
#         TPR = []
#         FPR = []
#         for i in alpha:
#
#             idxs_t1 = []
#
#             # Calculate Maximum possible entropy value
#
#             for idx, val in enumerate(IDR):
#                 if val.score <= i:
#                     idxs_t1.append(idx)
#
#             t1_unions = []
#             for j in range(len(ranks[0])):
#                 if j in idxs_t1:
#                     # t1cnt += 1
#                     t1_unions.append(ranks[0][j])
#             runfile = bed.BedFile(t1_unions, 'Limited')
#             pos = 0
#             neg = 0
#             if len(t1_unions) != 0:
#                 for ind, cls in enumerate(classified_peaks):
#                     val = overlap(cls, runfile)
#                     if ind == 0:
#                         tpr = int(val[1]) / len(cls)
#                         TPR.append(tpr)
#                     else:
#                         fpr = int(val[1]) / len(cls)
#                         FPR.append(fpr)
#             else:
#                 TPR.append(0)
#                 FPR.append(0)
#
#     for ind, cls in enumerate(classified_peaks):
#         overs = []
#         for bedf in reps:
#             val = overlap(cls, bedf)
#             for v in val[0]:
#                 if v not in overs:
#                     overs.append(v)
#         if ind == 0:
#             tp = len(overs)
#             fn = len(cls) - tp
#             tpr = tp / (tp + fn)
#             TPR.append(tpr)
#         else:
#             fp = len(overs)
#             tn = len(cls) - fp
#             spec = tn / (tn + fp)
#             fpr = 1 - spec
#             FPR.append(fpr)
#
#     return TPR, FPR, ent
#
#

#
# def scrubBED(bedfile):
#     pos = []
#     neg = []
#     amb = []
#     for i in bedfile:
#         try:
#             if i.peak == 1:
#                 pos.append(i)
#             elif i.peak == -1:
#                 neg.append(i)
#         except AttributeError:
#             amb.append(i)
#     return pos, neg, amb
#
# MAX_CONT_snyderIDR = contTable(MAX_POS, maxlist, IDR = max_snyder_idr_run, bypeak=True, minentries=4 , Binom=False)
# df = pd.DataFrame(MAX_CONT_snyderIDR)
# df.to_csv('~/Documents/Honours/MAX_snyderIDR_BP.csv')
#
# MAX_CONT_myersIDR_NEG = contTable(MAX_NEG, maxlist, IDR = max_myers_idr_run, bypeak=True, minentries=4 , Binom=False)
# df = pd.DataFrame(MAX_CONT_myersIDR_NEG)
# df.to_csv('~/Documents/Honours/MAX_myersIDR_NEG_BP.csv')
#

#
# """
# ChIP-seq Analysis
# """
#
# """
# MNT Analysis
# """
# mnt1 = bed.BedFile('ChIPdata/MNT_K562/Snyder_ABENCAB887GAG/MNT_K562_Snyder2_1.bed', 'Peaks')
# mnt2 = bed.BedFile('ChIPdata/MNT_K562/Snyder_ABENCAB887GAG/MNT_K562_Snyder2_2.bed', 'Peaks')
#
# mnt3 = bed.BedFile('ChIPdata/MNT_K562/Snyder_ABENCAB541WYB/MNT_K562_Snyder1_1.bed', 'Peaks')
# mnt4 = bed.BedFile('ChIPdata/MNT_K562/Snyder_ABENCAB541WYB/MNT_K562_Snyder1_2.bed', 'Peaks')
# MNT = [mnt1, mnt2, mnt3, mnt4]
#
#
# mnt_rp = RankProd.performrankprod(MNT, 1, filename='MNT_min2.bed')
#
#
# """
# REST/NRSF Analysis
# """
#
# #REST Michael Snyder k562 hg19
#
#
rest_snyder_idr_opt = bed.BedFile(RankProd.collapse(bed.BedFile('ChIPdata/NRSF_K562/Snyder_hg19/IDR/IDR_final_optimal.narrowPeak', 'Peaks'), idr=True))
# # rest_snyder_idr_con = bed.BedFile('ChIPdata/NRSF_K562/Snyder_hg19/IDR/IDR_final_conservative.narrowPeak', 'Peaks')
# # rest_snyder_idr_run = bed.BedFile(RankProd.collapse(bed.BedFile('ChIPdata/NRSF_K562/Snyder_hg19/idrValues.txt', 'Peaks'), idr=True))
#
# # REST_so_idr = overlap(LO_REST_K562_Classified, rest_snyder_idr_opt)
# # REST_sc_idr = overlap(LO_REST_K562_Classified, rest_snyder_idr_con)
# #
# # REST_so_idr_POS = overlap(REST_POS, rest_snyder_idr_opt)
# # REST_so_idr_NEG = overlap(REST_NEG, rest_snyder_idr_opt)
# # REST_so_idr_AMB = overlap(REST_AMB, rest_snyder_idr_opt)
# #
# # REST_sc_idr_POS = overlap(REST_POS, rest_snyder_idr_con)
# # REST_sc_idr_NEG = overlap(REST_NEG, rest_snyder_idr_con)
# # REST_sc_idr_AMB = overlap(REST_AMB, rest_snyder_idr_con)
# #
# # REST_sr_idr_POS = overlap(REST_POS, rest_snyder_idr_run)
# # REST_sr_idr_NEG = overlap(REST_NEG, rest_snyder_idr_run)
# # REST_sr_idr_AMB = overlap(REST_AMB, rest_snyder_idr_run)
#
# #REST Richard Myers k562 hg19
#
#
# rest_myers_idr_opt = bed.BedFile(RankProd.collapse(bed.BedFile('ChIPdata/NRSF_K562/Myers_hg19/IDR/optimal.bed', 'Peaks'), idr=True))
# rest_myers_idr_con = bed.BedFile('ChIPdata/NRSF_K562/Myers_hg19/IDR/conservative.bed', 'Peaks')
# rest_myers_idr_run = bed.BedFile(RankProd.collapse(bed.BedFile('ChIPdata/NRSF_K562/Myers_hg19/idrValues.txt', 'Peaks'), idr=True))
#
# REST_POS = bed.BedFile('ChIPdata/NRSF_K562/REST_POS_hg19.bed', 'Limited')
# REST_NEG = bed.BedFile('ChIPdata/NRSF_K562/REST_NEG_hg19.bed', 'Limited')
# # REST_AMB = bed.BedFile('ChIPdata/NRSF_K562/REST_AMB_hg19.bed', 'Limited')
# #
# # rest_s1m1_idr = bed.BedFile('ChIPdata/NRSF_K562/NRSFM1_S1.narrowPeak', 'Peaks')
# # rest_s1m2_idr = bed.BedFile('ChIPdata/NRSF_K562/NRSFM2_S1.narrowPeak', 'Peaks')
# # rest_s2m1_idr = bed.BedFile('ChIPdata/NRSF_K562/NRSFM1_S2.narrowPeak', 'Peaks')
# # rest_s2m2_idr = bed.BedFile('ChIPdata/NRSF_K562/NRSFM2_S2.narrowPeak', 'Peaks')
# #
# #
# # bed.writeBedFile(REST_K562_Classified, 'hg18_limited_rest_classified_peaks.bed', 'Limited')
# #
# #
# def getTopIDR(idr, sc=0.05):
#     top = []
#     for i in idr:
#         if i.globalIDR <= sc:
#             if i.globalIDR == 5.0:
#                 i.pValue = 0.0000000000000000000000000001
#             else:
#                 i.pValue = 10**-float(i.globalIDR)
#                 top.append(i)
#     f = bed.BedFile(top, 'Peaks')
#     return f
# #
# # rest_s1m1_idr_540 = getTopIDR(rest_s1m1_idr)
# # rest_s1m2_idr_540 = getTopIDR(rest_s1m2_idr)
# # rest_s2m1_idr_540 = getTopIDR(rest_s2m1_idr)
# # rest_s2m2_idr_540 = getTopIDR(rest_s2m2_idr)
# #
# rest_snyder_idr_540 = getTopIDR(rest_snyder_idr_run)
# rest_myers_idr_540 = getTopIDR(rest_myers_idr_run)
# restOverlap1 = overlap(rest_snyder_idr_540, rest_myers_idr_540)[1:3]
# restOverlap2 = overlap(rest_myers_idr_540, rest_snyder_idr_540)[1:3]
#
#
# rest_idr = []
# for i in rest_snyder_idr_540:
#     rest_idr.append(i)
# for i in rest_myers_idr_540:
#     rest_idr.append(i)
# rest_idr = bed.BedFile(rest_idr, 'Peaks')
#
# # rest_idrOverlap = overlap(rest_snyder_idr_run, rest_myers_idr_run)
# #
# #
# results_restSny = RankProd.performrankprod(rest,rankmethod='pvalue', minentries=3, filename='REST_hg19_RP_snyder.bed')
# # ALL_restSny = bed.BedFile('ChIPdata/NRSF_K562/ALL_REST_hg19_RP_snyder.bed', 'Peaks')
# T2_restSny = bed.BedFile('T2_REST_hg19_RP_snyder.bed', 'RP')
# T1_restSny = bed.BedFile('T1_REST_hg19_RP_snyder.bed', 'RP')
# snyOverlap1 = overlap(rest_myers_idr_540, T2_restSny)[1:3]
# snyOverlap2 = overlap( rest_snyder_idr_540, T2_restSny)[1:3]
# snyOverlap = overlap(rest_idr, T2_restSny)[1:3]
#
# #
# #
# results_restMye = RankProd.performrankprod([rest3, rest4], minentries=2, filename='REST_hg19_RP_myers.bed')
# # ALL_restMye = bed.BedFile('ChIPdata/NRSF_K562/ALL_REST_hg19_RP_myers.bed', 'Peaks')
# # T2_restMye = bed.BedFile('T2_REST_hg19_RP_myers.bed', 'Peaks')
# # T1_restMye = bed.BedFile('ChIPdata/NRSF_K562/T1_REST_hg19_RP_myers.bed', 'Peaks')
# myeOverlap1 = overlap(rest_myers_idr_540, T2_restMye )[1:3]
# myeOverlap2 = overlap(rest_snyder_idr_540, T2_restMye )[1:3]
# myeOverlap = overlap(rest_idr, T2_restMye)[1:3]
#
#
# #RankProd.performrankprod(rest, minentries=2, filename='REST_hg19_RP_min2.bed')
# #RankProd.performrankprod(rest, minentries=3, filename='REST_hg19_RP_min3.bed')
# #RankProd.performrankprod(rest, minentries=4, filename='REST_hg19_RP_min4.bed')
# # #
# #
# # ALL_restMin4 = bed.BedFile('ChIPdata/NRSF_K562/ALL_REST_hg19_RP_min4.bed', 'Peaks')
# # T2_restMin4 = bed.BedFile('T2_REST_hg19_RP_min4.bed', 'Peaks')
# # T1_restMin4 = bed.BedFile('ChIPdata/NRSF_K562/T1_REST_hg19_RP_min4.bed', 'Peaks')
# restMin4 = overlap( rest_idr, T2_restMin4)[1:3]
# restMin41 = overlap(rest_myers_idr_540, T2_restMin4)[1:3]
# restMin42 = overlap(rest_snyder_idr_540, T2_restMin4)[1:3]
#
# #
# # ALL_restMin3 = bed.BedFile('ChIPdata/NRSF_K562/ALL_REST_hg19_RP_min3.bed', 'Peaks')
# # T2_restMin3 = bed.BedFile('T2_REST_hg19_RP_min3.bed', 'Peaks')
# # T1_restMin3 = bed.BedFile('ChIPdata/NRSF_K562/T1_REST_hg19_RP_min3.bed', 'Peaks')
# restMin3 = overlap(rest_idr, T2_restMin3)[1:3]
# restMin31 = overlap(rest_myers_idr_540, T2_restMin3)[1:3]
# restMin32 = overlap(rest_snyder_idr_540, T2_restMin3)[1:3]
#
# #
# #ALL_restMin2 = bed.BedFile('ALL_REST_hg19_RP_min2.bed', 'Peaks')
# # T2_restMin2 = bed.BedFile('T2_REST_hg19_RP_min2.bed', 'Peaks')
# # T1_restMin2 = bed.BedFile('ChIPdata/NRSF_K562/T1_REST_hg19_RP_min2.bed', 'Peaks')
# restMin2 = overlap(rest_idr, T2_restMin2)[1:3]
# restMin21 = overlap(rest_myers_idr_540, T2_restMin2)[1:3]
# restMin22 = overlap(rest_snyder_idr_540, T2_restMin2)[1:3]
#
#
# #
#
# """
# MAX analysis
# """
#
# #Max Myers K562 hg19
# max1 = bed.BedFile('ChIPdata/Max_K562/Myers_hg19/Max_K562_Myers_1.bed', 'Peaks')
# max2 = bed.BedFile('ChIPdata/Max_K562/Myers_hg19/Max_K562_Myers_2.bed', 'Peaks')
#
# # max_myers_idr_opt = bed.BedFile(RankProd.collapse(bed.BedFile('ChIPdata/Max_K562/Myers_hg19/IDR/optimal.bed', 'Peaks'), idr=True))
# #max_myers_idr_con = bed.BedFile('ChIPdata/Max_K562/Myers_hg19/IDR/conservative.bed', 'Peaks')
# # max_myers_idr_run = bed.BedFile(RankProd.collapse(bed.BedFile('ChIPdata/Max_K562/Myers_hg19/MAX_myers_idrValues.bed', 'Peaks'), idr=True))
#
#
# #Max Snyder K562 hg19
# max3 = bed.BedFile('ChIPdata/Max_K562/Snyder_hg19/Max_K562_Snyder_1.bed', 'Peaks')
# max4 = bed.BedFile('ChIPdata/Max_K562/Snyder_hg19/Max_K562_Snyder_2.narrowPeak', 'Peaks')
# maxlist = [max1, max2, max3, max4]
#
# #max_snyder_idr_opt = bed.BedFile(RankProd.collapse(bed.BedFile('ChIPdata/Max_K562/Snyder_hg19/IDR/optimal.bed', 'Peaks'), idr=True))
# # max_snyder_idr_con = bed.BedFile('ChIPdata/Max_K562/Snyder_hg19/IDR/conservative.bed', 'Peaks')
# max_snyder_idr_run = bed.BedFile(RankProd.collapse(bed.BedFile('ChIPdata/Max_K562/Snyder_hg19/MAX_snyder_idrValues.bed', 'Peaks'), idr=True))
#
# # max_idr_so = overlap(LO_MAX_K562_Classified, max_snyder_idr_opt)
# # max_idr_sc = overlap(LO_MAX_K562_Classified, max_snyder_idr_con)
#
# MAX_POS = bed.BedFile('ChIPdata/Max_K562/MAX_POS_hg19.bed', 'Limited')
# MAX_NEG = bed.BedFile('ChIPdata/Max_K562/MAX_NEG_hg19.bed', 'Limited')
# # MAX_AMB = bed.BedFile('ChIPdata/Max_K562/MAX_AMB_hg19.bed', 'Limited')
# #
# #
# # # results_max = RankProd.performrankprod(maxlist, minentries=2, alpha=0.05, filename=None)
# #
# max_snyder_idr_540 = getTopIDR(max_snyder_idr_run)
# max_myers_idr_540 = getTopIDR(max_myers_idr_run)
# maxidrOverlap1 = overlap(max_myers_idr_540, max_snyder_idr_540)[1:3]
# maxidrOverlap2 = overlap(max_snyder_idr_540, max_myers_idr_540)[1:3]
#
# max_idr = []
# for i in max_snyder_idr_540:
#     max_idr.append(i)
# for i in max_myers_idr_540:
#     max_idr.append(i)
# max_idr = bed.BedFile(max_idr, 'Peaks')
# # max_idrOverlap = overlap(max_myers_idr_540, max_snyder_idr_540)
# #
# #
# #results_maxSny = RankProd.performrankprod([max3, max4], minentries=2, filename='MAX_hg19_RP_snyder.bed')
# # ALL_maxSny = bed.BedFile('ChIPdata/Max_K562/ALL_MAX_hg19_RP_snyder.bed', 'Peaks')
# T2_maxSny = bed.BedFile('T2_MAX_hg19_RP_snyder.bed', 'Peaks')
# # T1_maxSny = bed.BedFile('ChIPdata/Max_K562/T1_MAX_hg19_RP_snyder.bed', 'Peaks')
# snyOverlap1 = overlap(max_myers_idr_540,T2_maxSny)[1:3]
# snyOverlap2 = overlap(max_snyder_idr_540,T2_maxSny)[1:3]
# snyOverlap = overlap(max_idr,T2_maxSny)[1:3]
#
#
# #
# #
# #results_maxMye = RankProd.performrankprod([max1, max2], minentries=2, filename='MAX_hg19_RP_myers.bed')
# # ALL_maxMye = bed.BedFile('ChIPdata/Max_K562/ALL_MAX_hg19_RP_myers.bed', 'Peaks')
# T2_maxMye = bed.BedFile('T2_MAX_hg19_RP_myers.bed', 'Peaks')
# # T1_maxMye = bed.BedFile('ChIPdata/Max_K562/T1_MAX_hg19_RP_myers.bed', 'Peaks')
# myeOverlap1 = overlap(max_myers_idr_540,T2_maxMye)[1:3]
# myeOverlap2 = overlap(max_snyder_idr_540,T2_maxMye)[1:3]
# myeOverlap = overlap(max_idr,T2_maxMye)[1:3]
#
# # maxRPmye2Overlap = overlap(T2_maxMye, max_snyder_idr_540)
# # maxRPMyeallIDR = overlap(T2_maxMye, max_snyder_idr_540, max_myers_idr_540)
# # maxRPOverlap = overlap(ALL_maxSny, ALL_maxMye)
# #
#
# #results_max = RankProd.performrankprod(maxlist, minentries=2, filename='MAX_hg19_RP_min2.bed')
# #results_max = RankProd.performrankprod(maxlist, minentries=3, filename='MAX_hg19_RP_min3.bed')
# #results_max = RankProd.performrankprod(maxlist, minentries=4, filename='MAX_hg19_RP_min4.bed')
# #
# #
# #
# #
# # ALL_maxMin4 = bed.BedFile('ChIPdata/Max_K562/ALL_MAX_hg19_RP_min4.bed', 'Peaks')
# #T2_maxMin4 = bed.BedFile('T2_MAX_hg19_RP_min4.bed', 'Peaks')
# # T1_maxMin4 = bed.BedFile('ChIPdata/Max_K562/T1_MAX_hg19_RP_min4.bed', 'Peaks')
# maxMin4 = overlap(max_idr, T2_maxMin4)[1:3]
# maxMin41 = overlap(max_myers_idr_540, T2_maxMin4)[1:3]
# maxMin42 = overlap(max_snyder_idr_540, T2_maxMin4)[1:3]
# # maxRPMyemin4IDR = overlap(T2_maxMin4, max_snyder_idr_540, max_myers_idr_540)
# #
# #
# # ALL_maxMin3 = bed.BedFile('ChIPdata/Max_K562/ALL_MAX_hg19_RP_min3.bed', 'Peaks')
# #T2_maxMin3 = bed.BedFile('T2_MAX_hg19_RP_min3.bed', 'Peaks')
# # T1_maxMin3 = bed.BedFile('ChIPdata/Max_K562/T1_MAX_hg19_RP_min3.bed', 'Peaks')
# maxMin3 = overlap( max_idr, T2_maxMin3)[1:3]
# maxMin31 = overlap(max_myers_idr_540, T2_maxMin3)[1:3]
# maxMin32 = overlap(max_snyder_idr_540, T2_maxMin3)[1:3]
#
# # maxRPMyemin3IDR = overlap(T2_maxMin3, max_snyder_idr_540, max_myers_idr_540)
# #
# #
# # ALL_maxMin2 = bed.BedFile('ChIPdata/Max_K562/ALL_MAX_hg19_RP_min2.bed', 'Peaks')
# #T2_maxMin2 = bed.BedFile('T2_MAX_hg19_RP_min2.bed', 'Peaks')
# # T1_maxMin2 = bed.BedFile('ChIPdata/Max_K562/T1_MAX_hg19_RP_min2.bed', 'Peaks')
# maxMin2 = overlap( max_idr, T2_maxMin2)[1:3]
# maxMin21 = overlap(max_myers_idr_540, T2_maxMin2)[1:3]
# maxMin22 = overlap(max_snyder_idr_540, T2_maxMin2)[1:3]
#
#
# # maxRPMyemin2IDR = overlap(T2_maxMin2, max_snyder_idr_540, max_myers_idr_540)
# #
# #
# #
#
#
# """
# SRF analysis
# """
#
# #SRF Myers GM12878 hg19
# # srf1 = bed.BedFile('ChIPdata/SRF_GM12878/Myers_hg19/SRF_GM12878_Myers_1.narrowPeak', 'Peaks')
# # srf2 = bed.BedFile('ChIPdata/SRF_GM12878/Myers_hg19/SRF_GM12878_Myers_2.narrowPeak', 'Peaks')
# #
# # # srf_myers_idr_opt = bed.BedFile(RankProd.collapse(bed.BedFile('ChIPdata/SRF_GM12878/Myers_hg19/IDR/optimal.bed', 'Peaks'), idr=True))
# # # srf_myers_idr_con = bed.BedFile('ChIPdata/SRF_GM12878/Myers_hg19/IDR/conservative.bed', 'Peaks')
# # srf_myers_idr_run = bed.BedFile(RankProd.collapse(bed.BedFile('ChIPdata/SRF_GM12878/Myers_hg19/idrValues.txt', 'Peaks'), idr=True))
#
# #
# #SRF Snyder GM12878 hg19
# srf3 = bed.BedFile('ChIPdata/SRF_GM12878/Snyder_hg19/SRF_GM12878_Snyder_1.regionPeak', 'Peaks')
# srf4 = bed.BedFile('ChIPdata/SRF_GM12878/Snyder_hg19/SRF_GM12878_Snyder_2.regionPeak', 'Peaks')
#
# #SRF Myers2 GM12878 hg19
# srf1 = bed.BedFile('ChIPdata/SRF_GM12878/Myers2_hg19/SRF_GM12878_Myers2_1.bed', 'Peaks')
# srf2 = bed.BedFile('ChIPdata/SRF_GM12878/Myers2_hg19/SRF_GM12878_Myers2_2.regionPeak', 'Peaks')
# # srf_myers_idr_run = bed.BedFile(RankProd.collapse(bed.BedFile('ChIPdata/SRF_GM12878/Myers2_hg19/SRF_GM12878_Myers2IDRValues.bed', 'Peaks'), idr=True))
#
#
# srf = [srf1, srf2, srf3, srf4]
#
#
#
# # srf_snyder_idr_opt = bed.BedFile(RankProd.collapse(bed.BedFile('ChIPdata/SRF_GM12878/Snyder_hg19/IDR/IDR_final_optimal.narrowPeak', 'Peaks'), idr=True))
# # srf_snyder_idr_con = bed.BedFile('ChIPdata/SRF_GM12878/Snyder_hg19/IDR/IDR_final_conservative.narrowPeak', 'Peaks')
# # srf_snyder_idr_run = bed.BedFile(RankProd.collapse(bed.BedFile('ChIPdata/SRF_GM12878/Snyder_hg19/idrValues.txt', 'Peaks'), idr=True))
#
#
#
# SRF_POS = bed.BedFile('ChIPdata/SRF_GM12878/SRF_POS_hg19.bed', 'Limited')
# SRF_NEG = bed.BedFile('ChIPdata/SRF_GM12878/SRF_NEG_hg19.bed', 'Limited')
# SRF_AMB = bed.BedFile('ChIPdata/SRF_GM12878/SRF_AMB_hg19.bed', 'Limited')
#
# # results_srf = RankProd.performrankprod([srf1, srf2], minentries=2, filename=None)
#
# srf_snyder_idr_540 = getTopIDR(srf_snyder_idr_run)
# srf_myers_idr_540 = getTopIDR(srf_myers_idr_run)
# srf_idrOverlap = overlap(srf_myers_idr_540, srf_snyder_idr_540)[1:3]
# srfidrOverlap2 = overlap(srf_snyder_idr_540, srf_myers_idr_540)[1:3]
#
# srf_idr = []
# for i in srf_snyder_idr_540:
#     srf_idr.append(i)
# for i in srf_myers_idr_540:
#     srf_idr.append(i)
# srf_idr = bed.BedFile(srf_idr, 'Peaks')
#
#
# # RankProd.performrankprod([srf1, srf2, srf3], minentries=2, filename='SRF_hg19_RP_3r_min2_myers1.bed')
# # RankProd.performrankprod([srf1, srf2, srf4], minentries=2, filename='SRF_hg19_RP_3r_min2_myers2.bed')
# # RankProd.performrankprod([srf3, srf4, srf1], minentries=2, filename='SRF_hg19_RP_3r_min2_snyder1.bed')
# # RankProd.performrankprod([srf3, srf4, srf2], minentries=2, filename='SRF_hg19_RP_3r_min2_snyder2.bed')
#
#
#
#
#
# RankProd.performrankprod([srf3, srf4], minentries=2, filename='SRF_hg19_RP_snyder.bed')
# # ALL_srfSny = bed.BedFile('ChIPdata/SRF_GM12878/ALL_SRF_hg19_RP_snyder.bed', 'Peaks')
# T2_srfSny = bed.BedFile('T2_SRF_hg19_RP_snyder.bed', 'Peaks')
# # T1_srfSny = bed.BedFile('ChIPdata/SRF_GM12878/T1_SRF_hg19_RP_snyder.bed', 'Peaks')
# srfsnyOverlap1 = overlap(srf_myers_idr_540, T2_srfSny)[1:3]
# srfsnyOverlap2 = overlap(srf_snyder_idr_540, T2_srfSny)[1:3]
# srfsnyOverlap = overlap(srf_idr, T2_srfSny)[1:3]
#
#
# RankProd.performrankprod([srf1, srf2], minentries=2, filename='SRF_hg19_RP_myers.bed')
# T2_srfMye = bed.BedFile('T2_SRF_hg19_RP_myers.bed', 'Peaks')
# srfmyeOverlap1 = overlap(srf_myers_idr_540, T2_srfMye)[1:3]
# srfmyeOverlap2 = overlap(srf_snyder_idr_540, T2_srfMye)[1:3]
# srfMyeOverlap = overlap(srf_idr, T2_srfMye)[1:3]
#
#
# RankProd.performrankprod(srf, minentries=2, filename='SRF_hg19_RP_min2.bed')
# RankProd.performrankprod(srf, minentries=3, filename='SRF_hg19_RP_min3.bed')
# RankProd.performrankprod(srf, minentries=4, filename='SRF_hg19_RP_min4.bed')
# # #
# #
# # T2_srfMin4 = bed.BedFile('T2_SRF_hg19_RP_min4.bed', 'Peaks')
# srfMin4 = overlap( srf_idr, T2_srfMin4)[1:3]
# srfMin41 = overlap(srf_myers_idr_540, T2_srfMin4)[1:3]
# srfMin42 = overlap(srf_snyder_idr_540, T2_srfMin4)[1:3]
#
# #
# # ALL_srfMin3 = bed.BedFile('ChIPdata/NRSF_K562/ALL_SRF_hg19_RP_min3.bed', 'Peaks')
# # T2_srfMin3 = bed.BedFile('T2_SRF_hg19_RP_min3.bed', 'Peaks')
# # T1_srfMin3 = bed.BedFile('ChIPdata/NRSF_K562/T1_SRF_hg19_RP_min3.bed', 'Peaks')
# srfMin3 = overlap( srf_idr, T2_srfMin3)[1:3]
# srfMin31 = overlap(srf_myers_idr_540, T2_srfMin3)[1:3]
# srfMin32 = overlap(srf_snyder_idr_540, T2_srfMin3)[1:3]
#
# #
# # ALL_srfMin2 = bed.BedFile('ChIPdata/NRSF_K562/ALL_SRF_hg19_RP_min2.bed', 'Peaks')
# # T2_srfMin2 = bed.BedFile('T2_SRF_hg19_RP_min2.bed', 'Peaks')
# # T1_srfMin2 = bed.BedFile('ChIPdata/NRSF_K562/T1_SRF_hg19_RP_min2.bed', 'Peaks')
# srfMin2 = overlap(srf_idr, T2_srfMin2)[1:3]
# srfMin21 = overlap(srf_myers_idr_540, T2_srfMin2)[1:3]
# srfMin22 = overlap(srf_snyder_idr_540, T2_srfMin2)[1:3]
#
#
# max_ranks = RankProd.performrankprod(maxlist, minentries=2, filename=None)
# maxPks = thresholdCalc(max_ranks[1], np.arange(0, 1, 0.01))
#
# rest_ranks = RankProd.performrankprod(rest, minentries=4, filename=None)
# restPs = [1-p for p in rest_ranks[1]]
# rPs = thresholdCalc(restPs)
#
# rr = retrieveRanks(rest_ranks[0][0])
# rrnorm = []
# for r in rr:
#     rrnorm.append([x/len(r) for x in r])
#
# # restmad = wilcoxonRankSign(rr[0], rr[1])
# restmad = MAD(zip(*rrnorm), n = 2)
# rankmds = scipy.stats.rankdata(restmad[1])
# mdsprod = [x*y for x,y in zip(rest_ranks[0][1], rankmds)]
# mdsdrpb = RankProd.rankprodbounds(mdsprod, len(mdsprod), 5, 'geometric')
# combps = [scipy.stats.combine_pvalues([x, y])[1] for x, y in zip(rest_ranks[1], mdsdrpb)]
# log10rest = [-math.log10(i) for i in rest_ranks[1]]
# log10MAD = [-math.log10(i) for i in restmad[1]]
# # log10SAD = [math.log10(i) for i in restmad[0]]
# plt.plot(log10MAD, log10rest, '.')
# plt.xlabel("Average aboslute deviation from the mean of the ranks")
# plt.ylabel("Log10(p-values)")
# madmean = np.mean(log10MAD)
# sumlogs = [x+y for x, y in zip(log10rest, log10MAD)]
# difflogs = [x-y for x, y in zip(log10rest, log10MAD)]
# limitedsumlogs = [-math.log10(i) for i in sumlogs]
# sumunlogged = [10**-i for i in sumlogs]
# diffunlogged = [10**-i for i in difflogs]
# mult = [x*y for x, y in zip(log10MAD, log10rest)]
# unlogmult = [10**-x for x in mult]
#
# """t1 rest"""
# t1ps = []
# t1w = []
# t1mads = []
# t1sads = []
# otsads = []
# indexes = []
# for i, v in enumerate(rest_ranks[1]):
#     if v <= 0.05:
#         t1ps.append(v)
#         t1w.append(restmad[0][i])
#         t1mads.append(restmad[1][i])
#         t1sads.append(restmad[0][i])
#         indexes.append(i)
#     else:
#         otsads.append(restmad[0][i])
#
# t1log10rest = [-math.log10(i) for i in t1ps]
# t1log10MAD = [-math.log10(i) for i in t1mads]
# t1sumlogs = [x+y for x, y in zip(t1log10rest, t1log10MAD)]
# t1difflogs = [x-y for x, y in zip(t1log10rest, t1log10MAD)]
# t1sumun = [10**-i for i in t1sumlogs]
# t1diffun = [10**-i for i in t1difflogs]
# # sumt1unlog =
#
# plt.plot(t1log10MAD, t1log10rest, '.')
# plt.show()
#
#
# srf_ranks = RankProd.performrankprod([srf1, srf2], minentries=2, filename=None)
# sr = retrieveRanks(srf_ranks[0][0])
# srnorm = []
# for r in sr:
#     srnorm.append([x/len(r) for x in r])
#
# # restmad = wilcoxonRankSign(rr[0], rr[1])
# srfmad = MAD(zip(*srnorm), n = 4)
# srfmds = scipy.stats.rankdata(srfmad[1])
# srfprod = [x*y for x,y in zip(srf_ranks[0][1], srfmds)]
# srfrpb = RankProd.rankprodbounds(srfprod, len(srfprod), 5, 'geometric')
# madcomb = [scipy.stats.combine_pvalues([x, y])[1] for x, y in zip(srf_ranks[1], srfmadrpb)]
#
# log10srf = [-math.log10(i) for i in srf_ranks[1]]
# log10MADsrf = [math.log10(i) for i in srfmad[1]]
# plt.plot(log10MADsrf, log10srf, '.')
# plt.show()
#
# """t1 srf"""
# t1ps = []
# t1w = []
# t1mads = []
# for i, v in enumerate(srf_ranks[1]):
#     if v <= 0.05:
#         t1ps.append(v)
#         t1w.append(srfmad[0][i])
#         t1mads.append(srfmad[1][i])
#
# t1log10srf = [-math.log10(i) for i in t1ps]
# t1log10MADsrf = [math.log10(i) for i in t1mads]
# plt.plot(t1log10MADsrf, t1log10srf, '.')
# plt.show()
#
# log10srf = [-math.log10(i) for i in srf_ranks[1]]
#
#
# alphas = [0.01, 0.05, 0.1, 0.15, 0.2, 0.5, 1]
#
# #
# # REST_CONT_min2Binom = contTable(REST_POS, rest, minentries=2 )
# # df = pd.DataFrame(REST_CONT_min2Binom)
# # df.to_csv('~/Documents/Honours/REST_RP_min2.csv')
# # REST_CONT_min3Binom = contTable(REST_POS, rest, minentries=3 )
# # df = pd.DataFrame(REST_CONT_min3Binom)
# # df.to_csv('~/Documents/Honours/REST_RP_min3.csv')
# # REST_CONT_min4Binom = contTable(REST_POS, rest, minentries=4 )
# # df = pd.DataFrame(REST_CONT_min4Binom)
# # df.to_csv('~/Documents/Honours/REST_RP_min4.csv')
# #
# REST_CONT_RP_snyder = contTable(REST_POS, [rest1, rest2], minentries=2 )
# # df = pd.DataFrame(REST_CONT_RP_snyder)
# # df.to_csv('~/Documents/Honours/REST_snyderRP.csv')
# # REST_CONT_RP_myers = contTable(REST_POS, [rest3, rest4], minentries=2 )
# # df = pd.DataFrame(REST_CONT_RP_myers)
# # df.to_csv('~/Documents/Honours/REST_myersRP.csv')
# #
# # REST_CONT_RP_snyder1_3r = contTable(REST_POS, [rest1, rest2, rest3], minentries=3 )
# # df = pd.DataFrame(REST_CONT_RP_snyder1_3r)
# # df.to_csv('~/Documents/Honours/REST_snyder1RP_3r.csv')
# # REST_CONT_RP_snyder2_3r = contTable(REST_POS, [rest1, rest2, rest4], minentries=3 )
# # df = pd.DataFrame(REST_CONT_RP_snyder2_3r)
# # df.to_csv('~/Documents/Honours/REST_snyder2RP_3r.csv')
# # REST_CONT_RP_myers1_3r = contTable(REST_POS, [rest3, rest4, rest1], minentries=3 )
# # df = pd.DataFrame(REST_CONT_RP_myers1_3r)
# # df.to_csv('~/Documents/Honours/REST_myers1RP_3r.csv')
# # REST_CONT_RP_myers2_3r = contTable(REST_POS, [rest3, rest4, rest2], minentries=3 )
# # df = pd.DataFrame(REST_CONT_RP_myers2_3r)
# # df.to_csv('~/Documents/Honours/REST_myers2RP_3r.csv')
# #
# # REST_CONT_RP_min2_snyder1_3r = contTable(REST_POS, [rest1, rest2, rest3], minentries=2 )
# # df = pd.DataFrame(REST_CONT_RP_min2_snyder1_3r)
# # df.to_csv('~/Documents/Honours/REST_snyder1RP_3r_min2.csv')
# # REST_CONT_RP_min2_snyder2_3r = contTable(REST_POS, [rest1, rest2, rest4], minentries=2 )
# # df = pd.DataFrame(REST_CONT_RP_min2_snyder2_3r)
# # df.to_csv('~/Documents/Honours/REST_snyder2RP_3r_min2.csv')
# # REST_CONT_RP_min2_myers1_3r = contTable(REST_POS, [rest3, rest4, rest1], minentries=2 )
# # df = pd.DataFrame(REST_CONT_RP_min2_myers1_3r)
# # df.to_csv('~/Documents/Honours/REST_myers1RP_3r_min2.csv')
# # REST_CONT_RP_min2_myers2_3r = contTable(REST_POS, [rest3, rest4, rest2], minentries=2 )
# # df = pd.DataFrame(REST_CONT_RP_min2_myers2_3r)
# # df.to_csv('~/Documents/Honours/REST_myers2RP_3r_min2.csv')
# #
# #
# # REST_CONT_min2Binom = contTable(REST_NEG, rest, minentries=2 )
# # df = pd.DataFrame(REST_CONT_min2Binom)
# # df.to_csv('~/Documents/Honours/REST_RP_min2_NEG.csv')
# # REST_CONT_min3Binom = contTable(REST_NEG, rest, minentries=3 )
# # df = pd.DataFrame(REST_CONT_min3Binom)
# # df.to_csv('~/Documents/Honours/REST_RP_min3_NEG.csv')
# # REST_CONT_min4Binom = contTable(REST_NEG, rest, minentries=4 )
# # df = pd.DataFrame(REST_CONT_min4Binom)
# # df.to_csv('~/Documents/Honours/REST_RP_min4_NEG.csv')
# #
# # REST_CONT_RP_snyder = contTable(REST_NEG, [rest1, rest2], minentries=2 )
# # df = pd.DataFrame(REST_CONT_RP_snyder)
# # df.to_csv('~/Documents/Honours/REST_snyderRP_NEG.csv')
# # REST_CONT_RP_myers = contTable(REST_NEG, [rest3, rest4], minentries=2 )
# # df = pd.DataFrame(REST_CONT_RP_myers)
# # df.to_csv('~/Documents/Honours/REST_myersRP_NEG.csv')
# #
# # REST_CONT_RP_snyder1_3r = contTable(REST_NEG, [rest1, rest2, rest3], minentries=3 )
# # df = pd.DataFrame(REST_CONT_RP_snyder1_3r)
# # df.to_csv('~/Documents/Honours/REST_snyder1RP_3r_NEG.csv')
# # REST_CONT_RP_snyder2_3r = contTable(REST_NEG, [rest1, rest2, rest4], minentries=3 )
# # df = pd.DataFrame(REST_CONT_RP_snyder2_3r)
# # df.to_csv('~/Documents/Honours/REST_snyder2RP_3r_NEG.csv')
# # REST_CONT_RP_myers1_3r = contTable(REST_NEG, [rest3, rest4, rest1], minentries=3 )
# # df = pd.DataFrame(REST_CONT_RP_myers1_3r)
# # df.to_csv('~/Documents/Honours/REST_myers1RP_3r_NEG.csv')
# # REST_CONT_RP_myers2_3r = contTable(REST_NEG, [rest3, rest4, rest2], minentries=3 )
# # df = pd.DataFrame(REST_CONT_RP_myers2_3r)
# # df.to_csv('~/Documents/Honours/REST_myers2RP_3r_NEG.csv')
# #
# # REST_CONT_RP_min2_snyder1_3r = contTable(REST_NEG, [rest1, rest2, rest3], minentries=2 )
# # df = pd.DataFrame(REST_CONT_RP_min2_snyder1_3r)
# # df.to_csv('~/Documents/Honours/REST_snyder1RP_3r_min2_NEG.csv')
# # REST_CONT_RP_min2_snyder2_3r = contTable(REST_NEG, [rest1, rest2, rest4], minentries=2 )
# # df = pd.DataFrame(REST_CONT_RP_min2_snyder2_3r)
# # df.to_csv('~/Documents/Honours/REST_snyder2RP_3r_min2_NEG.csv')
# # REST_CONT_RP_min2_myers1_3r = contTable(REST_NEG, [rest3, rest4, rest1], minentries=2 )
# # df = pd.DataFrame(REST_CONT_RP_min2_myers1_3r)
# # df.to_csv('~/Documents/Honours/REST_myers1RP_3r_min2_NEG.csv')
# # REST_CONT_RP_min2_myers2_3r = contTable(REST_NEG, [rest3, rest4, rest2], minentries=2 )
# # df = pd.DataFrame(REST_CONT_RP_min2_myers2_3r)
# # df.to_csv('~/Documents/Honours/REST_myers2RP_3r_min2_NEG.csv')
#
#
# REST_CONT_snyderIDR = contTable(REST_POS, rest, IDR = rest_snyder_idr_run, minentries=4 )
# df = pd.DataFrame(REST_CONT_snyderIDR)
# df.to_csv('~/Documents/Honours/REST_snyderIDR.csv')
# REST_CONT_myersIDR = contTable(REST_POS, rest, IDR = rest_myers_idr_run, minentries=4 )
# df = pd.DataFrame(REST_CONT_myersIDR)
# df.to_csv('~/Documents/Honours/REST_myersIDR.csv')
#
# REST_CONT_snyderIDR = contTable(REST_NEG, rest, IDR = rest_snyder_idr_run, minentries=4 )
# df = pd.DataFrame(REST_CONT_snyderIDR)
# df.to_csv('~/Documents/Honours/REST_snyderIDR_NEG.csv')
# REST_CONT_myersIDR = contTable(REST_NEG, rest, IDR = rest_myers_idr_run, minentries=4 )
# df = pd.DataFrame(REST_CONT_myersIDR)
# df.to_csv('~/Documents/Honours/REST_myersIDR_NEG.csv')
# # scipy.stats.chi2_contingency([REST_CONT_snyderIDR[0][0], REST_CONT_myersIDR[0][0]])
# #
# # REST_ROC = ROC([REST_POS, REST_NEG], rest, minentries=2, alpha=np.arange(0, 0.99, 0.01))
# # REST_ROC_snyder = ROC([REST_POS, REST_NEG], [rest1, rest2], minentries=2, alpha=np.arange(0.01, 0.5, 0.01))
# #
# # MAX_CONT_RP_min2Binom = contTable(MAX_POS, maxlist, minentries=2 )
# # df = pd.DataFrame(MAX_CONT_RP_min2Binom)
# # df.to_csv('~/Documents/Honours/MAX_RP_min2.csv')
# #
# #
# # MAX_CONT_RP_min3Binom = contTable(MAX_POS, maxlist, minentries=3 )
# # df = pd.DataFrame(MAX_CONT_RP_min3Binom)
# # df.to_csv('~/Documents/Honours/MAX_RP_min3.csv')
# # MAX_CONT_RP_min4Binom = contTable(MAX_POS, maxlist, minentries=4 )
# # df = pd.DataFrame(MAX_CONT_RP_min4Binom)
# # df.to_csv('~/Documents/Honours/MAX_RP_min4.csv')
# #
# #
# # MAX_CONT_RP_myers = contTable(MAX_POS, [max1, max2], minentries=2 )
# # df = pd.DataFrame(MAX_CONT_RP_myers)
# # df.to_csv('~/Documents/Honours/MAX_myersRP.csv')
# # MAX_CONT_RP_snyder = contTable(MAX_POS, [max3, max4], minentries=2 )
# # df = pd.DataFrame(MAX_CONT_RP_snyder)
# # df.to_csv('~/Documents/Honours/MAX_snyderRP.csv')
# #
# # MAX_CONT_RP_3r_myers1 = contTable(MAX_POS, [max1, max2, max3], minentries=3 )
# # df = pd.DataFrame(MAX_CONT_RP_3r_myers1)
# # df.to_csv('~/Documents/Honours/MAX_myers1RP_3r.csv')
# # MAX_CONT_RP_3r_myers2 = contTable(MAX_POS, [max1, max2, max4], minentries=3 )
# # df = pd.DataFrame(MAX_CONT_RP_3r_myers2)
# # df.to_csv('~/Documents/Honours/MAX_myers2RP_3r.csv')
# # MAX_CONT_RP_3r_snyder1 = contTable(MAX_POS, [max3, max4, max1], minentries=3 )
# # df = pd.DataFrame(MAX_CONT_RP_3r_snyder1)
# # df.to_csv('~/Documents/Honours/MAX_snyder1RP_3r.csv')
# # MAX_CONT_RP_3r_snyder2 = contTable(MAX_POS, [max3, max4, max2], minentries=3 )
# # df = pd.DataFrame(MAX_CONT_RP_3r_snyder2)
# # df.to_csv('~/Documents/Honours/MAX_snyder2RP_3r.csv')
# #
# MAX_CONT_RP_min2_3r_myers1 = contTable(MAX_POS, [max1, max2, max3], minentries=2 )
# df = pd.DataFrame(MAX_CONT_RP_min2_3r_myers1)
# df.to_csv('~/Documents/Honours/MAX_myers1RP_3r_min2.csv')
# # MAX_CONT_RP_min2_3r_myers2 = contTable(MAX_POS, [max1, max2, max4], minentries=2 )
# # df = pd.DataFrame(MAX_CONT_RP_min2_3r_myers2)
# # df.to_csv('~/Documents/Honours/MAX_myers2RP_3r_min2.csv')
# MAX_CONT_RP_min2_3r_snyder1 = contTable(MAX_POS, [max3, max4, max1], minentries=2 )
# df = pd.DataFrame(MAX_CONT_RP_min2_3r_snyder1)
# df.to_csv('~/Documents/Honours/MAX_snyder1RP_3r_min2.csv')
# MAX_CONT_RP_min2_3r_snyder2 = contTable(MAX_POS, [max3, max4, max2], minentries=2 )
# df = pd.DataFrame(MAX_CONT_RP_min2_3r_snyder2)
# df.to_csv('~/Documents/Honours/MAX_snyder2RP_3r_min2.csv')
# #
# #
# # ##MAX NEGATIVES
# # MAX_CONT_RP_min2_NEG = contTable(MAX_NEG, maxlist, minentries=2 )
# # df = pd.DataFrame(MAX_CONT_RP_min2_NEG)
# # df.to_csv('~/Documents/Honours/MAX_RP_min2_NEG.csv')
# #
# # MAX_CONT_RP_min3_NEG = contTable(MAX_NEG, maxlist, minentries=3 )
# # df = pd.DataFrame(MAX_CONT_RP_min3_NEG)
# # df.to_csv('~/Documents/Honours/MAX_RP_min3_NEG.csv')
# #
# # MAX_CONT_RP_min4_NEG = contTable(MAX_NEG, maxlist, minentries=4 )
# # df = pd.DataFrame(MAX_CONT_RP_min4_NEG)
# # df.to_csv('~/Documents/Honours/MAX_RP_min4_NEG.csv')
# #
# # MAX_CONT_RP_myers_NEG = contTable(MAX_NEG, [max1, max2], minentries=2 )
# # df = pd.DataFrame(MAX_CONT_RP_myers_NEG)
# # df.to_csv('~/Documents/Honours/MAX_myersRP_NEG.csv')
# # MAX_CONT_RP_snyder = contTable(MAX_NEG, [max3, max4], minentries=2 )
# # df = pd.DataFrame(MAX_CONT_RP_snyder)
# # df.to_csv('~/Documents/Honours/MAX_snyderRP_NEG.csv')
# #
# # MAX_CONT_RP_3r_myers1 = contTable(MAX_NEG, [max1, max2, max3], minentries=3 )
# # df = pd.DataFrame(MAX_CONT_RP_3r_myers1)
# # df.to_csv('~/Documents/Honours/MAX_myers1RP_3r_NEG.csv')
# # MAX_CONT_RP_3r_myers2 = contTable(MAX_NEG, [max1, max2, max4], minentries=3 )
# # df = pd.DataFrame(MAX_CONT_RP_3r_myers2)
# # df.to_csv('~/Documents/Honours/MAX_myers2RP_3r_NEG.csv')
# # MAX_CONT_RP_3r_snyder1 = contTable(MAX_NEG, [max3, max4, max1], minentries=3 )
# # df = pd.DataFrame(MAX_CONT_RP_3r_snyder1)
# # df.to_csv('~/Documents/Honours/MAX_snyder1RP_3r_NEG.csv')
# # MAX_CONT_RP_3r_snyder2 = contTable(MAX_NEG, [max3, max4, max2], minentries=3 )
# # df = pd.DataFrame(MAX_CONT_RP_3r_snyder2)
# # df.to_csv('~/Documents/Honours/MAX_snyder2RP_3r_NEG.csv')
# #
# MAX_CONT_RP_min2_3r_myers1 = contTable(MAX_NEG, [max1, max2, max3], minentries=2 )
# df = pd.DataFrame(MAX_CONT_RP_min2_3r_myers1)
# df.to_csv('~/Documents/Honours/MAX_myers1RP_3r_min2_NEG.csv')
# MAX_CONT_RP_min2_3r_myers2 = contTable(MAX_NEG, [max1, max2, max4], minentries=2 )
# df = pd.DataFrame(MAX_CONT_RP_min2_3r_myers2)
# df.to_csv('~/Documents/Honours/MAX_myers2RP_3r_min2_NEG.csv')
# MAX_CONT_RP_min2_3r_snyder1 = contTable(MAX_NEG, [max3, max4, max1], minentries=2 )
# df = pd.DataFrame(MAX_CONT_RP_min2_3r_snyder1)
# df.to_csv('~/Documents/Honours/MAX_snyder1RP_3r_min2_NEG.csv')
# MAX_CONT_RP_min2_3r_snyder2 = contTable(MAX_NEG, [max3, max4, max2], minentries=2 )
# df = pd.DataFrame(MAX_CONT_RP_min2_3r_snyder2)
# df.to_csv('~/Documents/Honours/MAX_snyder2RP_3r_min2_NEG.csv')
#
#
# # max_s1m1_idr_run = bed.BedFile(RankProd.collapse(bed.BedFile('ChIPdata/Max_K562/Snyder_hg19/MAX_s1m1_idrValues.bed', 'Peaks'), idr=True))
# # max_s1m2_idr_run = bed.BedFile(RankProd.collapse(bed.BedFile('ChIPdata/Max_K562/Snyder_hg19/MAX_s1m2_idrValues.bed', 'Peaks'), idr=True))
# # max_s2m1_idr_run = bed.BedFile(RankProd.collapse(bed.BedFile('ChIPdata/Max_K562/Snyder_hg19/MAX_s2m1_idrValues.bed', 'Peaks'), idr=True))
# # max_s2m2_idr_run = bed.BedFile(RankProd.collapse(bed.BedFile('ChIPdata/Max_K562/Snyder_hg19/MAX_s2m2_idrValues.bed', 'Peaks'), idr=True))
#
#
# MAX_CONT_myersIDR = contTable(MAX_POS, maxlist, IDR = max_myers_idr_run, minentries=4 , Binom=False)
# df = pd.DataFrame(MAX_CONT_myersIDR)
# df.to_csv('~/Documents/Honours/MAX_myersIDR.csv')
#
#
# MAX_CONT_snyderIDR = contTable(MAX_POS, maxlist, IDR = max_snyder_idr_run, minentries=4 , Binom=False)
# df = pd.DataFrame(MAX_CONT_snyderIDR)
# df.to_csv('~/Documents/Honours/MAX_snyderIDR.csv')
#
# MAX_CONT_myersIDR_NEG = contTable(MAX_NEG, maxlist, IDR = max_myers_idr_run, minentries=4 , Binom=False)
# df = pd.DataFrame(MAX_CONT_myersIDR_NEG)
# df.to_csv('~/Documents/Honours/MAX_myersIDR_NEG.csv')
#
#
# MAX_CONT_snyderIDR_NEG = contTable(MAX_NEG, maxlist, IDR = max_snyder_idr_run, minentries=4 , Binom=False)
# df = pd.DataFrame(MAX_CONT_snyderIDR_NEG)
# df.to_csv('~/Documents/Honours/MAX_snyderIDR_NEG.csv')
#
# # MAX_CONT_s1m1IDR = contTable(MAX_POS, maxlist, IDR = max_s1m1_idr_run, minentries=4 , Binom=False)
# # df = pd.DataFrame(MAX_CONT_s1m1IDR)
# # df.to_csv('~/Documents/Honours/MAX_s1m1IDR.csv')
# # MAX_CONT_s1m2IDR = contTable(MAX_POS, maxlist, IDR = max_s1m2_idr_run, minentries=4 , Binom=False)
# # df = pd.DataFrame(MAX_CONT_s1m2IDR)
# # df.to_csv('~/Documents/Honours/MAX_s1m2IDR.csv')
# # MAX_CONT_s2m1IDR = contTable(MAX_POS, maxlist, IDR = max_s2m1_idr_run, minentries=4 , Binom=False)
# # df = pd.DataFrame(MAX_CONT_s2m1IDR)
# # df.to_csv('~/Documents/Honours/MAX_s2m1IDR.csv')
# # MAX_CONT_s2m2IDR = contTable(MAX_POS, maxlist, IDR = max_s2m2_idr_run, minentries=4 , Binom=False)
# # df = pd.DataFrame(MAX_CONT_s2m2IDR)
# # df.to_csv('~/Documents/Honours/MAX_s2m2IDR.csv')
#
#
# #
# # SRF_CONT_RP_min2 = contTable(SRF_POS, srf, minentries=2 )
# # df = pd.DataFrame(SRF_CONT_RP_min2)
# # df.to_csv('~/Documents/Honours/SRF_RP_min2.csv')
# # SRF_CONT_RP_min3 = contTable(SRF_POS, srf, minentries=3 )
# # df = pd.DataFrame(SRF_CONT_RP_min3)
# # df.to_csv('~/Documents/Honours/SRF_RP_min3.csv')
# # SRF_CONT_RP_min4 = contTable(SRF_POS, srf, minentries=4 )
# # df = pd.DataFrame(SRF_CONT_RP_min4)
# # df.to_csv('~/Documents/Honours/SRF_RP_min4.csv')
# #
# # SRF_CONT_RP_3R_myers1 = contTable(SRF_POS, [srf1, srf2, srf3], minentries=3 )
# # df = pd.DataFrame(SRF_CONT_RP_3R_myers1)
# # df.to_csv('~/Documents/Honours/SRF_myers1RP_3r.csv')
# # SRF_CONT_RP_3R_myers2 = contTable(SRF_POS, [srf1, srf2, srf4], minentries=3 )
# # df = pd.DataFrame(SRF_CONT_RP_3R_myers2)
# # df.to_csv('~/Documents/Honours/SRF_myers2RP_3r.csv')
# # SRF_CONT_RP_3R_snyder1 = contTable(SRF_POS, [srf3, srf4, srf1], minentries=3 )
# # df = pd.DataFrame(SRF_CONT_RP_3R_snyder1)
# # df.to_csv('~/Documents/Honours/SRF_snyder1RP_3r.csv')
# # SRF_CONT_RP_3R_snyder2 = contTable(SRF_POS, [srf3, srf4, srf2], minentries=3 )
# # df = pd.DataFrame(SRF_CONT_RP_3R_snyder2)
# # df.to_csv('~/Documents/Honours/SRF_snyder2RP_3r.csv')
# #
# # SRF_CONT_RP_min2_3R_myers1 = contTable(SRF_POS, [srf1, srf2, srf3], minentries=2 )
# # df = pd.DataFrame(SRF_CONT_RP_min2_3R_myers1)
# # df.to_csv('~/Documents/Honours/SRF_myers1RP_3r_min2.csv')
# # SRF_CONT_RP_min2_3R_myers2 = contTable(SRF_POS, [srf1, srf2, srf4], minentries=2 )
# # df = pd.DataFrame(SRF_CONT_RP_min2_3R_myers2)
# # df.to_csv('~/Documents/Honours/SRF_myers2RP_3r_min2.csv')
# # SRF_CONT_RP_min2_3R_snyder1 = contTable(SRF_POS, [srf3, srf4, srf1], minentries=2 )
# # df = pd.DataFrame(SRF_CONT_RP_min2_3R_snyder1)
# # df.to_csv('~/Documents/Honours/SRF_snyder1RP_3r_min2.csv')
# # SRF_CONT_RP_min2_3R_snyder2 = contTable(SRF_POS, [srf3, srf4, srf2], minentries=2 )
# # df = pd.DataFrame(SRF_CONT_RP_min2_3R_snyder2)
# # df.to_csv('~/Documents/Honours/SRF_snyder2RP_3r_min2.csv')
# #
# # SRF_CONT_RP_myers = contTable(SRF_POS, [srf1, srf2], minentries=2 )
# # df = pd.DataFrame(SRF_CONT_RP_myers)
# # df.to_csv('~/Documents/Honours/SRF_myersRP.csv')
# # SRF_CONT_RP_snyder = contTable(SRF_POS, [srf3, srf4], minentries=2 )
# # df = pd.DataFrame(SRF_CONT_RP_snyder)
# # df.to_csv('~/Documents/Honours/SRF_snyderRP.csv')
# # # SRF_CONT_RP_myers = contTable(SRF_POS, [srf5, srf6], minentries=2 )
# # # df = pd.DataFrame(SRF_CONT_RP_myers)
# # # df.to_csv('~/Documents/Honours/SRF_myers2RP.csv')
# # #
# # # SRF_CONT_RP_myers = contTable(SRF_POS, [srf1, srf2, srf5, srf6], minentries=2 )
# # # df = pd.DataFrame(SRF_CONT_RP_myers)
# # # df.to_csv('~/Documents/Honours/SRF_myers2RP_min2.csv')
# # #
# # # SRF_CONT_RP_myers = contTable(SRF_POS, [srf3, srf4, srf5, srf6], minentries=2 )
# # # df = pd.DataFrame(SRF_CONT_RP_myers)
# # # df.to_csv('~/Documents/Honours/SRF_snyder2RP_min2.csv')
# #
# #
# # SRF_CONT_RP_min2 = contTable(SRF_NEG, srf, minentries=2 )
# # df = pd.DataFrame(SRF_CONT_RP_min2)
# # df.to_csv('~/Documents/Honours/SRF_RP_min2_NEG.csv')
# # SRF_CONT_RP_min3 = contTable(SRF_NEG, srf, minentries=3 )
# # df = pd.DataFrame(SRF_CONT_RP_min3)
# # df.to_csv('~/Documents/Honours/SRF_RP_min3_NEG.csv')
# # SRF_CONT_RP_min4 = contTable(SRF_NEG, srf, minentries=4 )
# # df = pd.DataFrame(SRF_CONT_RP_min4)
# # df.to_csv('~/Documents/Honours/SRF_RP_min4_NEG.csv')
# #
# # SRF_CONT_RP_3R_myers1 = contTable(SRF_NEG, [srf1, srf2, srf3], minentries=3 )
# # df = pd.DataFrame(SRF_CONT_RP_3R_myers1)
# # df.to_csv('~/Documents/Honours/SRF_myers1RP_3r_NEG.csv')
# # SRF_CONT_RP_3R_myers2 = contTable(SRF_NEG, [srf1, srf2, srf4], minentries=3 )
# # df = pd.DataFrame(SRF_CONT_RP_3R_myers2)
# # df.to_csv('~/Documents/Honours/SRF_myers2RP_3r_NEG.csv')
# # SRF_CONT_RP_3R_snyder1 = contTable(SRF_NEG, [srf3, srf4, srf1], minentries=3 )
# # df = pd.DataFrame(SRF_CONT_RP_3R_snyder1)
# # df.to_csv('~/Documents/Honours/SRF_snyder1RP_3r_NEG.csv')
# # SRF_CONT_RP_3R_snyder2 = contTable(SRF_NEG, [srf3, srf4, srf2], minentries=3 )
# # df = pd.DataFrame(SRF_CONT_RP_3R_snyder2)
# # df.to_csv('~/Documents/Honours/SRF_snyder2RP_3r_NEG.csv')
# #
# # SRF_CONT_RP_min2_3R_myers1 = contTable(SRF_NEG, [srf1, srf2, srf3], minentries=2 )
# # df = pd.DataFrame(SRF_CONT_RP_min2_3R_myers1)
# # df.to_csv('~/Documents/Honours/SRF_myers1RP_3r_min2_NEG.csv')
# # SRF_CONT_RP_min2_3R_myers2 = contTable(SRF_NEG, [srf1, srf2, srf4], minentries=2 )
# # df = pd.DataFrame(SRF_CONT_RP_min2_3R_myers2)
# # df.to_csv('~/Documents/Honours/SRF_myers2RP_3r_min2_NEG.csv')
# # SRF_CONT_RP_min2_3R_snyder1 = contTable(SRF_NEG, [srf3, srf4, srf1], minentries=2 )
# # df = pd.DataFrame(SRF_CONT_RP_min2_3R_snyder1)
# # df.to_csv('~/Documents/Honours/SRF_snyder1RP_3r_min2_NEG.csv')
# # SRF_CONT_RP_min2_3R_snyder2 = contTable(SRF_NEG, [srf3, srf4, srf2], minentries=2 )
# # df = pd.DataFrame(SRF_CONT_RP_min2_3R_snyder2)
# # df.to_csv('~/Documents/Honours/SRF_snyder2RP_3r_min2_NEG.csv')
# #
# # SRF_CONT_RP_myers = contTable(SRF_NEG, [srf1, srf2], minentries=2 )
# # df = pd.DataFrame(SRF_CONT_RP_myers)
# # df.to_csv('~/Documents/Honours/SRF_myersRP_NEG.csv')
# # SRF_CONT_RP_snyder = contTable(SRF_NEG, [srf3, srf4], minentries=2 )
# # df = pd.DataFrame(SRF_CONT_RP_snyder)
# # df.to_csv('~/Documents/Honours/SRF_snyderRP_NEG.csv')
# #
# # SRF_CONT_myersIDR = contTable(SRF_POS, srf, IDR = srf_myers_idr_run, minentries=4 )
# # df = pd.DataFrame(SRF_CONT_myersIDR)
# # df.to_csv('~/Documents/Honours/SRF_myersIDR.csv')
# # SRF_CONT_snyderIDR = contTable(SRF_POS, srf, IDR = srf_snyder_idr_run, minentries=4 )
# # df = pd.DataFrame(SRF_CONT_snyderIDR)
# # df.to_csv('~/Documents/Honours/SRF_snyderIDR.csv')
# #
# # SRF_CONT_myersIDR = contTable(SRF_NEG, srf, IDR = srf_myers_idr_run, minentries=4 )
# # df = pd.DataFrame(SRF_CONT_myersIDR)
# # df.to_csv('~/Documents/Honours/SRF_myersIDR_NEG.csv')
# # SRF_CONT_snyderIDR = contTable(SRF_NEG, srf, IDR = srf_snyder_idr_run, minentries=4 )
# # df = pd.DataFrame(SRF_CONT_snyderIDR)
# # df.to_csv('~/Documents/Honours/SRF_snyderIDR_NEG.csv')
#
# """BY PEAK"""
# # #
# # REST_CONT_min2Binom = contTable(REST_POS, rest, bypeak=True, minentries=2 )
# # df = pd.DataFrame(REST_CONT_min2Binom)
# # df.to_csv('~/Documents/Honours/REST_RP_min2_BP.csv')
# # REST_CONT_min3Binom = contTable(REST_POS, rest, bypeak=True, minentries=3 )
# # df = pd.DataFrame(REST_CONT_min3Binom)
# # df.to_csv('~/Documents/Honours/REST_RP_min3_BP.csv')
# # REST_CONT_min4Binom = contTable(REST_POS, rest, bypeak=True, minentries=4 )
# # df = pd.DataFrame(REST_CONT_min4Binom)
# # df.to_csv('~/Documents/Honours/REST_RP_min4_BP.csv')
# #
# # REST_CONT_RP_snyder = contTable(REST_POS, [rest1, rest2], bypeak=True, minentries=2 )
# # df = pd.DataFrame(REST_CONT_RP_snyder)
# # df.to_csv('~/Documents/Honours/REST_snyderRP_BP.csv')
# # REST_CONT_RP_myers = contTable(REST_POS, [rest3, rest4], bypeak=True, minentries=2 )
# # df = pd.DataFrame(REST_CONT_RP_myers)
# # df.to_csv('~/Documents/Honours/REST_myersRP_BP.csv')
# #
# # REST_CONT_RP_snyder1_3r = contTable(REST_POS, [rest1, rest2, rest3], bypeak=True, minentries=3 )
# # df = pd.DataFrame(REST_CONT_RP_snyder1_3r)
# # df.to_csv('~/Documents/Honours/REST_snyder1RP_3r_BP.csv')
# # REST_CONT_RP_snyder2_3r = contTable(REST_POS, [rest1, rest2, rest4], bypeak=True, minentries=3 )
# # df = pd.DataFrame(REST_CONT_RP_snyder2_3r)
# # df.to_csv('~/Documents/Honours/REST_snyder2RP_3r_BP.csv')
# # REST_CONT_RP_myers1_3r = contTable(REST_POS, [rest3, rest4, rest1], bypeak=True, minentries=3 )
# # df = pd.DataFrame(REST_CONT_RP_myers1_3r)
# # df.to_csv('~/Documents/Honours/REST_myers1RP_3r_BP.csv')
# # REST_CONT_RP_myers2_3r = contTable(REST_POS, [rest3, rest4, rest2], bypeak=True, minentries=3 )
# # df = pd.DataFrame(REST_CONT_RP_myers2_3r)
# # df.to_csv('~/Documents/Honours/REST_myers2RP_3r_BP.csv')
# #
# # REST_CONT_RP_min2_snyder1_3r = contTable(REST_POS, [rest1, rest2, rest3], bypeak=True, minentries=2 )
# # df = pd.DataFrame(REST_CONT_RP_min2_snyder1_3r)
# # df.to_csv('~/Documents/Honours/REST_snyder1RP_3r_min2_BP.csv')
# # REST_CONT_RP_min2_snyder2_3r = contTable(REST_POS, [rest1, rest2, rest4], bypeak=True, minentries=2 )
# # df = pd.DataFrame(REST_CONT_RP_min2_snyder2_3r)
# # df.to_csv('~/Documents/Honours/REST_snyder2RP_3r_min2_BP.csv')
# # REST_CONT_RP_min2_myers1_3r = contTable(REST_POS, [rest3, rest4, rest1], bypeak=True, minentries=2 )
# # df = pd.DataFrame(REST_CONT_RP_min2_myers1_3r)
# # df.to_csv('~/Documents/Honours/REST_myers1RP_3r_min2_BP.csv')
# # REST_CONT_RP_min2_myers2_3r = contTable(REST_POS, [rest3, rest4, rest2], bypeak=True, minentries=2 )
# # df = pd.DataFrame(REST_CONT_RP_min2_myers2_3r)
# # df.to_csv('~/Documents/Honours/REST_myers2RP_3r_min2_BP.csv')
# #
# #
# # REST_CONT_min2Binom = contTable(REST_NEG, rest, bypeak=True, minentries=2 )
# # df = pd.DataFrame(REST_CONT_min2Binom)
# # df.to_csv('~/Documents/Honours/REST_RP_min2_NEG_BP.csv')
# # REST_CONT_min3Binom = contTable(REST_NEG, rest, bypeak=True, minentries=3 )
# # df = pd.DataFrame(REST_CONT_min3Binom)
# # df.to_csv('~/Documents/Honours/REST_RP_min3_NEG_BP.csv')
# # REST_CONT_min4Binom = contTable(REST_NEG, rest, bypeak=True, minentries=4 )
# # df = pd.DataFrame(REST_CONT_min4Binom)
# # df.to_csv('~/Documents/Honours/REST_RP_min4_NEG_BP.csv')
# #
# # REST_CONT_RP_snyder = contTable(REST_NEG, [rest1, rest2], bypeak=True, minentries=2 )
# # df = pd.DataFrame(REST_CONT_RP_snyder)
# # df.to_csv('~/Documents/Honours/REST_snyderRP_NEG_BP.csv')
# # REST_CONT_RP_myers = contTable(REST_NEG, [rest3, rest4], bypeak=True, minentries=2 )
# # df = pd.DataFrame(REST_CONT_RP_myers)
# # df.to_csv('~/Documents/Honours/REST_myersRP_NEG_BP.csv')
# #
# # REST_CONT_RP_snyder1_3r = contTable(REST_NEG, [rest1, rest2, rest3], bypeak=True, minentries=3 )
# # df = pd.DataFrame(REST_CONT_RP_snyder1_3r)
# # df.to_csv('~/Documents/Honours/REST_snyder1RP_3r_NEG_BP.csv')
# # REST_CONT_RP_snyder2_3r = contTable(REST_NEG, [rest1, rest2, rest4], bypeak=True, minentries=3 )
# # df = pd.DataFrame(REST_CONT_RP_snyder2_3r)
# # df.to_csv('~/Documents/Honours/REST_snyder2RP_3r_NEG_BP.csv')
# # REST_CONT_RP_myers1_3r = contTable(REST_NEG, [rest3, rest4, rest1], bypeak=True, minentries=3 )
# # df = pd.DataFrame(REST_CONT_RP_myers1_3r)
# # df.to_csv('~/Documents/Honours/REST_myers1RP_3r_NEG_BP.csv')
# # REST_CONT_RP_myers2_3r = contTable(REST_NEG, [rest3, rest4, rest2], bypeak=True, minentries=3 )
# # df = pd.DataFrame(REST_CONT_RP_myers2_3r)
# # df.to_csv('~/Documents/Honours/REST_myers2RP_3r_NEG_BP.csv')
# #
# # REST_CONT_RP_min2_snyder1_3r = contTable(REST_NEG, [rest1, rest2, rest3], bypeak=True, minentries=2 )
# # df = pd.DataFrame(REST_CONT_RP_min2_snyder1_3r)
# # df.to_csv('~/Documents/Honours/REST_snyder1RP_3r_min2_NEG_BP.csv')
# # REST_CONT_RP_min2_snyder2_3r = contTable(REST_NEG, [rest1, rest2, rest4], bypeak=True, minentries=2 )
# # df = pd.DataFrame(REST_CONT_RP_min2_snyder2_3r)
# # df.to_csv('~/Documents/Honours/REST_snyder2RP_3r_min2_NEG_BP.csv')
# # REST_CONT_RP_min2_myers1_3r = contTable(REST_NEG, [rest3, rest4, rest1], bypeak=True, minentries=2 )
# # df = pd.DataFrame(REST_CONT_RP_min2_myers1_3r)
# # df.to_csv('~/Documents/Honours/REST_myers1RP_3r_min2_NEG_BP.csv')
# # REST_CONT_RP_min2_myers2_3r = contTable(REST_NEG, [rest3, rest4, rest2], bypeak=True, minentries=2 )
# # df = pd.DataFrame(REST_CONT_RP_min2_myers2_3r)
# # df.to_csv('~/Documents/Honours/REST_myers2RP_3r_min2_NEG_BP.csv')
#
#
# REST_CONT_snyderIDR = contTable(REST_POS, rest, IDR = rest_snyder_idr_run, bypeak=True, minentries=4 )
# df = pd.DataFrame(REST_CONT_snyderIDR)
# df.to_csv('~/Documents/Honours/REST_snyderIDR_BP.csv')
# REST_CONT_myersIDR = contTable(REST_POS, rest, IDR = rest_myers_idr_run, bypeak=True, minentries=4 )
# df = pd.DataFrame(REST_CONT_myersIDR)
# df.to_csv('~/Documents/Honours/REST_myersIDR_BP.csv')
#
# REST_CONT_snyderIDR = contTable(REST_NEG, rest, IDR = rest_snyder_idr_run, bypeak=True, minentries=4 )
# df = pd.DataFrame(REST_CONT_snyderIDR)
# df.to_csv('~/Documents/Honours/REST_snyderIDR_NEG_BP.csv')
# REST_CONT_myersIDR = contTable(REST_NEG, rest, IDR = rest_myers_idr_run, bypeak=True, minentries=4 )
# df = pd.DataFrame(REST_CONT_myersIDR)
# df.to_csv('~/Documents/Honours/REST_myersIDR_NEG_BP.csv')
# # scipy.stats.chi2_contingency([REST_CONT_snyderIDR[0][0], REST_CONT_myersIDR[0][0]])
# #
# # REST_ROC = ROC([REST_POS, REST_NEG], rest, bypeak=True, minentries=2, alpha=np.arange(0, 0.99, 0.01))
# # REST_ROC_snyder = ROC([REST_POS, REST_NEG], [rest1, rest2], bypeak=True, minentries=2, alpha=np.arange(0.01, 0.5, 0.01))
# # #
# # MAX_CONT_RP_min2Binom = contTable(MAX_POS, maxlist, bypeak=True, minentries=2 )
# # df = pd.DataFrame(MAX_CONT_RP_min2Binom)
# # df.to_csv('~/Documents/Honours/MAX_RP_min2_BP.csv')
# #
# #
# # MAX_CONT_RP_min3Binom = contTable(MAX_POS, maxlist, bypeak=True, minentries=3 )
# # df = pd.DataFrame(MAX_CONT_RP_min3Binom)
# # df.to_csv('~/Documents/Honours/MAX_RP_min3_BP.csv')
# # MAX_CONT_RP_min4Binom = contTable(MAX_POS, maxlist, bypeak=True, minentries=4 )
# # df = pd.DataFrame(MAX_CONT_RP_min4Binom)
# # df.to_csv('~/Documents/Honours/MAX_RP_min4_BP.csv')
# #
# #
# # MAX_CONT_RP_myers = contTable(MAX_POS, [max1, max2], bypeak=True, minentries=2 )
# # df = pd.DataFrame(MAX_CONT_RP_myers)
# # df.to_csv('~/Documents/Honours/MAX_myersRP_BP.csv')
# # MAX_CONT_RP_snyder = contTable(MAX_POS, [max3, max4], bypeak=True, minentries=2 )
# # df = pd.DataFrame(MAX_CONT_RP_snyder)
# # df.to_csv('~/Documents/Honours/MAX_snyderRP_BP.csv')
# #
# # MAX_CONT_RP_3r_myers1 = contTable(MAX_POS, [max1, max2, max3], bypeak=True, minentries=3 )
# # df = pd.DataFrame(MAX_CONT_RP_3r_myers1)
# # df.to_csv('~/Documents/Honours/MAX_myers1RP_3r_BP.csv')
# # MAX_CONT_RP_3r_myers2 = contTable(MAX_POS, [max1, max2, max4], bypeak=True, minentries=3 )
# # df = pd.DataFrame(MAX_CONT_RP_3r_myers2)
# # df.to_csv('~/Documents/Honours/MAX_myers2RP_3r_BP.csv')
# # MAX_CONT_RP_3r_snyder1 = contTable(MAX_POS, [max3, max4, max1], bypeak=True, minentries=3 )
# # df = pd.DataFrame(MAX_CONT_RP_3r_snyder1)
# # df.to_csv('~/Documents/Honours/MAX_snyder1RP_3r_BP.csv')
# # MAX_CONT_RP_3r_snyder2 = contTable(MAX_POS, [max3, max4, max2], bypeak=True, minentries=3 )
# # df = pd.DataFrame(MAX_CONT_RP_3r_snyder2)
# # df.to_csv('~/Documents/Honours/MAX_snyder2RP_3r_BP.csv')
# # #
# # MAX_CONT_RP_min2_3r_myers1 = contTable(MAX_POS, [max1, max2, max3], bypeak=True, minentries=2 )
# # df = pd.DataFrame(MAX_CONT_RP_min2_3r_myers1)
# # df.to_csv('~/Documents/Honours/MAX_myers1RP_3r_min2_BP.csv')
# # MAX_CONT_RP_min2_3r_myers2 = contTable(MAX_POS, [max1, max2, max4], bypeak=True, minentries=2 )
# # df = pd.DataFrame(MAX_CONT_RP_min2_3r_myers2)
# # df.to_csv('~/Documents/Honours/MAX_myers2RP_3r_min2_BP.csv')
# # MAX_CONT_RP_min2_3r_snyder1 = contTable(MAX_POS, [max3, max4, max1], bypeak=True, minentries=2 )
# # df = pd.DataFrame(MAX_CONT_RP_min2_3r_snyder1)
# # df.to_csv('~/Documents/Honours/MAX_snyder1RP_3r_min2_BP.csv')
# # MAX_CONT_RP_min2_3r_snyder2 = contTable(MAX_POS, [max3, max4, max2], bypeak=True, minentries=2 )
# # df = pd.DataFrame(MAX_CONT_RP_min2_3r_snyder2)
# # df.to_csv('~/Documents/Honours/MAX_snyder2RP_3r_min2_BP.csv')
# # #
# # #
# # # ##MAX NEGATIVES
# # MAX_CONT_RP_min2_NEG = contTable(MAX_NEG, maxlist, bypeak=True, minentries=2 )
# # df = pd.DataFrame(MAX_CONT_RP_min2_NEG)
# # df.to_csv('~/Documents/Honours/MAX_RP_min2_NEG_BP.csv')
# #
# # MAX_CONT_RP_min3_NEG = contTable(MAX_NEG, maxlist, bypeak=True, minentries=3 )
# # df = pd.DataFrame(MAX_CONT_RP_min3_NEG)
# # df.to_csv('~/Documents/Honours/MAX_RP_min3_NEG_BP.csv')
# #
# # MAX_CONT_RP_min4_NEG = contTable(MAX_NEG, maxlist, bypeak=True, minentries=4 )
# # df = pd.DataFrame(MAX_CONT_RP_min4_NEG)
# # df.to_csv('~/Documents/Honours/MAX_RP_min4_NEG_BP.csv')
# #
# # MAX_CONT_RP_myers_NEG = contTable(MAX_NEG, [max1, max2], bypeak=True, minentries=2 )
# # df = pd.DataFrame(MAX_CONT_RP_myers_NEG)
# # df.to_csv('~/Documents/Honours/MAX_myersRP_NEG_BP.csv')
# # MAX_CONT_RP_snyder = contTable(MAX_NEG, [max3, max4], bypeak=True, minentries=2 )
# # df = pd.DataFrame(MAX_CONT_RP_snyder)
# # df.to_csv('~/Documents/Honours/MAX_snyderRP_NEG_BP.csv')
# #
# # MAX_CONT_RP_3r_myers1 = contTable(MAX_NEG, [max1, max2, max3], bypeak=True, minentries=3 )
# # df = pd.DataFrame(MAX_CONT_RP_3r_myers1)
# # df.to_csv('~/Documents/Honours/MAX_myers1RP_3r_NEG_BP.csv')
# # MAX_CONT_RP_3r_myers2 = contTable(MAX_NEG, [max1, max2, max4], bypeak=True, minentries=3 )
# # df = pd.DataFrame(MAX_CONT_RP_3r_myers2)
# # df.to_csv('~/Documents/Honours/MAX_myers2RP_3r_NEG_BP.csv')
# # MAX_CONT_RP_3r_snyder1 = contTable(MAX_NEG, [max3, max4, max1], bypeak=True, minentries=3 )
# # df = pd.DataFrame(MAX_CONT_RP_3r_snyder1)
# # df.to_csv('~/Documents/Honours/MAX_snyder1RP_3r_NEG_BP.csv')
# # MAX_CONT_RP_3r_snyder2 = contTable(MAX_NEG, [max3, max4, max2], bypeak=True, minentries=3 )
# # df = pd.DataFrame(MAX_CONT_RP_3r_snyder2)
# # df.to_csv('~/Documents/Honours/MAX_snyder2RP_3r_NEG_BP.csv')
# # #
# # MAX_CONT_RP_min2_3r_myers1 = contTable(MAX_NEG, [max1, max2, max3], bypeak=True, minentries=2 )
# # df = pd.DataFrame(MAX_CONT_RP_min2_3r_myers1)
# # df.to_csv('~/Documents/Honours/MAX_myers1RP_3r_min2_NEG_BP.csv')
# # MAX_CONT_RP_min2_3r_myers2 = contTable(MAX_NEG, [max1, max2, max4], bypeak=True, minentries=2 )
# # df = pd.DataFrame(MAX_CONT_RP_min2_3r_myers2)
# # df.to_csv('~/Documents/Honours/MAX_myers2RP_3r_min2_NEG_BP.csv')
# # MAX_CONT_RP_min2_3r_snyder1 = contTable(MAX_NEG, [max3, max4, max1], bypeak=True, minentries=2 )
# # df = pd.DataFrame(MAX_CONT_RP_min2_3r_snyder1)
# # df.to_csv('~/Documents/Honours/MAX_snyder1RP_3r_min2_NEG_BP.csv')
# # MAX_CONT_RP_min2_3r_snyder2 = contTable(MAX_NEG, [max3, max4, max2], bypeak=True, minentries=2 )
# # df = pd.DataFrame(MAX_CONT_RP_min2_3r_snyder2)
# # df.to_csv('~/Documents/Honours/MAX_snyder2RP_3r_min2_NEG_BP.csv')
#
#
# # max_s1m1_idr_run = bed.BedFile(RankProd.collapse(bed.BedFile('ChIPdata/Max_K562/Snyder_hg19/MAX_s1m1_idrValues.bed', 'Peaks'), idr=True))
# # max_s1m2_idr_run = bed.BedFile(RankProd.collapse(bed.BedFile('ChIPdata/Max_K562/Snyder_hg19/MAX_s1m2_idrValues.bed', 'Peaks'), idr=True))
# # max_s2m1_idr_run = bed.BedFile(RankProd.collapse(bed.BedFile('ChIPdata/Max_K562/Snyder_hg19/MAX_s2m1_idrValues.bed', 'Peaks'), idr=True))
# # max_s2m2_idr_run = bed.BedFile(RankProd.collapse(bed.BedFile('ChIPdata/Max_K562/Snyder_hg19/MAX_s2m2_idrValues.bed', 'Peaks'), idr=True))
#
#
# MAX_CONT_myersIDR = contTable(MAX_POS, maxlist, IDR = max_myers_idr_run, bypeak=True, minentries=4 , Binom=False)
# df = pd.DataFrame(MAX_CONT_myersIDR)
# df.to_csv('~/Documents/Honours/MAX_myersIDR_BP.csv')
#
#
# MAX_CONT_snyderIDR = contTable(MAX_POS, maxlist, IDR = max_snyder_idr_run, bypeak=True, minentries=4 , Binom=False)
# df = pd.DataFrame(MAX_CONT_snyderIDR)
# df.to_csv('~/Documents/Honours/MAX_snyderIDR_BP.csv')
#
# MAX_CONT_myersIDR_NEG = contTable(MAX_NEG, maxlist, IDR = max_myers_idr_run, bypeak=True, minentries=4 , Binom=False)
# df = pd.DataFrame(MAX_CONT_myersIDR_NEG)
# df.to_csv('~/Documents/Honours/MAX_myersIDR_NEG_BP.csv')
#
#
# MAX_CONT_snyderIDR_NEG = contTable(MAX_NEG, maxlist, IDR = max_snyder_idr_run, bypeak=True, minentries=4 , Binom=False)
# df = pd.DataFrame(MAX_CONT_snyderIDR_NEG)
# df.to_csv('~/Documents/Honours/MAX_snyderIDR_NEG_BP.csv')
#
# # MAX_CONT_s1m1IDR = contTable(MAX_POS, maxlist, IDR = max_s1m1_idr_run, bypeak=True, minentries=4 , Binom=False)
# # df = pd.DataFrame(MAX_CONT_s1m1IDR)
# # df.to_csv('~/Documents/Honours/MAX_s1m1IDR_BP.csv')
# # MAX_CONT_s1m2IDR = contTable(MAX_POS, maxlist, IDR = max_s1m2_idr_run, bypeak=True, minentries=4 , Binom=False)
# # df = pd.DataFrame(MAX_CONT_s1m2IDR)
# # df.to_csv('~/Documents/Honours/MAX_s1m2IDR_BP.csv')
# # MAX_CONT_s2m1IDR = contTable(MAX_POS, maxlist, IDR = max_s2m1_idr_run, bypeak=True, minentries=4 , Binom=False)
# # df = pd.DataFrame(MAX_CONT_s2m1IDR)
# # df.to_csv('~/Documents/Honours/MAX_s2m1IDR_BP.csv')
# # MAX_CONT_s2m2IDR = contTable(MAX_POS, maxlist, IDR = max_s2m2_idr_run, bypeak=True, minentries=4 , Binom=False)
# # df = pd.DataFrame(MAX_CONT_s2m2IDR)
# # df.to_csv('~/Documents/Honours/MAX_s2m2IDR_BP.csv')
# MAX_CONT_1 = contTable(MAX_POS, maxlist, bypeak=True, minentries=1 , Binom=False)
# df = pd.DataFrame(MAX_CONT_s2m2IDR)
# df.to_csv('~/Documents/Honours/MAX_s2m2IDR_BP.csv')
#
# #
# # SRF_CONT_RP_min2 = contTable(SRF_POS, srf, bypeak=True, minentries=2 )
# # df = pd.DataFrame(SRF_CONT_RP_min2)
# # df.to_csv('~/Documents/Honours/SRF_RP_min2_BP.csv')
# # SRF_CONT_RP_min3 = contTable(SRF_POS, srf, bypeak=True, minentries=3 )
# # df = pd.DataFrame(SRF_CONT_RP_min3)
# # df.to_csv('~/Documents/Honours/SRF_RP_min3_BP.csv')
# # SRF_CONT_RP_min4 = contTable(SRF_POS, srf, bypeak=True, minentries=4 )
# # df = pd.DataFrame(SRF_CONT_RP_min4)
# # df.to_csv('~/Documents/Honours/SRF_RP_min4_BP.csv')
# #
# # SRF_CONT_RP_3R_myers1 = contTable(SRF_POS, [srf1, srf2, srf3], bypeak=True, minentries=3 )
# # df = pd.DataFrame(SRF_CONT_RP_3R_myers1)
# # df.to_csv('~/Documents/Honours/SRF_myers1RP_3r_BP.csv')
# # SRF_CONT_RP_3R_myers2 = contTable(SRF_POS, [srf1, srf2, srf4], bypeak=True, minentries=3 )
# # df = pd.DataFrame(SRF_CONT_RP_3R_myers2)
# # df.to_csv('~/Documents/Honours/SRF_myers2RP_3r_BP.csv')
# # SRF_CONT_RP_3R_snyder1 = contTable(SRF_POS, [srf3, srf4, srf1], bypeak=True, minentries=3 )
# # df = pd.DataFrame(SRF_CONT_RP_3R_snyder1)
# # df.to_csv('~/Documents/Honours/SRF_snyder1RP_3r_BP.csv')
# # SRF_CONT_RP_3R_snyder2 = contTable(SRF_POS, [srf3, srf4, srf2], bypeak=True, minentries=3 )
# # df = pd.DataFrame(SRF_CONT_RP_3R_snyder2)
# # df.to_csv('~/Documents/Honours/SRF_snyder2RP_3r_BP.csv')
# #
# # SRF_CONT_RP_min2_3R_myers1 = contTable(SRF_POS, [srf1, srf2, srf3], bypeak=True, minentries=2 )
# # df = pd.DataFrame(SRF_CONT_RP_min2_3R_myers1)
# # df.to_csv('~/Documents/Honours/SRF_myers1RP_3r_min2_BP.csv')
# # SRF_CONT_RP_min2_3R_myers2 = contTable(SRF_POS, [srf1, srf2, srf4], bypeak=True, minentries=2 )
# # df = pd.DataFrame(SRF_CONT_RP_min2_3R_myers2)
# # df.to_csv('~/Documents/Honours/SRF_myers2RP_3r_min2_BP.csv')
# # SRF_CONT_RP_min2_3R_snyder1 = contTable(SRF_POS, [srf3, srf4, srf1], bypeak=True, minentries=2 )
# # df = pd.DataFrame(SRF_CONT_RP_min2_3R_snyder1)
# # df.to_csv('~/Documents/Honours/SRF_snyder1RP_3r_min2_BP.csv')
# # SRF_CONT_RP_min2_3R_snyder2 = contTable(SRF_POS, [srf3, srf4, srf2], bypeak=True, minentries=2 )
# # df = pd.DataFrame(SRF_CONT_RP_min2_3R_snyder2)
# # df.to_csv('~/Documents/Honours/SRF_snyder2RP_3r_min2_BP.csv')
# #
# # SRF_CONT_RP_myers = contTable(SRF_POS, [srf1, srf2], bypeak=True, minentries=2 )
# # df = pd.DataFrame(SRF_CONT_RP_myers)
# # df.to_csv('~/Documents/Honours/SRF_myersRP_BP.csv')
# # SRF_CONT_RP_snyder = contTable(SRF_POS, [srf3, srf4], bypeak=True, minentries=2 )
# # df = pd.DataFrame(SRF_CONT_RP_snyder)
# # df.to_csv('~/Documents/Honours/SRF_snyderRP_BP.csv')
# # # SRF_CONT_RP_myers = contTable(SRF_POS, [srf5, srf6], bypeak=True, minentries=2 )
# # # df = pd.DataFrame(SRF_CONT_RP_myers)
# # # df.to_csv('~/Documents/Honours/SRF_myers2RP_BP.csv')
# # #
# # # SRF_CONT_RP_myers = contTable(SRF_POS, [srf1, srf2, srf5, srf6], bypeak=True, minentries=2 )
# # # df = pd.DataFrame(SRF_CONT_RP_myers)
# # # df.to_csv('~/Documents/Honours/SRF_myers2RP_min2_BP.csv')
# # #
# # # SRF_CONT_RP_myers = contTable(SRF_POS, [srf3, srf4, srf5, srf6], bypeak=True, minentries=2 )
# # # df = pd.DataFrame(SRF_CONT_RP_myers)
# # # df.to_csv('~/Documents/Honours/SRF_snyder2RP_min2_BP.csv')
# # #
# #
# # SRF_CONT_RP_min2 = contTable(SRF_NEG, srf, bypeak=True, minentries=2 )
# # df = pd.DataFrame(SRF_CONT_RP_min2)
# # df.to_csv('~/Documents/Honours/SRF_RP_min2_NEG_BP.csv')
# # SRF_CONT_RP_min3 = contTable(SRF_NEG, srf, bypeak=True, minentries=3 )
# # df = pd.DataFrame(SRF_CONT_RP_min3)
# # df.to_csv('~/Documents/Honours/SRF_RP_min3_NEG_BP.csv')
# # SRF_CONT_RP_min4 = contTable(SRF_NEG, srf, bypeak=True, minentries=4 )
# # df = pd.DataFrame(SRF_CONT_RP_min4)
# # df.to_csv('~/Documents/Honours/SRF_RP_min4_NEG_BP.csv')
# #
# # SRF_CONT_RP_3R_myers1 = contTable(SRF_NEG, [srf1, srf2, srf3], bypeak=True, minentries=3 )
# # df = pd.DataFrame(SRF_CONT_RP_3R_myers1)
# # df.to_csv('~/Documents/Honours/SRF_myers1RP_3r_NEG_BP.csv')
# # SRF_CONT_RP_3R_myers2 = contTable(SRF_NEG, [srf1, srf2, srf4], bypeak=True, minentries=3 )
# # df = pd.DataFrame(SRF_CONT_RP_3R_myers2)
# # df.to_csv('~/Documents/Honours/SRF_myers2RP_3r_NEG_BP.csv')
# # SRF_CONT_RP_3R_snyder1 = contTable(SRF_NEG, [srf3, srf4, srf1], bypeak=True, minentries=3 )
# # df = pd.DataFrame(SRF_CONT_RP_3R_snyder1)
# # df.to_csv('~/Documents/Honours/SRF_snyder1RP_3r_NEG_BP.csv')
# # SRF_CONT_RP_3R_snyder2 = contTable(SRF_NEG, [srf3, srf4, srf2], bypeak=True, minentries=3 )
# # df = pd.DataFrame(SRF_CONT_RP_3R_snyder2)
# # df.to_csv('~/Documents/Honours/SRF_snyder2RP_3r_NEG_BP.csv')
# #
# # SRF_CONT_RP_min2_3R_myers1 = contTable(SRF_NEG, [srf1, srf2, srf3], bypeak=True, minentries=2 )
# # df = pd.DataFrame(SRF_CONT_RP_min2_3R_myers1)
# # df.to_csv('~/Documents/Honours/SRF_myers1RP_3r_min2_NEG_BP.csv')
# # SRF_CONT_RP_min2_3R_myers2 = contTable(SRF_NEG, [srf1, srf2, srf4], bypeak=True, minentries=2 )
# # df = pd.DataFrame(SRF_CONT_RP_min2_3R_myers2)
# # df.to_csv('~/Documents/Honours/SRF_myers2RP_3r_min2_NEG_BP.csv')
# # SRF_CONT_RP_min2_3R_snyder1 = contTable(SRF_NEG, [srf3, srf4, srf1], bypeak=True, minentries=2 )
# # df = pd.DataFrame(SRF_CONT_RP_min2_3R_snyder1)
# # df.to_csv('~/Documents/Honours/SRF_snyder1RP_3r_min2_NEG_BP.csv')
# # SRF_CONT_RP_min2_3R_snyder2 = contTable(SRF_NEG, [srf3, srf4, srf2], bypeak=True, minentries=2 )
# # df = pd.DataFrame(SRF_CONT_RP_min2_3R_snyder2)
# # df.to_csv('~/Documents/Honours/SRF_snyder2RP_3r_min2_NEG_BP.csv')
# #
# # SRF_CONT_RP_myers = contTable(SRF_NEG, [srf1, srf2], bypeak=True, minentries=2 )
# # df = pd.DataFrame(SRF_CONT_RP_myers)
# # df.to_csv('~/Documents/Honours/SRF_myersRP_NEG_BP.csv')
# # SRF_CONT_RP_snyder = contTable(SRF_NEG, [srf3, srf4], bypeak=True, minentries=2 )
# # df = pd.DataFrame(SRF_CONT_RP_snyder)
# # df.to_csv('~/Documents/Honours/SRF_snyderRP_NEG_BP.csv')
#
# SRF_CONT_myersIDR = contTable(SRF_POS, srf, IDR = srf_myers_idr_run, bypeak=True, minentries=4 )
# df = pd.DataFrame(SRF_CONT_myersIDR)
# df.to_csv('~/Documents/Honours/SRF_myersIDR_BP.csv')
# SRF_CONT_snyderIDR = contTable(SRF_POS, srf, IDR = srf_snyder_idr_run, bypeak=True, minentries=4 )
# df = pd.DataFrame(SRF_CONT_snyderIDR)
# df.to_csv('~/Documents/Honours/SRF_snyderIDR_BP.csv')
#
# SRF_CONT_myersIDR = contTable(SRF_NEG, srf, IDR = srf_myers_idr_run, bypeak=True, minentries=4 )
# df = pd.DataFrame(SRF_CONT_myersIDR)
# df.to_csv('~/Documents/Honours/SRF_myersIDR_NEG_BP.csv')
# SRF_CONT_snyderIDR = contTable(SRF_NEG, srf, IDR = srf_snyder_idr_run, bypeak=True, minentries=4 )
# df = pd.DataFrame(SRF_CONT_snyderIDR)
# df.to_csv('~/Documents/Honours/SRF_snyderIDR_NEG_BP.csv')
# # rest=[]
# # maxlist=[]
# # srf=[]
# #
# # SRF_ROC_myers = ROC([SRF_POS, SRF_NEG], [srf1, srf2], minentries=2, alpha=np.arange(0.001, 0.01, 0.001))
# # SRF_ROC_snyder = ROC([SRF_POS, SRF_NEG], [srf3, srf4], minentries=2, alpha=np.arange(0.001, 0.01, 0.001))
# # SRF_ROC_4_min2 = ROC([SRF_POS, SRF_NEG], srf, minentries=2, alpha=np.arange(0.01, 0.4, 0.01))
# #
# # x1, y1 = ds.make_classification(10000, 4, 1, 3, 0, 2, 1, (0.5, 0.5), 0, 2, True, 15, 3, True)
#
#
#
#
#
#
# #JUNB ChIP-seq
# rep1 = bed.BedFile('ChIPdata/GSM2437872_ENCFF626BPP_peaks_GRCh38.bed', 'Peaks')
# rep2 = bed.BedFile('ChIPdata/GSM2437873_ENCFF056FWG_peaks_GRCh38.bed', 'Peaks')
# rep3 = bed.BedFile('ChIPdata/GSM2437874_ENCFF072SWK_peaks_GRCh38.bed', 'Peaks')
# rep4 = bed.BedFile('ChIPdata/GSM2437831_ENCFF895YDP_peaks_GRCh38.bed', 'Peaks')
# rep5 = bed.BedFile('ChIPdata/GSM2437832_ENCFF305VSD_peaks_GRCh38.bed', 'Peaks')
# rep6 = bed.BedFile('ChIPdata/GSM2437833_ENCFF245HYM_peaks_GRCh38.bed', 'Peaks')
# reps = [rep1, rep2, rep3, rep4, rep5, rep6]
# #junb = testcode.rankReps(rep1, rep2, rep3)
#
# results_reps = RankProd.performrankprod(reps, rankmethod='signalvalue',
#                                         minentries=2, duphandling='random',
#                                         alpha=0.05, filename='JUNB_unions.bed')
#
# """
# ChIA-PET Analysis
# """
#
#
# """
# MCF-7 CTCF hg19 ChIP and ChIA
# """
# ctcf_mcf7_chip1 = bed.BedFile('ChIPdata/MCF7_CTCF/ENCSR000DMR_VI/MCF7_CTCF_ChIP_hg19_r1.bed', 'Peaks')
# ctcf_mcf7_chip2 = bed.BedFile('ChIPdata/MCF7_CTCF/ENCSR000DMR_VI/MCF7_CTCF_ChIP_hg19_r2.bed', 'Peaks')
# ctcf_mcf7_chip3 = bed.BedFile('ChIPdata/MCF7_CTCF/ENCSR000AHD_RM/MCF7_CTCF_ChIP_hg19_r3.bed', 'Peaks')
# ctcf_mcf7_chip4 = bed.BedFile('ChIPdata/MCF7_CTCF/ENCSR000AHD_RM/MCF7_CTCF_ChIP_hg19_r4.bed', 'Peaks')
# ctcf = [ctcf_mcf7_chip1, ctcf_mcf7_chip2, ctcf_mcf7_chip3, ctcf_mcf7_chip4]
# ctcf_cp2_1 = bed.BedFile('CTCF_results/MCF7_CTCF_CP2_rep1_peaks.narrowPeak', 'Peaks')
# ctcf_cp2_2 = bed.BedFile('CTCF_results/MCF7_CTCF_CP2_rep2_peaks.narrowPeak', 'Peaks')
# ctcf_cp2_rp = RankProd.performrankprod([ctcf_cp2_1, ctcf_cp2_2], minentries=2, filename="MCF7_CTCF_CP2.bed")
#
# ctcf_results = RankProd.performrankprod(ctcf, minentries=4, filename="MCF7_CTCF_4.bed")
# ctcf_results_bed = bed.BedFile(ctcf_cp2_rp[0][0][0], 'Peaks')
#
# A549_TAD = bed.BedFile('TADs/A549_TADs.bed')
#
# T1_MCF7_CTCF_CP2 = bed.BedFile("T1_MCF7_CTCF_CP2.bed", "Peaks")
# MCF7_CTCF_CP2 = bed.BedFile("MCF7_CTCF_CP2.bed", "Peaks")
# T1_MCF7_CTCF_2 = bed.BedFile("T1_MCF7_CTCF_2.bed", "Peaks")
# T1_MCF7_CTCF_3 = bed.BedFile("T1_MCF7_CTCF_3.bed", "Peaks")
# T1_MCF7_CTCF_4 = bed.BedFile("T1_MCF7_CTCF_4.bed", "Peaks")
# ctcf_mcf7_chia_encodebed12 = bed.BedFile("ChIA_PET/MCF7_CTCF/MCF7_CTCF_ChIA_combined.bed", "BED12")
# ctcf_mcf7_micc = bed.BedFile("ChIA_PET/MCF7_CTCF/MCF7_CTCF_MICC.MICC", "BEDPE")
# ctcf_mcf7_chia_encodebedpe = bed.BED12toBEDPE(ctcf_mcf7_chia_encodebed12)
# ctcfmiccbed12 = bed.BEDPEtoBED12(ctcf_mcf7_micc)
#
# ctcf_encode_cp2 = RankProd.verifyChIAPeaks(ctcf_mcf7_chia_encodebed12, ctcf_results_bed)
# verified_ints = []
# for i in ctcf_encode_cp2:
#     if i[1] <= 0.05:
#         verified_ints.append(i[0])
#
# verifiedbed12 = bed.BEDPEtoBED12(verified_ints)
# bed.writeBedFile(verifiedbed12, 'verified_CTCF_MCF7.bed', 'BED12')
#
# ctcf_encode_cp2_bed12 = bed.BEDPEtoBED12(ctcf_encode_cp2)
# ctcf_encode_cp2_bed12 = bed.BedFile(ctcf_encode_cp2_bed12, 'BED12')
# bed.writeBedFile(ctcf_encode_cp2_bed12, 'CTCF_MCF7_CP2.bed', 'BED12')
# ctcf_bedpair = BedPairFile.BedPairFile(ctcf_encode_cp2_bed12)
#
# overlaps = []
# for i in ctcf_mcf7_chia_encodebedpe:
#     o = T1_MCF7_CTCF_CP2.getOverlap(i)
#     if o[0] is not None and o[1] is not None:
#         if len(o[0]) != 0 and len(o[1]) != 0:
#             overlaps.append(i)
#
# overlapsBED12 = bed.BEDPEtoBED12(overlaps)
# bed.writeBedFile(overlapsBED12, "MCF7_CTCF_ENCODEChIA_CP2.bed", "BED12")
#
# ctcfbed12_2 = bed.BedFile('CTCF_results/MCF7_CTCF_ENCODEChIA_2.bed', 'BED12')
# ctcfbedpe_2 = bed.BED12toBEDPE(ctcfbed12_2)
# ctcfbed12_3 = bed.BedFile('CTCF_results/MCF7_CTCF_ENCODEChIA_3.bed', 'BED12')
# ctcfbedpe_3 = bed.BED12toBEDPE(ctcfbed12_3)
# ctcfbed12_4 = bed.BedFile('CTCF_results/MCF7_CTCF_ENCODEChIA_4.bed', 'BED12')
# ctcfbedpe_4 = bed.BED12toBEDPE(ctcfbed12_4)
#
# bed.TADBoundary(K562_TADs, ctcfbedpe_2)
# bed.TADBoundary(K562_TADs, ctcfbedpe_3)
# bed.TADBoundary(K562_TADs, ctcfbedpe_4)
# bed.TADBoundary(K562_TADs, ctcf_mcf7_micc)
#
# c_2 = motif_discovery.bounded(ctcfbedpe_2)
# c_3 = motif_discovery.bounded(ctcfbedpe_3)
# c_4 = motif_discovery.bounded(ctcfbedpe_4)
# c_micc = motif_discovery.bounded(ctcf_mcf7_micc)
# t = [c_2[0], c_3[0], c_4[0], c_micc[0]]
# f = [c_2[1], c_3[1], c_4[1], c_micc[1]]
#
# p1 = plt.bar(np.arange(4), t, 0.35)
# p2 = plt.bar(np.arange(4), f, 0.35, bottom=t)
# plt.ylabel('Counts')
# plt.xticks(np.arange(4), ('C_2', 'C_3', 'C_4', 'C_micc'))
# plt.yticks(np.arange(0, 50000, 5000))
# plt.legend((p1[0], p2[0]), ('Bounded', 'Unbounded'))
# plt.title('CTCF MCF7 interactions bounded by K562 TADs')
#
# plt.show()
#
# p1 = plt.bar(np.arange(4), t, 0.35)
# # p2 = plt.bar(np.arange(4), f, 0.35)
# plt.ylabel('Counts')
# plt.xticks(np.arange(4), ('C_2', 'C_3', 'C_4', 'C_micc'))
# plt.yticks(np.arange(0, max(t)+5000, 5000))
# plt.legend('Bounded')
# plt.title('CTCF MCF7 interactions bounded by K562 TADs')
#
# plt.show()
#
# p2 = plt.bar(np.arange(4), f, 0.35)
# plt.ylabel('Counts')
# plt.xticks(np.arange(4), ('C_2', 'C_3', 'C_4', 'C_micc'))
# plt.yticks(np.arange(0, max(f)+500, 1000))
# plt.legend('Unbounded')
# plt.title('CTCF MCF7 interactions unbounded by K562 TADs')
#
# plt.show()
# """
# POLR2A K562 Analysis
# """
#
# POLR2A_chip1 = bed.BedFile('ChIPdata/compare/K562_POLR2A_ENCODE_hg19_rep1.bed', 'Peaks')
# POLR2A_chip2 = bed.BedFile('ChIPdata/compare/K562_POLR2A_ENCODE_hg19_rep2.bed', 'Peaks')
# POLR2A_chip3 = bed.BedFile('ChIPdata/K562_POLR2A_MSpe_rep1.regionPeak', 'Peaks')
# POLR2A_chip4 = bed.BedFile('ChIPdata/K562_POLR2A_MSpe_rep2.regionPeak', 'Peaks')
# polr2a_chip = [POLR2A_chip1, POLR2A_chip2, POLR2A_chip3, POLR2A_chip4]
#
# results2_POLR2A = RankProd.performrankprod(polr2a_chip, minentries=2, alpha=0.05, entAlpha=0.5, filename='POLR2A_K562_2.bed')
# results3_POLR2A = RankProd.performrankprod(polr2a_chip, minentries=3, alpha=0.05, entAlpha=0.5, filename='POLR2A_K562_3.bed')
# results4_POLR2A = RankProd.performrankprod(polr2a_chip, minentries=4, alpha=0.05, entAlpha=0.5, filename='POLR2A_K562_4.bed')
#
# T1_POLR2A_K562_2 = bed.BedFile('T1_POLR2A_K562_2.bed', 'Peaks')
# T1_POLR2A_K562_3 = bed.BedFile('T1_POLR2A_K562_3.bed', 'Peaks')
# T1_POLR2A_K562_4 = bed.BedFile('T1_POLR2A_K562_4.bed', 'Peaks')
# polr2a_micc = bed.BedFile('ChIPdata/POLR2A_MICC.bedpe', 'BEDPE')
# polr2a_k562_encodebed12 = bed.BedFile('ChIA_PET/K562_POLR2A_ChIA.bed', 'BED12')
# polr2a_k562_bedpe = bed.BED12toBEDPE(polr2a_k562_encodebed12)
# K562_TADs = bed.BedFile("K562_TADS_properchroms.bed")
#
# overlaps_2 = []
# for i in polr2a_k562_bedpe:
#     o = T1_POLR2A_K562_2.getOverlap(i)
#     if o[0] is not None and o[1] is not None:
#         if len(o[0]) != 0 and len(o[1]) != 0:
#             overlaps_2.append(i)
#
# bed.writeBedFile(overlaps_2, 'POLR2A_Overlapped_2', 'BEDPE')
# o12_2 = bed.BEDPEtoBED12(overlaps_2)
# bed.writeBedFile(o12_2, "POLR2A_K562_ENCODEChIA_2.bed", 'BED12')
#
# o12_3 = bed.BEDPEtoBED12(overlaps_3)
# bed.writeBedFile(o12_3, "POLR2A_K562_ENCODEChIA_3.bed", 'BED12')
#
# o12_4 = bed.BEDPEtoBED12(overlaps_4)
# bed.writeBedFile(o12_4, "POLR2A_K562_ENCODEChIA_4.bed", 'BED12')
#
#
# bed.TADBoundary(K562_TADs, overlaps_2)
# bed.TADBoundary(K562_TADs, overlaps_3)
# bed.TADBoundary(K562_TADs, overlaps_4)
# bed.TADBoundary(K562_TADs, polr2a_micc)
#
# p_2 = motif_discovery.bounded(overlaps_2)
# p_3 = motif_discovery.bounded(overlaps_3)
# p_4 = motif_discovery.bounded(overlaps_4)
# p_micc = motif_discovery.bounded(polr2a_micc)
#
# t = [p_2[0], p_3[0], p_4[0], p_micc[0]]
# f = [p_2[1], p_3[1], p_4[1], p_micc[1]]
#
# p1 = plt.bar(np.arange(4), t, 0.35)
# p2 = plt.bar(np.arange(4), f, 0.35, bottom=t)
# plt.ylabel('Counts')
# plt.xticks(np.arange(4), ('C_2', 'C_3', 'C_4', 'C_micc'))
# plt.yticks(np.arange(0, 35000, 5000))
# plt.legend((p1[0], p2[0]), ('Bounded', 'Unbounded'))
# plt.title('POLR2A K562 interactions bounded by K562 TADs')
#
# plt.show()
#
# p2 = plt.bar(np.arange(4), f, 0.35)
# plt.ylabel('Counts')
# plt.xticks(np.arange(4), ('C_2', 'C_3', 'C_4', 'C_micc'))
# plt.yticks(np.arange(0, max(f)+1000, 1000))
# plt.legend('Unbounded')
# plt.title('POLR2A K562 interactions unbounded by K562 TADs')
#
# plt.show()
#
#
#
# mm10_1 = bed.readBedFile('TH_mm10/ENCFF055WIO.bed', 'Peaks')
# mm10_2 = bed.readBedFile('TH_mm10/ENCFF736MPD.bed', 'Peaks')
#
# CP_peaks_1 = bed.BedFile('/home/rhys/Documents/Data/SNY_K562_POL2/SNY_K562_POL2_1_All_Peaks.narrowPeak', 'Peaks')
# CP_peaks_2 = bed.BedFile('/home/rhys/Documents/Data/SNY_K562_POL2/SNY_K562_POL2_2_All_Peaks.narrowPeak', 'Peaks')
# CP_peaks_3 = bed.BedFile('/home/rhys/Documents/Data/SNY_K562_POL2/SNY_K562_POL2_peaks.narrowPeak', 'Peaks')
# CP_peaks_4 = bed.BedFile('/home/rhys/Documents/Data/SNY_K562_POL2/R2.bed', 'Peaks')
#
# newpeaks1 = []
# for peak in CP_peaks_1:
#     if peak.qValue >= 1.30102999566:
#         newpeaks1.append(peak)
# newpeaks1 = bed.BedFile(newpeaks1, 'peaks')
# newpeaks2 = []
# for peak in CP_peaks_2:
#     if peak.qValue >= 1.30102999566:
#         newpeaks2.append(peak)
# newpeaks2 = bed.BedFile(newpeaks2, 'peaks')
#
# RPCP = RankProd.performrankprod([CP_peaks_1, CP_peaks_2, CP_peaks_3], minentries=2, filename="SNY_K562_POL2_RP.bed")
# T2_RPCP = bed.BedFile('T2_SNY_K562_POL2_RP.bed', 'RP')
# T1_RPCP = bed.BedFile('T1_SNY_K562_POL2_RP.bed', 'Peaks')
#
# ALL_RPCP = bed.BedFile('ALL_SNY_K562_POL2_RP.bed', 'RP')
#
# CP_idr = bed.BedFile('/home/rhys/Documents/Data/SNY_K562_POL2/idrValues.txt', 'idr')
# CP_idr = RankProd.collapse(CP_idr, idr=True)
# CP_idr_540 = getTopIDR(CP_idr, sc=5)
#
#
# CP_interactions_1 = bed.BedFile('/home/rhys/Documents/Data/SNY_K562_POL2/SNY_K562_POL2_1.interactions.all.mango', 'BEDPE')
# CP_interactions_2 = bed.BedFile('/home/rhys/Documents/Data/SNY_K562_POL2/SNY_K562_POL2_2.interactions.all.mango', 'BEDPE')
# CP_interactions = bed.BedFile('/home/rhys/Documents/Data/SNY_K562_POL2/SNY_K562_POL2.interactions.all.mango', 'BEDPE')
#
# CP_fdr_1 = bed.BedFile('/home/rhys/Documents/Data/SNY_K562_POL2/SNY_K562_POL2_1.interactions.fdr.mango', 'BEDPE')
# CP_fdr_2 = bed.BedFile('/home/rhys/Documents/Data/SNY_K562_POL2/SNY_K562_POL2_2.interactions.fdr.mango', 'BEDPE')
# CP_fdr = bed.BedFile('/home/rhys/Documents/Data/SNY_K562_POL2/SNY_K562_POL2.interactions.fdr.mango', 'BEDPE')
#
#
# CP_int_join = RankProd.joinChIA([CP_interactions_1, CP_interactions_2])
#
# CP_int_join = RankProd.joinChIA([CP_fdr_1, CP_fdr_1])
#
# distsAll = []
# for i in CP_int_join:
#     i.addOption(itemRgb='255,0,0', comment='All Mango')
#     distsAll.append(i.getDistance())
#
# # distsAllFdr = []
# # for i in CP_fdr_join:
# #     i.addOption(itemRgb='255,0,0', comment='All Mango')
# #     distsAll.append(i.getDistance())
#
#
# verifyIDR = RankProd.verifyChIAPeaks(CP_fdr, CP_idr_540, alpha=0.05, filename="IDR", mindist=10000)
# verifyIDRAll = verifyIDR[0]
# verifyIDRFDR = verifyIDR[1]
#
# distsIDR = []
# for i in verifyIDR[1]:
#     i.addOption(itemRgb='0,0,255', comment='verified IDR')
#     distsIDR.append(i.getDistance())
#
#
# verified = RankProd.ChIAreproducibility([CP_interactions], [ALL_RPCP], alpha=0.05)
# verifiedAll = verified[0]
# verifiedFDR = verified[1]
#
# verified4 = RankProd.ChIAreproducibility([CP_fdr_1, CP_fdr_2], [CP_peaks_1, CP_peaks_2], minentries=2)
#
#
# distsRP = []
# for i in verifiedFDR:
#     i.addOption(itemRgb='0,255,0', comment='verified RP')
#     distsRP.append(i.getDistance())
#
#
# verOverlap = overlap(verified[1], verifyIDR[1], ChIA=True)
#
#
# CP_int_join_100k = []
# for inter in CP_fdr:
#     if inter.getDistance() >= 20000 and inter not in CP_int_join_100k:
#         CP_int_join_100k.append(inter)
#
# RP4_100k = []
# for inter in verified4[1]:
#     if inter.getDistance() >= 20000:
#         RP4_100k.append(inter)
#
# RP2_100k = []
# RP2pv = []
# for inter in verified[1]:
#     if inter.getDistance() >= 100000:
#         RP2_100k.append(inter)
#         RP2pv.append(inter.pValue)
# sortRP = [x for _,x in sorted(zip(RP2pv, RP2_100k),key=lambda pair: pair[0])]
#
#
# IDR_100k = []
# IDRpv = []
# for inter in verifyIDR[1]:
#     if inter.getDistance() >= 20000:
#         IDR_100k.append(inter)
#         IDRpv.append(inter.pValue)
# sortIDR = [x for _,x in sorted(zip(IDRpv, IDR_100k),key=lambda pair: pair[0])]
#
#
# bed.writeBedFile(CP_int_join_100k, "All_POL2_contacts1.txt", "contact1", "chr1\tx1\tx2\tchr2\ty1\ty2\tcolor\tcomment")
# bed.writeBedFile(IDR_100k, "IDR_POL2_contacts1.txt", "contact1", "chr1\tx1\tx2\tchr2\ty1\ty2\tcolor\tcomment")
# bed.writeBedFile(RP2_100k, "RP2_POL2_contacts1.txt", "contact1", "chr1\tx1\tx2\tchr2\ty1\ty2\tcolor\tcomment")
# bed.writeBedFile(RP4_100k, "RP4_POL2_contacts1.txt", "contact1", "chr1\tx1\tx2\tchr2\ty1\ty2\tcolor\tcomment")
#
#
# CP_int_join = bed.BedFile(CP_int_join_100k, "BEDPE")
# idr = bed.BedFile(sortIDR, "BEDPE")
# rp = bed.BedFile(sortRP, "BEDPE")
# bed.writeBedFile(CP_int_join, "All_POL2_contacts2.txt", "contact2", "chr1\tx1\tx2\tchr2\ty1\ty2")
# bed.writeBedFile(idr, "IDR_POL2_contacts2.txt", "contact2", "chr1\tx1\tx2\tchr2\ty1\ty2")
# bed.writeBedFile(RP4_100k, "RP4_POL2_contacts2.txt", "contact2", "chr1\tx1\tx2\tchr2\ty1\ty2")
# bed.writeBedFile(rp, "RP2_POL2_contacts2.txt", "contact2", "chr1\tx1\tx2\tchr2\ty1\ty2")
#
# bed.writeBedFile(CP_int_join_100k[0:2000], "All_POL2_contacts2.txt", "contact2", "chr1\tx1\tx2\tchr2\ty1\ty2")
# bed.writeBedFile(verifyIDR[0][0:2000], "IDR_POL2_contacts2.txt", "contact2", "chr1\tx1\tx2\tchr2\ty1\ty2")
# bed.writeBedFile(verified4[2], "RP4_POL2_contacts2.txt", "contact2", "chr1\tx1\tx2\tchr2\ty1\ty2")
# bed.writeBedFile(verified[1], "RP2_POL2_contacts2.txt", "contact2", "chr1\tx1\tx2\tchr2\ty1\ty2")
#
#
