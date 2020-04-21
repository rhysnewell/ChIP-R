import unittest
from chipr import bed
from chipr.rankprod import *

class TestRankProdMethods(unittest.TestCase):

    def test_bed(self):
        bf = bed.BedFile("test/data/med/med1_peaks.broadPeak", "Peaks")
        print(bf.chroms.keys())
        g = bf.generate('chr1')
        print(next(g))
        print(next(g))
        print(next(g))
        cnt = 0
        for entry in bf:
            cnt += 1
            print(str(cnt) + '\t' + str(entry))
            if cnt == 100:
                break
        entry1 = bed.BedEntry('chrX', 3858266, 3858530)
        print(entry1 in bf)
        entry2 = bed.BedEntry('chrX', 10047550, 10067694)
        for x in bf.getOverlap(entry2):
            print(x)
        entry3 = bed.BedEntry('chr9', 102699903, 102700167)
        for x in bf.getClosest(entry3):
            print(x)
            for y in x:
                print(y)

    def test_entry_creation(self):
        entry1_1 = bed.BedEntry('X', 8000, 9000)
        entry1_1.signalValue = 10
        entry1_2 = bed.BedEntry('X', 80, 900)
        entry1_2.signalValue = 3
        bed1 = bed.BedFile([entry1_1, entry1_2])

        entry2_1 = bed.BedEntry('X', 8500, 9000)
        entry2_2 = bed.BedEntry('X', 80, 900)
        entry2_1.signalValue = 10
        entry2_2.signalValue = 3
        bed2 = bed.BedFile([entry2_1, entry2_2])

        entry3_1 = bed.BedEntry('X', 7500, 9000)
        entry3_2 = bed.BedEntry('X', 80, 900)
        entry3_1.signalValue = 10
        entry3_2.signalValue = 3
        bed3 = bed.BedFile([entry3_1, entry3_2])

        entry4_1 = bed.BedEntry('X', 5000, 8999)
        entry4_2 = bed.BedEntry('X', 80, 900)

        entry4_1.signalValue = 10
        entry4_2.signalValue = 3
        bed4 = bed.BedFile([entry4_1, entry4_2])

        unions = union([bed1, bed2, bed3, bed4], 2)
        for fragment in unions[0]:
            print(fragment)

    def test_single_chrom(self):
        med1 = bed.BedFile("test/data/med/med1_peaks.broadPeak", "Peaks").getChrom("chr17")
        med2 = bed.BedFile("test/data/med/med2_peaks.broadPeak", "Peaks").getChrom("chr17")
        med3 = bed.BedFile("test/data/med/med3_peaks.broadPeak", "Peaks").getChrom("chr17")

        bedf = [med1, med2, med3]
        minentries = 2

        # First create intersection and rank the entries in each replicate and return the rankproduct values
        ranks = rankreps(bedf, minentries, rankmethod='signalValue')

        # Calculate rank product for each entry that contributes to a union entry
        # Calculate the pvalues of the rank product values
        print('Calculating rank product probabilities...')
        rpb_up = rankprodbounds(ranks[1], len(ranks[1]), len(bedf), 'geometric')
        print('Calculating binomial threshold...')
        # Calculate rpb and binomial intersection point
        Pks = thresholdCalc(rpb_up, k=len(bedf) - (minentries - 1))
        if len(Pks[2]) != 0:
            binomAlpha = round(min(Pks[2]), 3)
        else:
            print('No binomial convergence, defaulting to 0.1')
            binomAlpha = 0.1

        # Perform multiple hypothesis testing correction upon the pvals
        fdr = multipletesting.fdrcorrection(rpb_up)

        # Determine whether to remove entries that are called significant
        print('Cleaning up output...')
        for i, v in enumerate(ranks[0][0]):
            p = rpb_up[i]

            if p != 0.0:
                ranks[0][0][i].addOption(name='TBD',
                                         score=min([abs(int(125 * math.log2(rpb_up[i]))), 1000]),
                                         strand='.',
                                         pValue=rpb_up[i],
                                         qValue=fdr[1][i])
            else:
                ranks[0][0][i].addOption(name='TBD',
                                         score=1000,
                                         strand='.',
                                         pValue=2.5e-20,
                                         qValue=2.5e-20)
        collapsed = bed.BedFile(ranks[0][0], 'IDR')

        connected = connect_entries(collapsed, bedf, 20, True)
        self.assertNotEqual(len(collapsed), len(connected))


    def test_from_bed(self):
        med1 = bed.BedFile("test/data/med/med1_peaks.broadPeak", "Peaks")
        med2 = bed.BedFile("test/data/med/med2_peaks.broadPeak", "Peaks")
        med3 = bed.BedFile("test/data/med/med3_peaks.broadPeak", "Peaks")

        bedf = [med1, med2, med3]

        performrankprod(bedf, minentries=2, rankmethod="pvalue", specifyMax=None,
                        duphandling='average', random_seed=0.5,
                        alpha=0.05,
                        filename="test_fragments_true",
                        default_min_peak=20,
                        print_pvals=True,
                        fragment=True)

        performrankprod(bedf, minentries=2, rankmethod="pvalue", specifyMax=None,
                        duphandling='average', random_seed=0.5,
                        alpha=0.05,
                        filename="test_fragments_false",
                        default_min_peak=20,
                        print_pvals=True,
                        fragment=False)

if __name__ == '__main__':
    unittest.main()