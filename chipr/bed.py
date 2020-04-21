from chipr import ival
import numpy as np
"""
Defines class for reading BED entries in the form of BED Paired End entries
"""


class BedPE:

    def __init__(self, chrom1, chromStart1, chromEnd1, chrom2, chromStart2, chromEnd2):
        self.chrom1 = chrom1
        self.chromStart1 = chromStart1
        self.chromEnd1 = chromEnd1
        self.chrom2 = chrom2
        self.chromStart2 = chromStart2
        self.chromEnd2 = chromEnd2
        self.blockCount = None
        self.usesstrand = False
        self.strand1 = None
        self.strand2 = None
        self.name1 = ''
        self.name2 = ''
        self.bounded = None
        self.partner1 = BedEntry(self.chrom1, self.chromStart1, self.chromEnd1)
        self.partner2 = BedEntry(self.chrom2, self.chromStart2, self.chromEnd2)

        if chrom1 == chrom2:
            self.loop = BedEntry(self.chrom1, self.chromEnd1, self.chromStart2)
            self.range = BedEntry(self.chrom1, self.chromStart1, self.chromEnd2)

        else:
            self.loop = -1
            self.range = -1

    def addOption(self,
                  name1=None,
                  name2=None,
                  score=None,
                  strand1=None,
                  strand2=None,
                  thickStart=None,
                  thickEnd=None,
                  itemRgb=None,
                  blockCount=None,
                  blockSizes=None,
                  blockStarts=None,
                  depth1=None,
                  depth2=None,
                  PETs=None,
                  signalValue=None,
                  pValue=None,
                  qValue=None,
                  peak=None,
                  tags=None,
                  summit1=None,
                  summit2=None,
                  fold=None,
                  fdr=None,
                  zscore=None,
                  bg=None, rank=None, comment=None):

        if name1: self.name = name1
        if name2: self.name = name2
        if score: self.score = score
        if strand1:
            self.strand1 = strand1
            self.usestrand = True  # use reverse complement when sequence is requested from genome
        if strand2:
            self.strand2 = strand2
            self.usestrand = True  # use reverse complement when sequence is requested from genome
        if thickStart: self.thickStart = thickStart
        if thickEnd: self.thickEnd = thickEnd
        if itemRgb: self.itemRgb = itemRgb
        if blockCount:
            self.blockCount = max(0, blockCount)
            if blockCount > 0:
                tblockSizes = [int(sizeword) for sizeword in blockSizes.rstrip(',').split(',')]
                tblockStarts = [int(startword) for startword in blockStarts.rstrip(',').split(',')]
                if len(tblockSizes) != blockCount or len(tblockStarts) != blockCount:
                    raise RuntimeError('Blockcount is incorrect in BED entry \"%s\"' % str(self))
                else:
                    self.blockSizes = tblockSizes
                    self.blockStarts = tblockStarts
        if depth1: self.depth1 = depth1
        if depth2: self.depth2 = depth2
        if PETs: self.PETs = PETs
        if signalValue: self.signalValue = signalValue
        if pValue: self.pValue = pValue
        if qValue: self.qValue = qValue
        if peak: self.peak = peak
        if tags: self.tags = tags
        if summit1: self.summit1 = summit1
        if summit2: self.summit2 = summit2
        if fold: self.fold = fold
        if fdr: self.fdr = fdr
        if bg: self.bg = bg
        if zscore: self.zscore = zscore
        if rank: self.rank = rank
        if comment: self.comment = comment

    def __str__(self):
        return str((self.chrom1, self.chromStart1, self.chromEnd1, self.chrom2, self.chromStart2, self.chromEnd2))

    def __getitem__(self, i):
        if self.blockCount:
            return (self.chrom1, self.chromStart1 + self.blockStarts[i],
                    self.chromStart1 + self.blockStarts[i] + self.blockSizes[i])

    def __iter__(self):
        if self.blockCount:
            for i in range(self.blockCount):
                if self.blockSizes[i] > 0:
                    yield (self.chrom1, self.chromStart1 + self.blockStarts[i],
                           self.chromStart1 + self.blockStarts[i] + self.blockSizes[i])

    def __len__(self):
        return self.blockCount

    def loc(self, genome=None, fixedwidth=None, usesummit=False, useshift=None):
        """ Retrieve the genomic location for BED entry, or sequence if genome is provided
            genome: a dictionary with keys for sequence names, e.g. 'chr1', 'chrX', etc, and values with indexed/sliceable strings
            fixedwidth: the width of the location/sequence if the width in the BED entry is ignored, and only its centre is used
            usesummit: centre a fixedwidth window around an assigned "summit"
            useshift: centre a fixedwidth window around a shifted centre point, e.g. useshift=-125 will shiftcentre point 125bp upstream,
            to say capture a fixedwidth=350bp window with 350/2-125=50bp downstream
        """
        otherstrand1 = False
        otherstrand2 = False
        if self.usesstrand:
            if self.strand1 == '-':
                otherstrand1 = True
            if self.strand2 == '-':
                otherstrand2 = True

        if otherstrand1 == False and otherstrand2 == False:
            end1 = self.chromEnd1
            start1 = self.chromStart1
            end2 = self.chromEnd2
            start2 = self.chromStart2
            mywidth1 = fixedwidth or (self.chromEnd1 - self.chromStart1)
            mycentre1 = start1 + (self.chromEnd1 - self.chromStart1) // 2
            mywidth2 = fixedwidth or (self.chromEnd2 - self.chromStart2)
            mycentre2 = start2 + (self.chromEnd2 - self.chromStart2) // 2
            if usesummit:
                mycentre1 = self.summit1
                mycentre2 = self.summit2
            if useshift:
                mycentre1 = mycentre1 + useshift
                mycentre2 = mycentre2 + useshift
            if fixedwidth:  # we need to re-calculate start and end
                if genome:
                    end1 = min(len(genome[self.chrom1]), mycentre1 + (mywidth1 // 2))
                    end2 = min(len(genome[self.chrom2]), mycentre2 + (mywidth2 // 2))
                else:
                    end1 = mycentre1 + (mywidth1 // 2)
                    end2 = mycentre2 + (mywidth2 // 2)
                start1 = max(0, mycentre1 - (mywidth1 // 2))
                start2 = max(0, mycentre2 - (mywidth2 // 2))

        #NOTE: TBC, need to add all strand configs
        else: # other strand
            start1 = self.chromEnd1
            end1 = self.chromStart1
            start2 = self.chromEnd2
            end2 = self.chromStart2
            mywidth = fixedwidth or (self.chromEnd1 - self.chromStart1)
            mycentre = self.chromStart1 + (self.chromEnd1 - self.chromStart1) // 2
            if usesummit:
                mycentre = self.summit
            if useshift:
                mycentre = mycentre - useshift  # shift is reversed on other strand
            if fixedwidth:  # we need to re-calculate start and end
                end = max(0, mycentre - (mywidth // 2))
                if genome:
                    start = min(len(genome[self.chrom1]), mycentre + (mywidth // 2))
                else:
                    start = mycentre + (mywidth // 2)

        if genome:  # refer to the genome sequence
            return [genome[self.chrom1][start1: end1], genome[self.chrom2][start2: end2]]
        else:
            return [self.chrom1, start1, end1, self.chrom2, start2, end2]

    def getInterval(self):
        return [ival.Interval(self.chromStart1, self.chromEnd1), ival.Interval(self.chromStart2, self.chromEnd2)]

    def getLength(self):
        return [self.chromEnd1 - self.chromStart1, self.chromEnd2 - self.chromStart2]

    def getDistance(self):
        if self.chrom1 == self.chrom2:
            d = abs(((self.chromStart1+self.chromEnd1)//2) - ((self.chromStart2+self.chromEnd2)//2))
            return d
        else:
            return None



"""
Defines class for reading in BED entries in a variety of formats
"""

class BedEntry:

    def __init__(self, chrom, chromStart, chromEnd):
        self.chrom = chrom
        self.chromStart = chromStart
        self.chromEnd = chromEnd
        self.blockCount = None
        self.usestrand = False
        self.strand = None
        self.name = ''

    def addOption(self,
                  name=None,
                  score=None,
                  strand=None,
                  thickStart=None,
                  thickEnd=None,
                  itemRgb=None,
                  blockCount=None,
                  blockSizes=None,
                  blockStarts=None,
                  strBlockSizes=None,
                  strBlockStarts=None,
                  signalValue=None,
                  pValue=None,
                  qValue=None,
                  peak=None,
                  localIDR=None,
                  globalIDR=None,
                  tags=None,
                  summit=None,
                  fold=None,
                  fdr=None,
                  zscore=None,
                  bg=None, rank=None,
                  rankProduct=None,
                  gene=None,
                  depth=None,
                  usestrand=False):
        if depth: self.depth = depth
        if name: self.name = name
        if score: self.score = score
        if strand:
            self.strand = strand
            self.usestrand = True # use reverse complement when sequence is requested from genome
        if usestrand: self.usestrand = usestrand
        if thickStart:
            if thickStart == 0:
                self.thickStart = "0"
            else:
                self.thickStart = thickStart
        if thickEnd: self.thickEnd = thickEnd
        if itemRgb: self.itemRgb = itemRgb
        if gene: self.gene = gene
        if blockCount:
            self.blockCount = max(0, blockCount)
            if blockCount > 0:
                tblockSizes = [int(sizeword) for sizeword in blockSizes.rstrip(',').split(',')]
                tblockStarts = [int(startword) for startword in blockStarts.rstrip(',').split(',')]
                if len(tblockSizes) != blockCount or len(tblockStarts) != blockCount:
                    raise RuntimeError('Blockcount is incorrect in BED entry \"%s\"' % str(self))
                else:
                    self.strBlockSizes = blockSizes
                    self.strBlockStarts = blockStarts
                    self.blockSizes = tblockSizes
                    self.blockStarts = tblockStarts

        if signalValue: self.signalValue = signalValue
        if pValue: self.pValue = pValue
        if qValue: self.qValue = qValue
        if peak: self.peak = peak
        if localIDR: self.localIDR = localIDR
        if globalIDR: self.globalIDR = globalIDR
        if tags: self.tags = tags
        if summit: self.summit = summit
        if fold: self.fold = fold
        if fdr: self.fdr = fdr
        if bg: self.bg = bg
        if zscore: self.zscore = zscore
        if rank: self.rank = rank
        if rankProduct: self.rankProduct = rankProduct

    def __str__(self):
        return str((self.chrom, self.chromStart, self.chromEnd))

    def __getitem__(self, i):
        if self.blockCount:
            return (self.chrom, self.chromStart + self.blockStarts[i], self.chromStart + self.blockStarts[i] + self.blockSizes[i])

    def __iter__(self):
        if self.blockCount:
            for i in range(self.blockCount):
                if self.blockSizes[i] > 0:
                    yield (self.chrom, self.chromStart + self.blockStarts[i], self.chromStart + self.blockStarts[i] + self.blockSizes[i])

    def __len__(self):
        return self.blockCount

    def loc(self, genome = None, fixedwidth = None, usesummit = False, useshift = None):
        """ Retrieve the genomic location for BED entry, or sequence if genome is provided
            genome: a dictionary with keys for sequence names, e.g. 'chr1', 'chrX', etc, and values with indexed/sliceable strings
            fixedwidth: the width of the location/sequence if the width in the BED entry is ignored, and only its centre is used
            usesummit: centre a fixedwidth window around an assigned "summit"
            useshift: centre a fixedwidth window around a shifted centre point, e.g. useshift=-125 will shiftcentre point 125bp upstream,
            to say capture a fixedwidth=350bp window with 350/2-125=50bp downstream
        """
        otherstrand = False
        if (self.usestrand):
            if (self.strand == '-'):
                otherstrand = True

        if (otherstrand == False):
            end = self.chromEnd
            start = self.chromStart
            mywidth = fixedwidth or (self.chromEnd - self.chromStart)
            mycentre = start + (self.chromEnd - self.chromStart) // 2
            if usesummit:
                mycentre = self.summit
            if useshift:
                mycentre = mycentre + useshift
            if fixedwidth: # we need to re-calculate start and end
                if genome:
                    end = min(len(genome[self.chrom]), mycentre + (mywidth // 2))
                else:
                    end = mycentre + (mywidth // 2)
                start = max(0, mycentre - (mywidth // 2))

        else: # other strand
            start = self.chromEnd
            end = self.chromStart
            mywidth = fixedwidth or (self.chromEnd - self.chromStart)
            mycentre = self.chromStart + (self.chromEnd - self.chromStart) // 2
            if usesummit:
                mycentre = self.summit
            if useshift:
                mycentre = mycentre - useshift # shift is reversed on other strand
            if fixedwidth: # we need to re-calculate start and end
                end = max(0, mycentre - (mywidth // 2))
                if genome:
                    start = min(len(genome[self.chrom]), mycentre + (mywidth // 2))
                else:
                    start = mycentre + (mywidth // 2)

        if genome: # refer to the genome sequence
            return genome[self.chrom][start : end]
        else:
            return (self.chrom, start, end)

    def setwidth(self, fixedwidth = None, usesummit = False):
        if fixedwidth:
            if usesummit:
                diff = self.summit - fixedwidth // 2
            else:
                diff = (self.chromEnd - self.chromStart) // 2 - fixedwidth // 2
            self.chromStart += diff
            self.chromStart += diff + fixedwidth
        return (self.chrom, self.chromStart, self.chromEnd)

    def getInterval(self):
        return ival.Interval(self.chromStart, self.chromEnd)

    def getLength(self):
        return self.chromEnd - self.chromStart




def _getLength(interval):
    return interval.chromEnd - interval.chromStart


def dist(entry1, entry2, signed = False, centre2centre = False):
    """ Calculate and return the BedEntry with the closest distance
        (from one end of the interval of this to the end of the interval of that).
        If centre2centre is True, use the centre-to-centre distance instead.
        If signed is True, the distance is negative if this interval is after the that.
    """
    if isinstance(entry1, BedEntry) and isinstance(entry2, BedEntry):
        if entry1.chrom == entry2.chrom:
            return ival.dist(entry1.getInterval(), entry2.getInterval(), signed, centre2centre)
    elif isinstance(entry1, BedPE) and isinstance(entry2, BedEntry):
        if entry1.chrom1 == entry1.chrom2 == entry2.chrom:
            return [ival.dist(entry1.getInterval()[0], entry2.getInterval(), signed, centre2centre),
                    ival.dist(entry1.getInterval()[1], entry2.getInterval(), signed, centre2centre)]
    elif isinstance(entry1, BedEntry) and isinstance(entry2, BedPE):
        if entry2.chrom1 == entry2.chrom2 == entry1.chrom:
            return [ival.dist(entry1.getInterval(), entry2.getInterval()[0], signed, centre2centre),
                    ival.dist(entry1.getInterval(), entry2.getInterval()[1], signed, centre2centre)]
    elif isinstance(entry1, BedPE) and isinstance(entry2, BedPE):
        if entry1.chrom1 == entry1.chrom2 == entry2.chrom1 == entry2.chrom2:
            return [ival.dist(entry1.getInterval()[0], entry2.getInterval()[0], signed, centre2centre),
                    ival.dist(entry1.getInterval()[0], entry2.getInterval()[1], signed, centre2centre),
                    ival.dist(entry1.getInterval()[1], entry2.getInterval()[0], signed, centre2centre),
                    ival.dist(entry1.getInterval()[1], entry2.getInterval()[1], signed, centre2centre)]
    return None


class BedFile:
    """ Read BED file.

        See http://genome.ucsc.edu/FAQ/FAQformat#format1

        The first three required BED fields are (part of all supported sub-formats):

        chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671).
        chromStart - The starting position of the feature in the chromosome or scaffold.
                     The first base in a chromosome is numbered 0.
        chromEnd - The ending position of the feature in the chromosome or scaffold.
                   The chromEnd base is not included in the display of the feature. .

        The 9 additional optional BED fields are (part of sub-format "Optional"):

        name - Defines the name of the BED line.
        score - A score between 0 and 1000.
                If the track line useScore attribute is set to 1 for this annotation data set,
                the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray).
        strand - Defines the strand - either '+' or '-'.
        thickStart - The starting position at which the feature is drawn thickly (for example, the start codon in gene displays).
        thickEnd - The ending position at which the feature is drawn thickly (for example, the stop codon in gene displays).
        itemRgb - An RGB value of the form R,G,B (e.g. 255,0,0).
                  If the track line itemRgb attribute is set to "On",
                  this RBG value will determine the display color of the data contained in this BED line.
        blockCount - The number of blocks (exons) in the BED line.
        blockSizes - A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount.
        blockStarts - A comma-separated list of block starts. All of the blockStart positions should be calculated relative to chromStart.
                      The number of items in this list should correspond to blockCount.

        ENCODE also defines broadpeaks and narrowpeaks format (part of our "Peaks" sub-format):

        name - Defines the name of the BED line.
        score - Indicates how dark the peak will be displayed in the browser (0-1000).
        strand - +/- to denote strand or orientation (whenever applicable). Use '.' if no orientation is assigned.
        signalValue - Measurement of overall (usually, average) enrichment for the region.
        pValue - Measurement of statistical significance (-log10). Use -1 if no pValue is assigned.
        qValue - Measurement of statistical significance using false discovery rate (-log10). Use -1 if no qValue is assigned.
        peak - Point-source called for this peak; 0-based offset from chromStart. Use -1 if no point-source called.

        MACS also defines a "summit" peaks format (part of our "Summit" sub-format)
        It contains the peak summits locations for every peaks. The 5th column in this file is the .
        In addition to the required three, the following fields follow:
        length         [redundant, ignored]
        summit         summit height of fragment pileup
        tags
        pValue         [-10*log10(pvalue)]
        fold           [enrichment]
        FDR            [%; optional]

        "CCAT" BED-like file format:
        chromosome,
        peakcenter    [converted to summit],
        regionstart,
        regionend,
        tags          [tagcount],
        bg            [bgcount],
        zscore,
        fdr

    """

    def __init__(self, entries, format = 'Limited'):
        """
        Create a BedFile instance.
        :param entries: an iterable of entries or a filename
        :param format: the format of the BED file
        """
        self.format = format

        if isinstance(entries, str): # filename
            self.chroms = readBedFile(entries, format)
        else:
            self.chroms = dict()
            if format.lower().startswith('bedpe'):
                for entry in entries:
                    for num in range(1, 3):
                        if num == 1:
                            tree = self.chroms.get(entry.chrom1)
                            if not tree:
                                tree = ival.IntervalTree()
                                self.chroms[entry.chrom1] = tree
                            iv = ival.Interval(entry.chromStart1, entry.chromEnd1)
                            tree.put(iv, entry)

                        elif num == 2:
                            tree = self.chroms.get(entry.chrom2)
                            if not tree:
                                tree = ival.IntervalTree()
                                self.chroms[entry.chrom2] = tree
                                # put the entry in the interval tree for the appropriate chromosome
                                iv = ival.Interval(entry.chromStart2, entry.chromEnd2)
                                tree.put(iv, entry)
            else:
                for entry in entries:
                    # check if the chromosome has been seen before
                    tree = self.chroms.get(entry.chrom)
                    if not tree:
                        tree = ival.IntervalTree()
                        self.chroms[entry.chrom] = tree
                    # put the entry in the interval tree for the appropriate chromosome
                    iv = ival.Interval(entry.chromStart, entry.chromEnd)
                    tree.put(iv, entry)

    def __len__(self):
        n = 0
        for c in self.chroms:
            n += len(self.chroms[c])
        return n

    def generate(self, chrom):
        mytree = self.chroms.get(chrom)
        if mytree != None:
            for e in mytree:
                for entry in e.values:
                    yield entry

    def getChrom(self, chrom):
        mytree = self.chroms.get(chrom)
        entries = []
        if mytree != None:
            for e in mytree:
                for entry in e.values:
                    entries.append(entry)
        return BedFile(entries, self.format)

    def __iter__(self):
        self.chromqueue = ival.Stack()
        for c in sorted(self.chroms.keys())[::-1]:
            self.chromqueue.push(self.generate(c))
        self.current = self.chromqueue.pop()
        return self

    def __next__(self):
        try:
            ret = next(self.current)
        except StopIteration:
            if not self.chromqueue.isEmpty():
                self.current = self.chromqueue.pop()
                ret = next(self.current)
            else:
                raise StopIteration
        return ret

    def __contains__(self, item):
        if isinstance(item, BedEntry):
            tree = self.chroms.get(item.chrom)
            if tree is None:
                return False
            else:
                return ival.Interval(item.chromStart, item.chromEnd) in tree

        elif isinstance(item, BedPE):
            tree1 = self.chroms.get(item.chrom1)
            tree2 = self.chroms.get(item.chrom2)
            container = []
            if tree1 is None:
                container.append(False)
            else:
                container.append(ival.Interval(item.chromStart1, item.chromEnd1) in tree1)
            if tree2 is None:
                container.append(False)
            else:
                container.append(ival.Interval(item.chromStart2, item. chromEnd2) in tree2)
            return container
        else:
            return False

    def getOverlap(self, item):
        if isinstance(item, BedEntry):
            tree = self.chroms.get(item.chrom)
            if tree is None: return None
            else:
                iv = ival.Interval(item.chromStart, item.chromEnd)
                res = tree.isectall(iv)
                ret = []
                for r in res:
                    ret.extend(r.values)
                return ret
        elif isinstance(item, BedPE):
            tree1 = self.chroms.get(item.chrom1)
            tree2 = self.chroms.get(item.chrom2)
            container = []
            if tree1 is None:
                container.append(None)
            else:
                iv = ival.Interval(item.chromStart1, item.chromEnd1)
                res = tree1.isectall(iv)
                ret = []
                for r in res:
                    ret.extend(r.values)
                container.append(ret)
            if tree2 is None:
                container.append(None)
            else:
                iv = ival.Interval(item.chromStart2, item.chromEnd2)
                res = tree2.isectall(iv)
                ret = []
                for r in res:
                    ret.extend(r.values)
                container.append(ret)
            return container
        else:
            return None

    def getClosest(self, item):
        if isinstance(item, BedEntry):
            tree = self.chroms.get(item.chrom)
            if tree is None:
                return None
            else:
                iv = ival.Interval(item.chromStart, item.chromEnd)
                node = tree.closest(iv)
                if node is not None:
                    return node.values
                else:
                    return None
        elif isinstance(item, BedPE):
            tree1 = self.chroms.get(item.chrom1)
            tree2 = self.chroms.get(item.chrom2)
            container = []
            if tree1 is None:
                container.append(None)
            else:
                iv = ival.Interval(item.chromStart1, item.chromEnd1)
                node = tree1.closest(iv)
                if node is not None:
                    container.append(node.values)
                else:
                    container.append(None)
            if tree2 is None:
                container.append(None)
            else:
                iv = ival.Interval(item.chromStart1, item.chromEnd1)
                node = tree2.closest(iv)
                if node is not None:
                    container.append(node.values)
                else:
                    container.append(None)
            return container
        else:
            return None


    def getOneOfClosest(self, item):
        all = self.getClosest(item)
        if all == None: return None
        else: return next(iter(all))

    def getOneOfOverlap(self, item):
        all = self.getOverlap(item)
        if all == None: return None
        elif len(all) == 0: return None
        else: return next(iter(all))

    """
    getMetrics()
    Method for BedFile class that returns the mean, standard deviation, median, and median absolute deviation
    of the peak sizes for each chromosome of the bedfile
    """
    def getMetrics(self):
        chroms = self.chroms.keys()
        metrics = {}
        all_widths = []
        for chr in chroms:
            chromarray = self.generate(chr)
            entry = next(chromarray, 0)
            widthlist = []
            while entry != 0:
                widthlist.append(entry.chromEnd - entry.chromStart)
                all_widths.append(entry.chromEnd - entry.chromStart)
                entry = next(chromarray, 0)

            std = np.std(widthlist)
            widthmean = np.mean(widthlist)
            # abdiff = []
            m = np.median(widthlist)
            # Median Absolute Deviation is not necessary
            # for j in widthlist:
            #     abdiff.append(abs(j - m))
            # MAD = np.median(abdiff)

            metrics[chr] = [widthmean, std, m]
        # Stats for all peaks
        std = np.std(all_widths)
        widthmean = np.mean(all_widths)
        m = np.median(all_widths)
        metrics['total'] = [widthmean, std, m]

        return metrics

    def poolBED(self, bedfile):
        for entry in bedfile:
            # check if the chromosome has been seen before
            tree = self.chroms.get(entry.chrom)
            if not tree:
                tree = ival.IntervalTree()
                self.chroms[entry.chrom] = tree
            # put the entry in the interval tree for the appropriate chromosome
            iv = ival.Interval(entry.chromStart, entry.chromEnd)
            tree.put(iv, entry)




def readBedFile(filename, format = 'Limited'):
    """ Read a BED file.
        format: specifies the format of the file,
        "Limited", e.g.
            chr22 1000 5000
            chr22 2000 6000
        "Optional", e.g.
            track name=pairedReads description="Clone Paired Reads" useScore=1
            chr22 1000 5000 cloneA 960 + 1000 5000 0 2 567,488, 0,3512
            chr22 2000 6000 cloneB 900 - 2000 6000 0 2 433,399, 0,3601
            ...
            (also handles the Limited + score, and BED6 format)
        "Peaks", e.g.
            chr1    569780    569930    .    0    .    19    6.07811    -1    -1
            chr1    713300    713450    .    0    .    54    49.1167    -1    -1
        "Strand", e.g.
            chr4    185772359    185772424    -
            chr18    20513381    20513401    +
        also supports a 5th label field
            chr5    20611949        20611949        +       ENSG00000251629_20611949
            chr3    42187863        42187863        -       ENSG00000234562_42187863
        "Summit", e.g.
            # d = 130
            chr      start    end   length summit  tags -10*log10(pvalue)    fold_enrichment    FDR(%)
            chr1     8250     8671    422    286    46    145.84    11.68    0.51
            chr1    36382    36984    603    405    46    315.23    27.05    0.24
        "CCAT", e.g.
            chr8    94747805    94747070    94749250    525     3    21.519196    0.002000
            chr17   55277895    55277070    55279280    560    18    21.283333    0.002000
        "Cropped", e.g.
            chr1    851602    10
            chr1    921184    18
            chr1    931838    9
        "BedPE", e.g.
            chrom1  chromStart1 chromEnd1   chrom2  chromStart2 chromEnd2 + any number of additional fields
            chr1    85617       86100       chr1    120030      125039
            chr2    73891       74871       chr5    12709       12990
    """
    f = open(filename)
    row = 0
    acceptHeaderRows = 1
    headerRow = None
    sissrs = False
    gem = False
    start = False
    chroms = dict()
    for line in f:
        row += 1
        words = line.strip().split()
        if len(words) == 0:
            continue
        if words[0].strip().startswith('='):
            sissrs = True
            continue
        if words[0].strip().startswith('Position'):
            gem = True
            continue
        if sissrs:
            if words[0].strip().startswith('-'):
                start = True
                continue
            elif start:
                chrom = words[0]
                chromStart = int(words[1])
                chromEnd = int(words[2])
                entry = BedEntry(chrom, chromStart, chromEnd)
                entry.addOption(signalValue=int(words[3]), name = " ", score = int(words[3]), strand = '.',
                                pValue = float(-1), qValue = float(-1), peak = int(-1))

                # check if the chromosome has been seen before
                tree = chroms.get(chrom)
                if not tree:
                    tree = ival.IntervalTree()
                    chroms[chrom] = tree
                # put the entry in the interval tree for the appropriate chromosome
                iv = ival.Interval(entry.chromStart, entry.chromEnd)
                tree.put(iv, entry)
            else:
                continue
        elif gem:
            chrom, centre = words[0].split(':')
            centre = int(centre)
            chromStart = centre - 50
            chromEnd = centre + 50
            entry = BedEntry(chrom, chromStart, chromEnd)
            entry.addOption(signalValue=float(words[1]), name=" ", score=float(words[7]), strand=words[13], peak=centre,
                            pValue=float(words[6]), qValue=float(words[5]))
            tree = chroms.get(chrom)
            if not tree:
                tree = ival.IntervalTree()
                chroms[chrom] = tree
            # put the entry in the interval tree for the appropriate chromosome
            iv = ival.Interval(entry.chromStart, entry.chromEnd)
            tree.put(iv, entry)
        else:
            if len(words) == 0:
                continue # ignore empty lines
            if words[0].strip().startswith('#'):
                continue # comment
            if words[0].strip().startswith('browser'):
                continue # ignore
            if words[0].strip().startswith('track'):
                continue # ignore
            if words[1].strip().startswith('start'):
                continue # ignore
            try:
                if format.lower().startswith('bedpe'):
                    chrom1 = words[0]
                    chromStart1 = int(words[1])
                    chromEnd1 = int(words[2])
                    chrom2 = words[3]
                    chromStart2 = int(words[4])
                    chromEnd2 = int(words[5])

                    entry = BedPE(chrom1, chromStart1, chromEnd1, chrom2, chromStart2, chromEnd2)
                    if len(words) == 8:
                        entry.addOption(PETs=int(words[6]), pValue=float(words[7]))
                    if len(words) == 13:
                        entry.addOption(name1=words[6], name2=words[7], depth1=int(words[8]), depth2=int(words[9]),
                                        PETs=int(words[10]), pValue=float(words[11]), fdr=float(words[12]))
                    if chrom1 == chrom2:
                        tree = chroms.get(chrom1)
                        if not tree:
                            tree = ival.IntervalTree()
                            chroms[chrom1] = tree
                        iv1 = ival.Interval(entry.chromStart1, entry.chromEnd1)
                        iv2 = ival.Interval(entry.chromStart2, entry.chromEnd2)
                        tree.put(iv1, entry)
                        tree.put(iv2, entry)

                    else:
                        tree1 = chroms.get(chrom1)
                        tree2 = chroms.get(chrom2)
                        if not tree1:
                            tree1 = ival.IntervalTree()
                            chroms[chrom1] = tree1
                        if not tree2:
                            tree2 = ival.IntervalTree()
                            chroms[chrom2] = tree2
                        # put the entry in the interval tree for the appropriate chromosome
                        iv1 = ival.Interval(entry.chromStart1, entry.chromEnd1)
                        iv2 = ival.Interval(entry.chromStart2, entry.chromEnd2)
                        tree1.put(iv1, entry)
                        tree2.put(iv2, entry)
                else:
                    chrom = words[0]
                    if format.lower().startswith('ccat'):
                        chromStart = int(words[2])
                        chromEnd = int(words[3])
                    else: # all other standard BED formats
                        try:
                            chromStart = int(words[1])
                            chromEnd = int(words[2])
                        except ValueError:
                            print(words)
                            continue
                    entry = BedEntry(chrom, chromStart, chromEnd)
                    if format.lower().startswith('opt'):
                        if len(words) >= 9:
                            entry.addOption(name = words[3], score = float(words[4]), strand = words[5],
                                            thickStart = int(words[6]), thickEnd = int(words[7]), itemRgb = words[8])
                        elif len(words) >= 6:
                            entry.addOption(name = words[3], score = float(words[4]), strand = words[5])
                        elif len(words) >= 5:
                            entry.addOption(name = words[3], score = float(words[4]))
                        elif len(words) >= 4:
                            entry.addOption(name = words[3])
                        else:
                            entry.addOption(name = '.', score = int(words[3]), strand = '.')
                    elif format.lower().startswith('bed6'):
                        entry.addOption(name=words[3], score=float(words[4]), strand=words[5])
                    elif format.lower().startswith('strand'):
                        if len(words) >= 4: # properly formatted
                            entry.addOption(strand = words[3])
                        if len(words) >= 5:
                            entry.addOption(name = words[4])
                    elif format.lower().startswith('peak'):
                        if len(words) >= 10: # narrowpeaks
                            entry.addOption(name = words[3], score = int(words[4]), strand = words[5],
                                            signalValue = float(words[6]), pValue = float(words[7]),
                                            qValue = float(words[8]), peak = int(words[9]))
                        else: # broadpeaks
                            entry.addOption(name = words[3], score = int(words[4]), strand = words[5],
                                            signalValue = float(words[6]), pValue = float(words[7]),
                                            qValue = float(words[8]))
                    elif format.lower().startswith('rp'):
                        entry.addOption(name=words[3], score=int(words[4]), strand=words[5],
                                        signalValue=float(words[6]), pValue=float(words[7]),
                                        qValue=float(words[8]), rank=[float(r) for r in list(words[9].split(","))])
                    elif format.lower().startswith('idr'):
                        entry.addOption(name=words[3], score=int(words[4]), strand=words[5],
                                        signalValue=float(words[6]), pValue=float(words[7]),
                                        qValue=float(words[8]))
                    elif format.lower().startswith('2idr'):
                        #For IDR input with actual IDR values
                        entry.addOption(name=words[3], pValue=float(words[9]),
                                        qValue=float(words[10]))
                    elif format.lower().startswith('summit'):
                        if len(words) >= 9:
                            entry.addOption(summit = int(words[4]), tags = int(words[5]), pValue = float(words[6]),
                                            fold = float(words[7]), fdr = float(words[8]))
                        else:
                            entry.addOption(summit = int(words[4]), tags = int(words[5]),
                                            pValue = float(words[6]), fold = float(words[7]))
                    elif format.lower().startswith('ccat'):
                        entry.addOption(summit = int(words[1]) - entry.chromStart, tags = int(words[4]), bg = int(words[5]),
                                        zscore = float(words[6]), fdr = float(words[7]), name = '.',
                                        score = int(words[4]), strand = '.')
                    elif format.lower().startswith('crop'):
                        entry.addOption(score = int(words[2]), name = '.', strand = '.')
                        entry.chromEnd = entry.chromStart + 1
                    elif format.lower().startswith('bed12'):
                        entry.addOption(name=words[3], score=float(words[4]), strand=words[5], thickStart=int(words[6]),
                                        thickEnd=int(words[7]), itemRgb=words[8], blockCount=int(words[9]),
                                        blockSizes=words[10], blockStarts=words[11])
                    elif format.lower().startswith('TSS'):
                        entry.addOption(name=str(words[3]), gene=str(words[4]), strand=words[5])
                    elif format.lower().startswith('mspc'):
                        entry.addOption(name=str(words[3]), signalValue=float(words[4]))

                    # check if the chromosome has been seen before
                    tree = chroms.get(chrom)
                    if not tree:
                        tree = ival.IntervalTree()
                        chroms[chrom] = tree
                    # put the entry in the interval tree for the appropriate chromosome
                    iv = ival.Interval(entry.chromStart, entry.chromEnd)
                    tree.put(iv, entry)
            except RuntimeError as e:
                if not acceptHeaderRows:
                    raise RuntimeError('Error in BED file at row %d (%s)' % (row, e.strerror))
                else:
                    headerRow = words
                    acceptHeaderRows -= 1 # count down the number of header rows that can occur
    f.close()
    return chroms


'''
BEDPEtoBED12 and BED12toBEDPE
Set of helper functions to convert bed files between bedpe and BED12 format
'''


def BEDPEtoBED12(bedpe):
    bed12 = []
    for peEnt in bedpe:
        if peEnt.chrom1 == peEnt.chrom2:
            start = min(peEnt.chromStart1, peEnt.chromStart2)
            if start == 0:
                start = 1
            end = max(peEnt.chromEnd1, peEnt.chromEnd2)
            bed12Ent = BedEntry(peEnt.chrom1, start, end)
            bSizes = [peEnt.chromEnd1 - peEnt.chromStart1, peEnt.chromEnd2 - peEnt.chromStart2]
            bStarts = [0, max(peEnt.chromStart1, peEnt.chromStart2) - start]
            name = peEnt.chrom1+":"+str(peEnt.chromStart1)+".."+str(peEnt.chromEnd1)+"-"+peEnt.chrom2+":"\
                   +str(peEnt.chromStart2)+".."+str(peEnt.chromEnd2)
            bed12Ent.addOption(name=name, thickStart=int(start), thickEnd=int(end), strand=".", blockCount=2, score=500,
                               itemRgb="255,0,0", blockSizes=str(bSizes[0])+","+str(bSizes[1]),
                               blockStarts=str(bStarts[0])+","+str(bStarts[1]))
            bed12.append(bed12Ent)
        else:
            bed12Ent = BedEntry(peEnt.chrom1, peEnt.chromStart1, peEnt.chromEnd1)
            name = peEnt.chrom1 + ":" + str(peEnt.chromStart1) + ".." + str(peEnt.chromEnd1) + "-" + peEnt.chrom2 + ":" \
                    + str(peEnt.chromStart2) + ".." + str(peEnt.chromEnd2)
            bed12Ent.addOption(name=name, score=5, blockCount=1, thickStart=peEnt.chromStart1, thickEnd=peEnt.chromEnd1,
                               strand=".", itemRgb="255,0,0", blockSizes=str(peEnt.chromEnd1-peEnt.chromStart1),
                               blockStarts=str(0))
            bed12.append(bed12Ent)
    return bed12


def BED12toBEDPE(bed12):
    bedpe = []
    for ent in bed12:

        if ent.blockCount > 1:
            blockSizes = ent.blockSizes
            chrom1 = ent.chrom
            chrom2 = ent.chrom
            start1 = ent.chromStart
            end1 = start1+blockSizes[0]

            end2 = ent.chromEnd
            start2 = ent.chromEnd - blockSizes[1]
            peEnt = BedPE(chrom1, start1, end1, chrom2, start2, end2)
            bedpe.append(peEnt)
        else:
            if len(ent.name.split(':')) == 2:
                name = ent.name
                name = name.split('-')
                chrom1 = name.split(':')[0][0]
                chrom2 = name.split(':')[1][0]
                start1 = name.split(':')[0][1].split('..')[0]
                end1 = name.split(':')[0][1].split('..')[1]
                start2 = name.split(':')[1][1].split('..')[0]
                end2 = name.split(':')[1][1].split('..')[1]
                peEnt = BedPE(chrom1, start1, end1, chrom2, start2, end2)
                bedpe.append(peEnt)
    return bedpe


"""
Function that determines whether the entries of a BedPE file are bounded within the entries of another file.
Only useful when comparing BedPE files to TADs BedFiles
"""
def TADBoundary(BedFile, BedPE):
    boundedBedPE = []
    for ent in BedPE:
        if ent.loop is not None:
            searchChrom = BedFile.generate(ent.loop.chrom)
            for tad in searchChrom:
                if tad.chromStart <= ent.loop.chromStart and ent.loop.chromEnd <= tad.chromEnd:
                    ent.bounded = True
                    boundedBedPE.append(ent)
                    break
                else:
                    ent.bounded = False
    return boundedBedPE

def writeBedFile(entries, filename, format = 'BED6', header = None):
    """ Save the BED entries to a BED file.
        format - the format to use for WRITING
    """
    f = open(filename, 'w+')
    if header:
        f.write(header + '\n')
    for row in entries:
        if format == 'Peaks':
            f.write("%s\t%d\t%d\t%s\t%d\t%s\t%f\t%.16f\t%.16f" %
                    (row.chrom, row.chromStart, row.chromEnd, row.name,
                     row.score, row.strand, row.signalValue, row.pValue, row.qValue))
        elif format == 'RP':
            f.write("%s\t%d\t%d\t%s\t%d\t%s\t%f\t%.16f\t%.16f\t%s" %
                    (row.chrom, row.chromStart, row.chromEnd, row.name,
                     row.score, row.strand, row.signalValue, row.pValue, row.qValue, ','.join(str(r) for r in row.rank)))
        elif format == 'Limited':
            f.write("%s\t%d\t%d" % (row.chrom, row.chromStart, row.chromEnd))
        elif format == 'Strand':
            f.write("%s\t%d\t%d\t%s\t$s\t%s" % (row.chrom, row.chromStart, row.chromEnd, row.strand, row.name))
        elif format == 'BED12':
            f.write("%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%s\t%d\t%s\t%s"
                                % (row.chrom, row.chromStart, row.chromEnd,
                                    row.name, row.score, row.strand,
                                    row.thickStart, row.thickEnd, row.itemRgb,
                                    row.blockCount, row.strBlockSizes, row.strBlockStarts))
        elif format == 'BEDPE':
            f.write("%s\t%d\t%d\t%s\t%d\t%d\t%16f" % (row.chrom1, row.chromStart1, row.chromStart1,
                                                      row.chrom2, row.chromStart2, row.chromEnd2, row.pValue))
        elif format == 'contact1':
            f.write("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s" % (row.chrom1, (row.chromStart1+row.chromEnd1)//2,
                                                        (row.chromStart2+row.chromEnd2)//2,
                                                      row.chrom2, (row.chromStart1+row.chromEnd1)//2,
                                                        (row.chromStart2+row.chromEnd2)//2, '0,255,0', ' '))
        elif format == 'contact2':
            f.write("%s\t%d\t%d\t%s\t%d\t%d" % (row.chrom1, row.chromStart1, row.chromEnd1,
                                                      row.chrom2, row.chromStart2,  row.chromEnd2))
        else:
            f.write("%s\t%d\t%d\t%s\t%d\t%s" %
                    (row.chrom, row.chromStart, row.chromEnd, row.name, row.score, row.strand))
        f.write("\n")
    f.close()