import ival
import bed
import BedPair

"""
Python translation of A. Essebier's BedPairFile Java module. 16-Jan-18

To load a file of DNA loop pairs based on chromosome conformation capture data
"""


class BedPairFile:

    # def __init__(self, collection, networks, graphs, rangeMap, pairMap):
    #     self.collection = []  # All pairs in original BedFile
    #     self.networks = set()  # All networks discovered in BedFile
    #     self.graphs = set()  # All graphs discovered in BedFile (corresponds to networks)
    #
    #     self.rangeMap = ival.IntervalTree  # Interval search tree for all pair ranges
    #     self.pairMap = ival.IntervalTree  # Interval search tree for all partners in pair

    collection = []  # All pairs in original BedFile
    networks = set()  # All networks discovered in BedFile
    graphs = set()  # All graphs discovered in BedFile (corresponds to networks)

    rangeMap = {}  # Interval search tree for all pair ranges
    pairMap = {}  # Interval search tree for all partners in pair

    LOGPSEUDO = 0.001

    genomeSizeHG19 = 3137161264
    genomeSizeHG38 = 3209286105
    genomeSizeMM9 = 2725765481
    genomeSizeMM10 = 2730871774
    supportedGenomes = ["hg19", "hg38", "mm9", "mm10"]

    defaultRGB = "0, 0, 128"
    rgbOptions = ["230,25,75", "60,180,75", "255,225,25",
                  "0,130,200", "245,130,48", "145,30,180",
                  "70,240,240", "240,50,230", "210,245,60", "250,190,190",
                  "0,128,128", "230,190,255", "170,110,40",
                  "170,255,195", "128,128,0", "255,215,180"]

    """
    Read a new collection of pairs from an existing bed file
    Excludes entries that don't have exactly 2 blocks
    Exclude entries spaced further apart than maxSize
    @param file filename
    @param maxSize the maximum genomic size (start to end) of a pair to include
    """
    def __init__(self, file, maxSize=9999999999, padding=0):

        data = bed.BedFile(file, 'bed12')
        for e in data:
            try:
                entry = BedPair.BedPair(e, padding)
                entry.partner1.name = entry.id+"_p1"
                entry.partner2.name = entry.id+"_p2"
                if entry.range.getLength() < maxSize:
                    self.collection.append(entry)
                    self.addEntry(entry)
            except IndexError:
                print("Entry %s:%d-%d does not have exactly two blocks and cannot build a pair. "
                      % (e.chrom, e.chromStart, e.chromEnd))
                pass

        # Had incorrect code for constructor
        # cnt = 0
        # for pair in self.collection:
        #     if pair.range.getLength() < maxSize:
        #         pair.setPartner1(
        #             bed.BedEntry(pair.chrom, pair.partner1.chromStart - padding, pair.partner1.chromEnd + padding))
        #         pair.setPartner2(
        #             bed.BedEntry(pair.chrom, pair.partner2.chromStart - padding, pair.partner2.chromEnd + padding))
        #         self.collection.append(pair)
        #         self.addEntry(pair)
        #         cnt += 1
        #         if cnt%1000 == 0:
        #             print(cnt)

    def __len__(self):
        return len(self.collection)

    """
    Given pair, add each partner to pair intervalST and range to range intervalST
    @param entry pair to be added to search trees
    """
    def addEntry(self, entry):
        self.addEntryPair(entry.partner1, entry)
        self.addEntryPair(entry.partner2, entry)
        self.addEntryRange(entry.range, entry)

    """
    Add partner as BedEntry with values as BedPair to track origin. This updates the pair tree which stores regions
    covered by partners of pair
    @param entry partner from a pair
    @param source BedPair which partner originates from
    """
    def addEntryPair(self, entry, source):
        chrom = entry.chrom
        tree = self.pairMap.get(chrom)
        if tree is None:
            tree = ival.IntervalTree()
            self.pairMap[chrom] = tree
        tree.put(ival.Interval(entry.chromStart, entry.chromEnd), source)

    """
    Add range as BedEntry with value as BedPair to track origin. This updates the range tree which tracks
    regions covered by a pair
    @param entry range of a pair
    @param source BedPair which partner originates from
    """
    def addEntryRange(self, entry, source):
        chrom = entry.chrom
        tree = self.rangeMap.get(chrom)
        if tree is None:
            tree = ival.IntervalTree()
            self.rangeMap[chrom] = tree
        tree.put(ival.Interval(entry.chromStart, entry.chromEnd), source)

    """
    Determine whether there are pairs in the bed pair file intersecting the query BedPair entry
    Use: given a query, determine whether other pairs intersect range, partner1 and partner2
    @param entry query entry
    @param padding used up and downstream to extend all features
    @return list of booleans describing types of intersects/overlaps
    Range overlap -> entry overlaps another pair in file
    P1 overlap -> P1 overlaps another partner from pair in file
    P2 overlap -> P2 overlaps another partner from pair in file
    """
    # def isIntersecting(self, entry, padding=0):
    #     pairTree = self.pairMap.get(entry.chrom)
    #     p1 =

    """
    Get all overlapping pairs from BedPairFile for each partner and range of query pair. Exclude self overlaps.
    Padding allows extension of features as margin of error. Where perfect alignment does not occur,
    select intersects within padding margin of error/tolerance.
    Use: to determine network paths, each pair must have record of other overlapping pairs.
    Use: find all pairs that could interact with query pair in some way.
    @param entry BedPair that overlaps are discovered
    @return List containing three sets representing:
    pairs overlapping range
    pairs overlapping partner 1
    pairs overlapping partner 2
    """
    def getIntersecting(self, entry, padding=0):
        pairTree = self.pairMap.get(entry.chrom)
        p1 = self.getSearchAll(pairTree, entry.partner1.chromStart - padding, entry.partner1.chromEnd + padding)
        p2 = self.getSearchAll(pairTree, entry.partner2.chromStart - padding, entry.partner2.chromEnd + padding)

        rangeTree = self.rangeMap.get(entry.chrom)
        r = self.getSearchAll(rangeTree, entry.range.chromStart - padding, entry.range.chromEnd + padding)

        r.remove(entry)
        p1.remove(entry)
        p2.remove(entry)

        return r, p1, p2


    """
    Given a search tree, range, and BedPair, determin whether range overlaps interval in search tree.
    Used to detect intersects/overlaps
    @param tree - tree to search in
    @param start - start point of search interval
    @param end - end point of search interval
    @param entry - identifier of self excluded from results - only care about overlaps with something other than self
    @return boolean - true if intersect/overlap detected
    """
    def treeSearchAll(self, tree, start, end, entry):
        ivals = tree.isectall(ival.Interval(start, end))
        if ivals is not None:
            foundall = None
            for interval in ivals:
                foundall = tree.get(interval)
            if entry is not None:
                foundall.remove(entry)
            if len(foundall) > 0:
                return True
        else:
            return False

    """
    Given a search tree, range, and BedPair, determine whether range overlaps interval in search tree.
    Used to detect intersects/overlaps
    @param tree - tree to search
    @param start - start point of search interval
    @param end - end point of search interval
    @return set of values from search tree overlapping provided interval
    """
    def getSearchAll(self, tree, start, end):
        ivals = tree.isectall(ival.Interval(start, end))
        if ivals is not None:
            ret = set()
            for i in ivals:
                foundall = tree.get(i)
                ret.add(foundall)
            return ret
        else:
            return set()

    """
    Calculate how much of given genome is covered by provided BedFile
    @param feature - BedFile containing regions to calculate coverage of
    @param genome - options: "hg19", "hg38", "mm9", "mm10"
    @return double - proportion of genome covered by feature
    """
    def getCoverage(self, feature, genome):
        bpCount = 0
        prop = 0
        rIter = iter(feature)
        while True:
            try:
                bpCount += bed._getLength(next(rIter))
            except StopIteration:
                break
        if genome == "hg19":
            prop = bpCount/self.genomeSizeHG19

        if genome == "hg38":
            prop = bpCount/self.genomeSizeHG38

        if genome == "mm9":
            prop = bpCount/self.genomeSizeMM9

        if genome == "mm10":
            prop = bpCount/self.genomeSizeMM10

        return prop

    def getExpressionInfo(self, expressionFile, ranges, maxDist, threshold, type="Closest"):
        line = ""
        hits = set()
        expressed = False
        maxScore = 0.0

        if type == "Closest":
            hits = expressionFile.getClosest(ranges)
            distance = bed.dist(hits, next(iter(ranges))) + self.LOGPSEUDO
            genes = list()
            for e in hits:
                if e.score > threshold:
                    expressed = True
                    genes.append(e)
                    if e.score > maxScore:
                        maxScore = e.score
                else:
                    genes.append("NA")
            if distance < maxDist:
                line += line.format("\t%f\t%s\t%f\t%s" % (distance, expressed, maxScore + self.LOGPSEUDO, ", ".join(genes)))
            else:
                line += "\tNA\tNA\tNA\tNA"
            return line

        if type == "Intersect":
            intersect = expressionFile.getOverlap(ranges)
            if len(intersect) == 0:
                line += "\t0\tNA\tNA\tNA"
            else:
                line += line.format("\t%d", len(intersect))
                for e in intersect:
                    if e.score > 0.0:
                        expressed = True
                        if e.score > maxScore:
                            maxScore = e.score
                line += line.format("\t%s\t%f" % (expressed, maxScore + self.LOGPSEUDO))
            return line

        if type == "Edge":
            closest = expressionFile.getOverlap(ranges)
            countExpressed = 0
            sumExpr = 0.0
            for e in closest:
                if e.score > threshold:
                    countExpressed += 1
                    sumExpr += e.score
            avgExpr = 0.0
            if countExpressed > 0:
                avgExpr = sumExpr/countExpressed
            line += line.format("\t%d\t%d\t%f" % (len(closest), countExpressed, avgExpr + self.LOGPSEUDO))
            return line
        else:
            print("Invalid expression info type")
            return None



