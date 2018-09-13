import ival
import bed
import statsmodels


"""
Python translation of A. Essebier's BedPair Java module. 16-Jan-18

To capture features in a pair of location representing a DNA loop secondary structure
"""
class BedPair:

    # def __init__(self, partner1, partner2, loop, range, chrom, id, intersects, pairedOverlaps, network):
    #     self.partner1 = partner1
    #     self.partner2 = partner2
    #     self.loop = loop
    #     self.range = range
    #     self.chrom = chrom
    #     self.id = id
    #     self.intersects = intersects
    #     self.pairedOverlaps = pairedOverlaps
    #     self.network = network

    def __init__(self, entry, padding=0):
        self.starts = [int(startword) for startword in entry.strBlockStarts.rstrip(',').split(',')]
        self.sizes = [int(sizeword) for sizeword in entry.strBlockSizes.rstrip(',').split(',')]
        if len(self.starts) != 2:
            raise IndexError("Entry %s:%d-%d does not have exactly two blocks and cannot build a pair. "
                             % (entry.chrom, entry.chromStart,entry.chromEnd))
        self.chrom = entry.chrom

        self.partner1 = bed.BedEntry(entry.chrom, ((entry.chromStart+self.starts[0])-padding),
                                ((entry.chromStart+self.starts[0]+self.sizes[0])+padding))
        self.partner2 = bed.BedEntry(entry.chrom, ((entry.chromStart+self.starts[1])-padding),
                                ((entry.chromStart+self.starts[1]+self.sizes[1])+padding))

        self.loop = bed.BedEntry(entry.chrom, self.partner1.chromEnd, self.partner2.chromStart)
        self.range = bed.BedEntry(entry.chrom, self.partner1.chromStart, self.partner2.chromEnd)

        self.id = entry.name

        # self.boundedTAD = None

    def setOverlap(self, intersects):
        self.intersects = intersects

    def setPairedOverlaps(self, state):
        self.pairedOverlaps = state

    def setPartner1(self, p1):
        self.partner1 = p1

    def setPartner2(self, p2):
        self.partner2 = p2

    # def TADBoundarySearch(self, TAD):
    #     if isinstance(TAD, bed.BedFile):
    #         searchChrom = TAD.generate(self.chrom)
    #         for i in searchChrom:





