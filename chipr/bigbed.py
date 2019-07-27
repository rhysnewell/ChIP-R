import pyBigWig
from chipr import bed, ival


def readBigBed(filename):

    file = pyBigWig.open(filename)
    chroms = dict()

    if file.isBigBed():
        if file.SQL().decode('utf8').lower().find('narrowpeak') != -1 or file.SQL().decode('utf8').lower().find(
                'broadpeak') != -1:
            for chrom in file.chroms():
                entries = file.entries(chrom, 0, file.chroms(chrom))
                for entry in entries:
                    try:
                        words = entry[2].strip().split()
                        chromStart = int(entry[0])
                        chromEnd = int(entry[1])
                        bed_entry = bed.BedEntry(chrom, chromStart, chromEnd)

                        if len(words) >= 7:  # narrowpeaks
                            bed_entry.addOption(name=words[0], score=int(words[1]), strand=words[2],
                                            signalValue=float(words[3]), pValue=float(words[4]),
                                            qValue=float(words[5]), peak=int(words[6]))
                        else:  # broadpeaks
                            bed_entry.addOption(name=words[0], score=int(words[1]), strand=words[2],
                                            signalValue=float(words[3]), pValue=float(words[4]),
                                            qValue=float(words[5]))

                        # check if the chromosome has been seen before
                        tree = chroms.get(chrom)
                        if not tree:
                            tree = ival.IntervalTree()
                            chroms[chrom] = tree
                        # put the entry in the interval tree for the appropriate chromosome
                        iv = ival.Interval(bed_entry.chromStart, bed_entry.chromEnd)
                        tree.put(iv, bed_entry)
                    except RuntimeError as e:
                        raise RuntimeError('Error in BIGBED file (%s)' % (e.strerror))

        else:
            print("BigBed file not ENCODE narrowPeak or broadPeak")
    file.close()
    return chroms