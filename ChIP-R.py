#!/usr/bin/python3

import argparse
import sys
import RankProd
import bed
import petprocessing

class RankProdAnalysis(object):

    def __init__(self):
        parser = argparse.ArgumentParser(prog='ChIP-R',
                                         description="Combine multiple ChIP-seq files and return a union of all peak "
                                                     "locations and a set confident, reproducible peaks as determined by "
                                                     "rank product analysis")
        parser.add_argument("-i", "--input",
                            help="ChIP-seq input files. These files must be in either narrowPeak, broadPeak, "
                                 "or regionPeak format. Multiple inputs a separeted by a single space",
                            dest="input",
                            type=argparse.FileType('r'),
                            nargs='+',
                            required=True)
        parser.add_argument("-o", "--output",
                            help="ChIP-seq output filename",
                            dest="output",
                            type=str,
                            default="rankprod",
                            required=False)
        parser.add_argument("-m", "--minentries",
                            help="The minimum peaks between replicates required to form an "
                                 "intersection of the peaks \n"
                                 "Default: 1",
                            dest="minentries",
                            default=1,
                            type=int,
                            required=False)
        parser.add_argument("--rankmethod",
                            help="The ranking method used to rank peaks within replicates. "
                                 "Options: 'signalvalue', 'pvalue', 'qvalue'. \n"
                                 "Default: pvalue",
                            dest="rankmethod",
                            default='pvalue',
                            type=str,
                            required=False)
        parser.add_argument("--duphandling",
                            help="Specifies how to handle entries that are ranked equally within a replicate "
                                 "Can either take the 'average' ranks or a 'random' rearrangement of the ordinal ranks \n"
                                 "Options: 'average', 'random' \n"
                                 "Default: 'average'",
                            dest="duphandling",
                            default='average',
                            type=str,
                            required=False)
        parser.add_argument("--seed",
                            help="Specify a seed to be used in conjunction with the 'random' option for -duphandling "
                                 "Must be between 0 and 1 \n"
                                 "default: 0.5",
                            dest="random_seed",
                            default=0.5,
                            type=float,
                            required=False)
        parser.add_argument("-a","--alpha",
                            help="Alpha specifies the user cut-off value for set of reproducible peaks "
                                 "The analysis will still produce results including peaks within the threshold calculated"
                                 "using the binomial method \n"
                                 "default: 0.05",
                            dest="alpha",
                            default=0.05,
                            type=float,
                            required=False)
        args = parser.parse_args(sys.argv[1:])
        self.run_1(args)

    def run_1(self, args):
        for i in args.input:
            print(str(i.name))

        print('Processing BedFiles...')
        bedfs = [bed.BedFile(str(i.name), 'Peaks') for i in args.input]

        RankProd.performrankprod(bedfs,
                                 args.minentries,
                                 args.rankmethod,
                                 'all',
                                 args.duphandling,
                                 args.random_seed,
                                 args.alpha,
                                 args.output)



if __name__ == "__main__":
    RankProdAnalysis()
