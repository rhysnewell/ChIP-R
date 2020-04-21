#!/usr/bin/python3

import argparse
import sys
from chipr import rankprod, bed



class RankProdAnalysis(object):

    def __init__(self):
        parser = argparse.ArgumentParser(prog='chipr',
                                         description="Combine multiple ChIP-seq files and return a union of all peak "
                                                     "locations and a set confident, reproducible peaks as determined by "
                                                     "rank product analysis")
        parser.add_argument("-i", "--input",
                            help="ChIP-seq input files. These files must be in either narrowPeak, broadPeak, "
                                 "or regionPeak format. Multiple inputs are separeted by a single space",
                            dest="input",
                            type=argparse.FileType('r'),
                            nargs='+',
                            required=True)
        parser.add_argument("-o", "--output",
                            help="ChIP-seq output filename prefix",
                            dest="output",
                            type=str,
                            default="rankprod",
                            required=False)
        parser.set_defaults(bigbed=False)
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
        parser.set_defaults(broadpeaks=False)
        parser.add_argument("--fragment",
                            help="Specifies whether the input peaks will be subject to high levels of fragmentation",
                            dest="fragment",
                            action="store_true",
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
                                 "Default: 0.5",
                            dest="random_seed",
                            default=0.5,
                            type=float,
                            required=False)
        parser.add_argument("-a","--alpha",
                            help="Alpha specifies the user cut-off value for set of reproducible peaks "
                                 "The analysis will still produce results including peaks within the threshold calculated"
                                 "using the binomial method \n"
                                 "Default: 0.05",
                            dest="alpha",
                            default=0.05,
                            type=float,
                            required=False)
        parser.add_argument("-s", "--size",
                            help="Sets the default minimum peak size when peaks are reconnected after fragmentation. \n"
                                 "Usually the minimum peak size is determined by the size of surrounding peaks, \n "
                                 "but in the case that there are no surrounding peaks this value will be used \n"
                                 "Default: 20",
                            dest="size",
                            default=20,
                            type=int,
                            required=False)
        args = parser.parse_args(sys.argv[1:])
        self.run_1(args)

    def run_1(self, args):
        for i in args.input:
            print(str(i.name))

        print('Processing Input...')

        bedfs = [bed.BedFile(str(i.name), 'Peaks') for i in args.input]

        rankprod.performrankprod(bedfs,
                                 args.minentries,
                                 args.rankmethod,
                                 'all',
                                 args.duphandling,
                                 args.random_seed,
                                 args.alpha,
                                 args.output,
                                 args.size,
                                 False,
                                 args.fragment)



def main():
    RankProdAnalysis()

if __name__ == "__main__":
    main()
