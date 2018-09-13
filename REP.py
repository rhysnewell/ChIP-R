#!/usr/bin/python3

import argparse
import sys
import RankProd
import bed

class RankProdAnalysis(object):

    def __init__(self):
        parser = argparse.ArgumentParser(
        description="Performs rank product analysis on ChIP data to either verify"
                                                     " reproducible ChIP-seq peaks or ChIA-PET interactions",
        usage='''REP <command> [<args>]
        chiprep    Rank Product reproducibility analysis for ChIP-seq
        chiarep    Rank Product reproducibility analysis for ChIA-PET
        ''')
        parser.add_argument('command', help='Subcommand to run')
        args = parser.parse_args(sys.argv[1:2])
        if not hasattr(self, args.command):
            print('Unrecognized command')
            parser.print_help()
            exit(1)
        getattr(self, args.command)()

    def run_1(self, args):
        for i in args.input:
            print(str(i.name))

        print('Processing BedFiles...')
        bedfs = [bed.BedFile(str(i.name), 'IDR') for i in args.input]

        RankProd.performrankprod(bedfs,
                                 args.minentries,
                                 args.rankmethod,
                                 'all',
                                 args.duphandling,
                                 args.random_seed,
                                 args.alpha,
                                 args.output)

    def run_2(self, args):
        RankProd.ChIAreproducibility(args.input1,
                                     args.input2,
                                     args.minentries,
                                     args.rankmethod,
                                     args.RPthreshold,
                                     args.alpha,
                                     args.output)


    def chiprep(self):
        parser = argparse.ArgumentParser(prog='chiprep',
                                         description="Combine multiple ChIP-seq files and return a union of all peak "
                                                     "locations and a set confident, reproducible peaks as determined by "
                                                     "rank product analysis")
        parser.add_argument("-input",
                            help="ChIP-seq input files. These files must be in either narrowPeak, broadPeak, "
                                 "or regionPeak format",
                            dest="input",
                            type=argparse.FileType('r'),
                            nargs='+',
                            required=True)
        parser.add_argument("-output",
                            help="ChIP-seq output filename",
                            dest="output",
                            type=str,
                            default="rankprod",
                            required=False)
        parser.add_argument("-minentries",
                            help="The minimum peaks between replicates required to form an "
                                 "intersection of the peaks \n"
                                 "Default: 1",
                            dest="minentries",
                            default=1,
                            type=int,
                            required=False)
        parser.add_argument("-rankmethod",
                            help="The ranking method used to rank peaks within replicates. "
                                 "Options: 'signalvalue', 'pvalue', 'qvalue'. \n"
                                 "Default: signalvalue",
                            dest="rankmethod",
                            default='pvalue',
                            type=str,
                            required=False)
        parser.add_argument("-duphandling",
                            help="Specifies how to handle entries that are ranked equally within a replicate "
                                 "Can either take the 'average' ranks or a 'random' rearrangement of the ordinal ranks \n"
                                 "Options: 'average', 'random' \n"
                                 "Default: 'average'",
                            dest="duphandling",
                            default='average',
                            type=str,
                            required=False)
        parser.add_argument("-randomseed",
                            help="Specify a seed to be used in conjunction with the 'random' option for -duphandling "
                                 "Must be between 0 and 1 \n"
                                 "default: 0.5",
                            dest="random_seed",
                            default=0.5,
                            type=float,
                            required=False)
        parser.add_argument("-alpha",
                            help="Alpha specifies the user cut-off value for set of reproducible peaks "
                                 "The analysis will still produce results including peaks within the threshold calculated"
                                 "using the binomial method \n"
                                 "default: 0.05",
                            dest="alpha",
                            default=0.05,
                            type=float,
                            required=False)
        args = parser.parse_args(sys.argv[2:])
        self.run_1(args)

    def chiarep(self):
        parser = argparse.ArgumentParser(prog='chiarep',
                                         description="Combine multiple ChIA-seq files and return a union of all peak "
                                                     "locations and a set confident, reproducible peaks as determined by "
                                                     "rank product analysis")
        parser.add_argument("-input1",
                            help="ChIA-PET input files.",
                            dest="input",
                            type=argparse.FileType('r'),
                            nargs='+',
                            required=True)
        parser.add_argument("-input2",
                            help="ChIP-seq input files. These files must be in either narrowPeak, broadPeak, "
                                 "or regionPeak format",
                            dest="input",
                            type=argparse.FileType('r'),
                            nargs='+',
                            required=True)
        parser.add_argument("-output",
                            help="ChIA-PET output filename",
                            dest="output",
                            type=str,
                            required=False)
        parser.add_argument("-minentries",
                            help="The minimum peaks between replicates required to form an "
                                 "intersection of the peaks \n"
                                 "Default: 1",
                            dest="minentries",
                            type=int,
                            required=False)
        parser.add_argument("-rankmethod",
                            help="The ranking method used to rank peaks within replicates. "
                                 "Options: 'signalvalue', 'pvalue', 'qvalue'. \n"
                                 "Default: signalvalue",
                            dest="rank",
                            type=str,
                            required=False)
        parser.add_argument("-RPthreshold",
                            help="Specifies which peak threshold group to use to verify the ChIA-PET interactions "
                                 "\n"
                                 "Options: 'alpha', 'binom', 'all' \n"
                                 "Default: 'all'",
                            dest="threshold",
                            type=str,
                            required=False)
        parser.add_argument("-alpha",
                            help="Alpha specifies the user cut-off value for set of reproducible peaks and interactions"
                                 "\n"
                                 "default: 0.05",
                            dest="alpha",
                            type=float,
                            required=False)
        args = parser.parse_args(sys.argv[2:])
        self.run_2(args)


if __name__ == "__main__":
    RankProdAnalysis()
