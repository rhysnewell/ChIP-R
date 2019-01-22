Rank Product Analysis
=========


Install
-------

- [Python3.x](https://www.python.org/getit/) with the following packages:
- Numpy
- Scipy

To install ChIA-PET2:

  git clone git@github.com:rhysnewell/ChIP-rep.git
  
We recommend you add the src folder to your pythonpath variable in your .bashrc profile


Usage
-----

In the command line, type in **' REP.py chiprep -h '** for detailed usage.

    $ REP.py chiprep -h
    
    usage: chiprep [-h] -input INPUT [INPUT ...] [-output OUTPUT]
                 [-minentries MINENTRIES] [-rankmethod RANKMETHOD]
                 [-duphandling DUPHANDLING] [-randomseed RANDOM_SEED]
                 [-alpha ALPHA]

    Combine multiple ChIP-seq files and return an intersection of all peak locations and a
    set confident, reproducible peaks as determined by rank product analysis
    
    optional arguments:
      -h, --help            show this help message and exit
      -input INPUT [INPUT ...]
                            ChIP-seq input files. These files must be in either
                            narrowPeak, broadPeak, or regionPeak format
      -output OUTPUT        ChIP-seq output filename
      -minentries MINENTRIES
                            The minimum peaks between replicates required to form
                            an intersection of the peaks Default: 2
      -rankmethod RANKMETHOD
                            The ranking method used to rank peaks within
                            replicates. Options: 'signalvalue', 'pvalue',
                            'qvalue'. Default: signalvalue
      -duphandling DUPHANDLING
                            Specifies how to handle entries that are ranked
                            equally within a replicate Can either take the
                            'average' ranks or a 'random' rearrangement of the
                            ordinal ranks Options: 'average', 'random' Default:
                            'average'
      -randomseed RANDOM_SEED
                            Specify a seed to be used in conjunction with the
                            'random' option for -duphandling Must be between 0 and
                            1 default: 0.5
      -alpha ALPHA          Alpha specifies the user cut-off value for set of
                            reproducible peaks The analysis will still produce
                            results including peaks within the threshold
                            calculatedusing the binomial method default: 0.05


    Performs rank product analysis on ChIP data to either verify reproducible
    ChIP-seq peaks or ChIA-PET interactions
    
    positional arguments:
      command     Subcommand to run
    
    optional arguments:
      -h, --help  show this help message and exit

Example
------
    $ REP.py chiprep -i input_prefix1.bed input_prefix2.bed input_prefix3.bed input_prefix4.bed -minentries 1 -o output_prefix   

Output
------

Important result files:

- **prefixname_ALL.bed**: All intersected peaks, ordered from most significant to least (10 columns)
- **prefixname_T2.bed**: The tier 2 intersected peaks, the peaks that fall within the binomial threshold (10 columns)
- **prefixname_T1.bed**: The tier 1 intersected peaks, the peaks that fall within the user defined threshold (10 columns)
- **prefixname_log.txt**: A log containing the number of peaks appearing in each tier.


prefixname.bed file has 10 columns. The output follows the standard peak format for bed files, with the addition of a 10th column that specifies the ranks of the peaks that produced this possible peak. See the toy example below.

|chr |start|end  |name |score |strand  |signalValue |p-value |q-value|
|----|-----|-----|----|------|-----|------|------|------|
|chr1|9118 |10409|T3_peak_87823|	491|	.	|1.000000	| 0.113938|0.712353	|


Citation
--------




Contact
-------

Author: Rhys Newell

Email:  rhys.newell(AT)uq.net.au
