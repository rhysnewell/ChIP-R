ChIP-R ("chipper")
=========

ChIP-R uses an adaptation of the rank product statistic to assess the reproducibility of ChIP-seq peaks by incorporating information from multiple ChIP-seq replicates and "fragmenting" peak locations to better combine the information present across the replicates.

Install
-------

- [Python3.x](https://www.python.org/getit/) with the following packages:
- Numpy
- Scipy

To install ChIP-R:

  git clone git@github.com:rhysnewell/ChIP-R.git
  
We recommend you add the src folder to your pythonpath variable in your .bashrc profile


Usage
-----

In the command line, type in **'python3 ChIP-R.py -h '** for detailed usage.

    $ python3 ChIP-R.py -h
    
    usage: ChIP-R [-h] -i INPUT [INPUT ...] [-o OUTPUT] [-m MINENTRIES]
                  [--rankmethod RANKMETHOD] [--duphandling DUPHANDLING]
                  [--seed RANDOM_SEED] [-a ALPHA]

    Combine multiple ChIP-seq files and return a union of all peak locations and a
    set confident, reproducible peaks as determined by rank product analysis

    optional arguments:
      -h, --help            show this help message and exit
      -i INPUT [INPUT ...], --input INPUT [INPUT ...]
                            ChIP-seq input files. These files must be in either
                            narrowPeak, broadPeak, or regionPeak format. Multiple
                            inputs a separeted by a single space
      -o OUTPUT, --output OUTPUT
                            ChIP-seq output filename
      -m MINENTRIES, --minentries MINENTRIES
                            The minimum peaks between replicates required to form
                            an intersection of the peaks Default: 1
      --rankmethod RANKMETHOD
                            The ranking method used to rank peaks within
                            replicates. Options: 'signalvalue', 'pvalue',
                            'qvalue'. Default: pvalue
      --duphandling DUPHANDLING
                            Specifies how to handle entries that are ranked
                            equally within a replicate Can either take the
                            'average' ranks or a 'random' rearrangement of the
                            ordinal ranks Options: 'average', 'random' Default:
                            'average'
      --seed RANDOM_SEED    Specify a seed to be used in conjunction with the
                            'random' option for -duphandling Must be between 0 and
                            1 default: 0.5
      -a ALPHA, --alpha ALPHA
                            Alpha specifies the user cut-off value for set of
                            reproducible peaks The analysis will still produce
                            results including peaks within the threshold
                            calculatedusing the binomial method default: 0.05


Example
------
    $ python3 ChIP-R.py -i input_prefix1.bed input_prefix2.bed input_prefix3.bed input_prefix4.bed -m 2 -o output_prefix   

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
|chr1|9118 |10409|T3_peak_87823|	491|	.	|15.000000	| 0.113938|0.712353	|


Citation
--------




Contact
-------

Authors: Rhys Newell, Alexandra Essebier, Michael Piper, Mikael Boden

Contact:  rhys.newell(AT)uq.net.au
