ChIP-R ("chipper")
==================

ChIP-R uses an adaptation of the rank product statistic to assess the reproducibility of ChIP-seq peaks by incorporating information from multiple ChIP-seq replicates and "fragmenting" peak locations to better combine the information present across the replicates.

Install
-------

- [Python3.x](https://www.python.org/getit/) with the following packages:
- Numpy
- Scipy

To install ChIP-R:
    
    pip install ChIP-R
    
OR if you want to install from source:

    git clone https://github.com/rhysnewell/ChIP-R.git
    cd ChIP-R
    python3 setup.py install



Usage
-----

ChIP-R requires only a single input type: A set of any number of BED file regions. Typically the output of peak calling from 
ChIP-seq peak calling on transcription factor or histone mark samples. Alternatively, ChIP-R can also be used on 
ATAC-seq peaks to retrieve reproducible peaks across ATAC-seq experiments.


#### Input

The input BED files must follow ENCODE narrowPeak or broadPeak format specifications. Typically, this format is the default
for peak callers such as MACS2. 

#### Peak calling

ChIP-R is compatible with the output peaks for any peak caller as long as the output is in the correct narrowPeak or broadPeak
format. Additionally, there is no need to call peaks with relaxed thresholds when using your chosen peak caller as is the suggested
by IDR.

#### Parameters

ChIP-R is fairly light on parameters that need to be chosen by the user. A couple of options that users may want to play with is
`minentries` and `size`. 

`minentries` determines the number of peak overlaps required to start calling a peak "reproducible". 
The default of 2 typically provides the best results in our benchmarks but there may be a case where a user requires 
ChIP-R to call peaks within a much stricter window.

`size` determines the minimum peak size during peak output. Transcription factors generally want more punctate peaks, and 
so the default value of 20 may be sufficient. However, histone marks may require a much larger value be set for this depending
on how broad you expect the histone mark to be. Generally, if you find ChIP-R produces too many small noisy peaks then this 
value can be increased to filter them out.

Example
------
    $ chipr -i sample1.bed sample2.bed sample3.bed sample4.bed -m 2 -o output_prefix   

In the command line, type in **'chipr -h '** for detailed usage.

    $ chipr -h
    
    usage: chipr [-h] -i INPUT [INPUT ...] [-o OUTPUT] [-m MINENTRIES]
             [--rankmethod RANKMETHOD] [--duphandling DUPHANDLING]
             [--seed RANDOM_SEED] [-a ALPHA]

    Combine multiple ChIP-seq files and return a union of all peak locations and a
    set confident, reproducible peaks as determined by rank product analysis

    optional arguments:
      -h, --help            show this help message and exit
      -i INPUT [INPUT ...], --input INPUT [INPUT ...]
                            ChIP-seq input files. These files must be in either
                            narrowPeak, broadPeak, or regionPeak format. Multiple
                            inputs are separeted by a single space
      -o OUTPUT, --output OUTPUT
                            ChIP-seq output filename prefix
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
                            1 Default: 0.5
      -a ALPHA, --alpha ALPHA
                            Alpha specifies the user cut-off value for set of
                            reproducible peaks The analysis will still produce
                            results including peaks within the threshold
                            calculated using the binomial method Default: 0.05
      -s SIZE, --size SIZE  Sets the default minimum peak size when peaks are
                            reconnected after fragmentation. Usually the minimum
                            peak size is determined by the size of surrounding
                            peaks, but in the case that there are no surrounding
                            peaks this value will be used Default: 20






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

**Preprint available on bioarxiv**
https://www.biorxiv.org/content/10.1101/2020.11.24.396960v1



Contact
-------

Authors: Rhys Newell, Michael Piper, Mikael Boden, Alexandra Essebier

Contact:  rhys.newell(AT)hdr.qut.edu.au
