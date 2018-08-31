# pygsi: a Python class to interrogate BIGISI 

As discussed in the preprint below, whilst the Short Read Archive (SRA) holds nearly half a million bacterial and virial genomes, it is extremely difficult to, given a nucleotide sequence X, answer the question "how many times has X been deposited in the SRA?".

> Bradley P, Den Bakker HC, Rocha EPC, McVean G, Iqbal Z. Real-time search of all bacterial and viral genomic data. 2017. [biorXiv](https://dx.doi.org/10.1101/234955)

This citation will be updated when this paper is published.

As part of their paper, the authors have made available [a web portal](https://bigsi.io), allowing individual nucleotide strings to be entered, with the results returned in a few seconds. This Python3 class allows this website to be interrogated programmatically, allowing one to systematically study the universe of small variations of the original nucleotide sequence X. Specifically, through the `permuate_positions()` method of the `NucleotideSequence` class, sequences with the three alternative nucleotides substituted can be queried against the BIGSI instance, allowing for example, all sequences with one SNP to be identified. Alternatively, the other 63 triplets from the codon table can be tried at each position. 

## Citation

If you use the pygsi package, please cite using the DOI below (in addition to the BIGSI paper referenced above).

[![DOI](https://zenodo.org/badge/117247246.svg)](https://zenodo.org/badge/latestdoi/117247246)

## Approach

It is important to know at a high level how the BIGSI index is interrogated as then you will understand some of the limitations of the approach. We pass to BIGSI the triplet we wish to examine for variation flanked by 30 bases on either side (or equivalently, a single amino acid flanked by 10 residues). Conceptually this can be represented as 

     ------------------------------ooo------------------------------

where '-' is a flanking base and 'o' is a base we are going to permute. An immediate consequence of this is that if you wish to examine the variation at the first residue in the protein sequence, then the nucleotide sequence must start 30 bases before the first residue (and end 30 bases afterwards). If you do not, then a shorter kmer will be passed to BIGSI (which has a minimum kmer length of 61 bases) and many false positives will be reported.  

Consider a trivial case where we wish to find out all the single amino acid (i.e. single triplet) differences in a stretch of 12 amino acids i.e. (amino acid pos 1-12 incl, or bases 1-36). We therefore need to provide pygsi with amino acids -9 to 22 incl. (bases -29 to 66) and so the first fishing sequence is

     base#                                  1         2         3
                                   123456789012345678901234567890123456
     ------------------------------ooo------------------------------................................

where '.' is a nucleotide that is not included in the fishing sequence defined by '-' and 'o' and we can find all variation using the permute_position() method of the pygsi class. Then we move on to consider the second amino acid position

     base#                                  1         2         3
                                   123456789012345678901234567890123456
     ...------------------------------ooo------------------------------..............................

and then the third

     base#                                  1         2         3
                                   123456789012345678901234567890123456
     ......------------------------------ooo------------------------------...........................

all the way up to the twelth

     base#                                  1         2         3
                                   123456789012345678901234567890123456
     .................................------------------------------ooo------------------------------

The other important consequence of fishing with a 63-kmer is that the approach is blind to double mutations that are both covered by the 63-kmer i.e. separated by twenty one or fewer amino acids, hence it cannot detect any double mutations in our test sequence as it is too short.


## Examples

For more information, please look at the included simple example that look for single nucleotide variants of the reference sequence of the important antibiotic resistance gene OXA-1. Simple commented python and Jupyter Notebook forms are included. To run the former:

     $ cd examples/ 
     $ python analyse-oxa1-variants.py 

Whilst to run the latter you need to have jupyter installed

     $ jupyter-notebook analyse-oxa1-variants.ipynb
     
and a window will be opened in your browser with the code. All the output files are written to examples/dat. It is easiest to look at the output by loading the CSV file into a spreadsheet package like Microsoft Excel or Apple Numbers.

## Installation

Download or clone the repository

     $ git clone https://github.com/philipwfowler/pygsi.git
     $ cd pygsi/

If you only wish to install the package in your $HOME directory (or don't have sudo access) issue the --user flag

     $ python setup.py install --user
     
Alternatively, to install system-wide

     $ sudo python setup.py install

## Pre-requisites

The `setup.py` script will automatically install the following Python modules. At present `setup.py` specifies minimum versions loosely based on those below. It is likely to work with earlier versions, but please exercise caution and you will need to edit the `setup.py` file to remove the version restriction. The code is `Python3` and was written and tested using version 3.5 on a Mac.

- `requests` (2.8.14)
- `numpy` (1.13.3)
- `pandas` (0.21.0)
- `tqdm` (4.19.4, not strictly part of the pygsi module but is used in the included OXA-1 example and besides, I like it.)



