# pygsi: a Python class to interrogate BIGISI 

As discussed in the preprint below, whilst the Short Read Archive (SRA) holds nearly half a million bacterial and virial genomes, it is extremely difficult to, given a nucleotide sequence X, answer the question "how many times has X been deposited in the SRA?".

> Bradley P, Den Bakker HC, Rocha EPC, McVean G, Iqbal Z. Real-time search of all bacterial and viral genomic data. 2017. [biorXiv](https://dx.doi.org/10.1101/234955)

This citation will be updated when this paper is published.

As part of their paper, the authors have made available [a web portal](https://bigsi.io), allowing individual nucleotide strings to be entered, with the results returned in a few seconds. This Python3 class allows this website to be interrogated programmatically, allowing one to systematically study the universe of small variations of the original nucleotide sequence X. Specifically, through the permuate_positions() method, sequences with the three alternative nucleotides substituted can be queried against the BIGSI instance, allowing for example, all sequences with one SNP to be identified. Alternatively, the other 63 triplets from the codon table can be tried at each position.

## Examples

For more information, please look at the included simple example that look for single nucleotide variants of the reference sequence of the important antibiotic resistance gene OXA-1. Simple python and Jupyter Notebook forms are included. To run the former:

     $ cd examples/ 
     $ python analyse-oxa1-variants.py 

Whilst to run the latter you need to have jupyter installed

     $ jupyter-notebook analyse-oxa1-variants.ipynb
     
and a window will be opened in your browser with the code.     

## Installation

Download or clone the repository

     $ git clone https://github.com/philipwfowler/pygsi.git
     $ cd pygsi/

If you only wish to install the package in your $HOME directory (or don't have sudo access) issue the --user flag

     $ python setup.py install --user
     
Alternatively, to install system-wide

     $ sudo python setup.py install

## Pre-requisites

The setup.py script will automatically install the following Python modules. At present setup.py specifies minimum versions loosely based on those below. It is likely to work with earlier versions, but please exercise caution and you will need to edit the setup.py file to remove the version restriction. The code is Python3 and was written and tested using version 3.5 on a Mac.

- json (2.09)
- requests (2.8.14)
- numpy (1.13.3)
- pandas (0.21.0)
- tqdm (4.19.4, optional i.e. not part of the pygsi module but is used in the included OXA-1 example and besides, I like it.)



