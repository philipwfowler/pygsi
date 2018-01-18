#! /usr/bin/env python

import pygsi
from tqdm import tqdm

# the oxa1 gene as specified on http://bigsi.io
oxa1_sequence="ATGAAAAACACAATACATATCAACTTCGCTATTTTTTTAATAATTGCAAATATTATCTACAGCAGCGCCAGTGCATCAACAGATATCTCTACTGTTGCATCTCCATTATTTGAAGGAACTGAAGGTTGTTTTTTACTTTACGATGCATCCACAAACGCTGAAATTGCTCAATTCAATAAAGCAAAGTGTGCAACGCAAATGGCACCAGATTCAACTTTCAAGATCGCATTATCACTTATGGCATTTGATGCGGAAATAATAGATCAGAAAACCATATTCAAATGGGATAAAACCCCCAAAGGAATGGAGATCTGGAACAGCAATCATACACCAAAGACGTGGATGCAATTTTCTGTTGTTTGGGTTTCGCAAGAAATAACCCAAAAAATTGGATTAAATAAAATCAAGAATTATCTCAAAGATTTTGATTATGGAAATCAAGACTTCTCTGGAGATAAAGAAAGAAACAACGGATTAACAGAAGCATGGCTCGAAAGTAGCTTAAAAATTTCACCAGAAGAACAAATTCAATTCCTGCGTAAAATTATTAATCACAATCTCCCAGTTAAAAACTCAGCCATAGAAAACACCATAGAGAACATGTATCTACAAGATCTGGATAATAGTACAAAACTGTATGGGAAAACTGGTGCAGGATTCACAGCAAATAGAACCTTACAAAACGGATGGTTTGAAGGGTTTATTATAAGCAAATCAGGACATAAATATGTTTTTGTGTCCGCACTTACAGGAAACTTGGGGTCGAATTTAACATCAAGCATAAAAGCCAAGAAAAATGCGATCACCATTCTAAACACACTAAATTTA"

first_nucleotide_position=0
last_nucleotide_position=len(oxa1_sequence)

print("OXA1 is %i bases long" % len(oxa1_sequence))

# create an instance of the class by giving it the nucleotide sequence as a string
# by setting species=None, all results will be considered
oxa1=pygsi.NucleotideStretch(nucleotide_sequence=oxa1_sequence,\
                             gene_name="oxa1",\
                             first_amino_acid_position=1,\
                             species_name=None)

# let's have a look at the summary
print(oxa1)

for position in tqdm(range(228,245)):
    for i in [1,2,3]:
        oxa1.permuate_position(position,triplet_position=i)

# view the dataframe
print(oxa1.df)

# save all the variables and arrays to a NPY file
oxa1.save("dat/oxa1.npy")

# save the Pandas dataset to a CSV and DTA file
oxa1.df.to_csv("dat/oxa1.csv")
oxa1.df.to_stata("dat/oxa1-oxa1.dta")
