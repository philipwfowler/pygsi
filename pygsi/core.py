#! /usr/bin/env python

import numpy
import json
import requests
import pandas

class NucleotideStretch():

    """ The NucleotideStretch class accepts a nucleotide sequence as a string and then allows you to ask the SRA how many of
    each (amino-acid based) permutation there are.

    Args:
        nucleotide_sequence (str): nucleotide sequence containing only (a,t,c,g). Can be upper- or lower-case. Must be divisible by three.
        filename (str): allows previously saved NPY files to be loaded, thereby re-creating the instance of the class. Either this on nucleotide_sequence must be specified.
        gene_name (str): name of the gene. Metadata, not used in search. (optional)
        species_name (str) and species_min_amount (float): only results that are predicted by Bracken to have at least this amount (range 0-1) of the
            specified species will be included. Defaults are None and 0.80 (i.e. 80%). A species_name of None means all results are included.
        first_amino_acid_position (int): the number of the first amino acid encoded by nucleotide_sequence (default is 0)
    """

    def __init__(self,filename=None,nucleotide_sequence=None,gene_name=None,species_name=None,species_min_amount=0.80,first_amino_acid_position=0):

        # insist that either a filename or a nucleotide sequence is specified (and NOT both)
        assert (filename or nucleotide_sequence), "either a nucleotide sequence must be given as a string, or the filename of a .npy containing saved results"
        assert not (filename and nucleotide_sequence), "you cannot specify BOTH a nucleotide sequence and a save file"

        # if the user has given a nucleotide sequence, initialise all the various arrays etc
        if nucleotide_sequence:

            # how many nucleotides have we been given?
            self.number_nucleotides=len(nucleotide_sequence)

            # insist that the number of nucleotides must be divisible by three so can be decoded into amino acids
            assert (self.number_nucleotides%3)==0, "string of nucleotides must be exactly divisible by three"

            # how many amino acids does that therefore correspond to?
            self.number_amino_acids=int(self.number_nucleotides/3)

            # remember the position of the first and last amino acid and the name of gene
            self.first_amino_acid_position=first_amino_acid_position
            self.last_amino_acid_position=first_amino_acid_position+self.number_amino_acids
            self.gene_name=gene_name

            # store the target species and minimum amount (between 0-1)
            self.species_name=species_name
            self.species_min_amount=species_min_amount

            # store the nucleotides in lower case to differentiate them from amino acids
            nucleotides_string=nucleotide_sequence.lower()

            # store the nucleotides as a 1D numpy array
            self.nucleotides=numpy.array([n for n in nucleotides_string])

            # store the triplets that will encode the amino acids
            self.triplets=numpy.array([nucleotides_string[i:i+3] for i in range(0,len(nucleotides_string),3)])

            # create a list of codons as well as a dictionary to work out what amino acids are encoded
            bases = ['t', 'c', 'a', 'g']
            aminoacids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
            self.codons = numpy.array([a+b+c for a in bases for b in bases for c in bases])
            self.triplet_to_amino_acid = dict(zip(self.codons, aminoacids))
            self.amino_acids_of_codons=numpy.array([self.triplet_to_amino_acid[i] for i in self.codons])
            self.number_codons=len(self.amino_acids_of_codons)

            # now translate the triplets into amino acids using this new dictionary
            self.amino_acids=numpy.array([self.triplet_to_amino_acid[i] for i in self.triplets])

            # store the positions of the amino acids, as setup by first_amino_acid_position
            self.amino_acid_position=numpy.array([i for i in range(self.first_amino_acid_position,self.last_amino_acid_position,1)])

            # now interrogate the bigsi web instance to see how many genomes exist with the specified reference sequence
            # (this can take 5-10 seconds)
            self.number_reference_genomes=self._interrogate_bigsi(nucleotide_sequence)

            # instantiate the arrays dictionary that is going to hold all the numpy 2D arrays (this maybe could be done with pandas)
            self.arrays={}

            # make some boring 2D arrays
            self.arrays["original_amino_acid"]=(numpy.broadcast_to(self.amino_acids,(self.number_codons,self.number_amino_acids)))
            self.arrays["amino_acid_position"]=(numpy.broadcast_to(self.amino_acid_position,(self.number_codons,self.number_amino_acids)))
            self.arrays["new_amino_acid"]=(numpy.rot90(numpy.broadcast_to(self.amino_acids_of_codons,(self.number_amino_acids,self.number_codons)),3))
            self.arrays["number_genomes"]=numpy.zeros((self.number_codons,self.number_amino_acids),dtype=int)
            self.arrays["number_reference_genomes"]=self.number_reference_genomes*numpy.ones((self.number_codons,self.number_amino_acids),dtype=int)

            # triplets here means the triplet encoding the original amino acid, whilst codons is all 64 possibilities
            self.arrays["original_triplet"]=numpy.broadcast_to(self.triplets,(self.number_codons,self.number_amino_acids))
            self.arrays["new_triplet"]=numpy.rot90(numpy.broadcast_to(self.codons,(self.number_amino_acids,self.number_codons)),3)

            # concatenate the original_amino_acid, position and new_amino_acid to make a mutation label
            tmp=numpy.core.defchararray.add(self.arrays["original_amino_acid"],self.arrays["amino_acid_position"].astype(str))
            foo=numpy.core.defchararray.add(tmp,self.arrays["new_amino_acid"])
            self.arrays["mutation"]=numpy.where(self.arrays["original_triplet"]!=self.arrays["new_triplet"],foo,"-")


            # construct simple Boolean arrays
            self.arrays["synonymous"]=(self.arrays["original_amino_acid"]==self.arrays["new_amino_acid"])
            self.arrays["non_synonymous"]=(self.arrays["original_amino_acid"]!=self.arrays["new_amino_acid"])

            # finally construct an array that tells you how many SNPs you need to do (1,2,3) to reach the new triplet
            self.arrays["number_nucleotide_changes"]=numpy.zeros((self.number_codons,self.number_amino_acids))
            for col in range(self.number_amino_acids):
                for row in range(self.number_codons):
                    self.arrays["number_nucleotide_changes"][(row,col)]=sum(1 for a,b in zip(self.codons[row],self.triplets[col]) if a!=b)

        # if a filename is given, simply load all the variables and arrays from the specified NPY file.
        elif filename:

            data=numpy.load(filename)
            self.nucleotides=data[()]["nucleotides"]
            self.amino_acids=data[()]["amino_acids"]
            self.triplets=data[()]["triplets"]
            self.first_amino_acid_position=data[()]["first_amino_acid_position"]
            self.amino_acid_position=data[()]["amino_acid_position"]
            self.last_amino_acid_position=data[()]["last_amino_acid_position"]
            self.number_amino_acids=data[()]["number_amino_acids"]
            self.number_nucleotides=data[()]["number_nucleotides"]
            self.gene_name=data[()]["gene_name"]
            self.species_name=data[()]["species_name"]
            self.species_min_amount=data[()]["species_min_amount"]
            self.codons=data[()]["codons"]
            self.triplet_to_amino_acid=data[()]["triplet_to_amino_acid"]
            self.amino_acids_of_codons=data[()]["amino_acids_of_codons"]
            self.number_codons=data[()]["number_codons"]
            self.number_reference_genomes=data[()]["number_reference_genomes"]
            self.arrays={}
            self.arrays["original_amino_acid"]=data[()]["arr_original_amino_acid"]
            self.arrays["new_amino_acid"]=data[()]["arr_new_amino_acid"]
            self.arrays["amino_acid_position"]=data[()]["arr_amino_acid_position"]
            self.arrays["mutation"]=data[()]["arr_mutation"]
            self.arrays["synonymous"]=data[()]["arr_synonymous"]
            self.arrays["non_synonymous"]=data[()]["arr_non_synonymous"]
            self.arrays["original_triplet"]=data[()]["arr_original_triplet"]
            self.arrays["new_triplet"]=data[()]["arr_new_triplet"]
            self.arrays["number_genomes"]=data[()]["arr_number_genomes"]
            self.arrays["number_reference_genomes"]=data[()]["arr_number_reference_genomes"]
            self.arrays["number_nucleotide_changes"]=data[()]["arr_number_nucleotide_changes"]

            # now that all the arrays are loaded, re-create the Pandas dataframe (bit hacky)
            self._create_dataframe()


    def save(self,filename):
        """ Save all the variables and arrays belonging to the class. This allows the instance to be re-created
        by using the load() method"""

        data={}
        data['nucleotides']=self.nucleotides
        data['amino_acids']=self.amino_acids
        data['triplets']=self.triplets
        data['first_amino_acid_position']=self.first_amino_acid_position
        data['amino_acid_position']=self.amino_acid_position
        data['last_amino_acid_position']=self.last_amino_acid_position
        data['number_amino_acids']=self.number_amino_acids
        data['number_nucleotides']=self.number_nucleotides
        data['gene_name']=self.gene_name
        data['species_name']=self.species_name
        data['species_min_amount']=self.species_min_amount
        data['codons']=self.codons
        data['triplet_to_amino_acid']=self.triplet_to_amino_acid
        data['amino_acids_of_codons']=self.amino_acids_of_codons
        data['number_codons']=self.number_codons
        data['number_reference_genomes']=self.number_reference_genomes
        for i in self.arrays.keys():
            data["arr_"+i]=self.arrays[i]

        # save the dictionary to the specified file
        numpy.save(filename, data)

    def __repr__(self):
        """ Change so that the print() function outputs a summary of the instance.
        """

        amino_acids_string=''.join(self.amino_acids)
        nucleotides_string=''.join(self.nucleotides)

        if self.species_name is None:
            species="All species considered"
        else:
            species=self.species_name

        return("%s\n%s gene\n%s to %s\n%s\n%s\n%i genomes found with this sequence" % (species,self.gene_name,amino_acids_string[0]+str(self.first_amino_acid_position),amino_acids_string[-1]+str(self.first_amino_acid_position+self.number_amino_acids),amino_acids_string,nucleotides_string,self.number_reference_genomes))

    def __len__(self):
        """ Change so len() simply returns the number of nucleotides stored within the instance of the class.
        """
        return(len(self.nucleotides))

    def permuate_position(self,aminoacid_number,triplet_position=False,output=False):
        '''For a given triplet encoding an amino acid in the sequence, count how many occurences of
        different permutations there are in the SRA.

        This method is written such that if the triplet_position is specified, the four permutations
        at that position will be calculated, but if triplet_position is NOT specified, all 64 codons
        will be tried.

        Args:
            triplet_position (int): which nucleotide to permute (1,2,3)
        '''

        # check that, if specified, the triplet_position can be only 1,2 or 3
        if triplet_position:
            assert triplet_position in (1,2,3), "triplet_position can only be one of {1,2,3}"

        # find out the 0-based index of the amino acid position being permutated
        triplet_idx=aminoacid_number-self.first_amino_acid_position

        # remember the reference triplet and corresponding amino acid
        original_triplet=self.triplets[triplet_idx]
        original_aminoacid=self.amino_acids[triplet_idx]

        if triplet_position:

            # from that, work out the 0-based nuclotide position
            nucleotide_idx=(3*(triplet_idx))+triplet_position-1

            nucleotides_string=''.join(self.nucleotides)

            # split the nucleotide sequence
            nucleotides_before=nucleotides_string[:nucleotide_idx]
            nucleotide_original=nucleotides_string[nucleotide_idx]
            nucleotides_after=nucleotides_string[nucleotide_idx+1:]

            permutations=['a','c','t','g']

            # don't consider the original nucleotide
            permutations.remove(nucleotide_original)

        else:

            # from the 0-based triplet index, work out the 0-based nuclotide index
            nucleotide_idx=(3*(triplet_idx))

            # split the nucleotide sequence
            nucleotides_before=''.join(self.nucleotides[:nucleotide_idx])
            tmp=self.nucleotides[nucleotide_idx:nucleotide_idx+3]
            triplet_original=''.join(i for i in tmp)
            nucleotides_after=''.join(self.nucleotides[nucleotide_idx+3:])

            # codons is a numpy array, so convert to a list so that remove will work
            permutations=self.codons.tolist()

            # remove the original triplet so only 63 permutations will be considered
            permutations.remove(triplet_original)

        first_pass=True

        for new_sequence in permutations:

            # create the query sequence
            query_sequence=nucleotides_before+new_sequence+nucleotides_after

            # find out the resulting new triplet
            new_triplet=[query_sequence[i:i+3] for i in range(0,len(query_sequence),3)][triplet_idx]

            # ..and therefore the new amino acid
            new_aminoacid=self.triplet_to_amino_acid[new_triplet]

            # count how many occurences are in the SRA by calling the BIGSI instance
            total=self._interrogate_bigsi(query_sequence)

            # work out which row this codon belongs to (we already know that the col==triplet_idx)
            row=numpy.where(self.codons==new_triplet)[0][0]

            # store the number of sequences in the final 2D array
            self.arrays["number_genomes"][(row,triplet_idx)]=total

            if output:
                if first_pass:
                    print("%1s%i encoded by triplet %3s." % (original_aminoacid,aminoacid_number,original_triplet))
                first_pass=False
                print("%3s : %7i " % (new_triplet,total))

        # finally, store all the results in a Pandas DataFrame
        self._create_dataframe()

    def _interrogate_bigsi(self,nucleotides_string):
        """ This is a private method that, given a string of nucleotides, calls the current BIGSI test site
        and parses the results, given the target species_name and species_min_amount.

        Args:
            nucleotides_string (str): nucleotide sequence containing only (a,t,c,g) to be queried. Can be upper- or lower-case. Must be divisible by three.
        """

        # create an empty list to store the list of sample numbers from the Short Read Archive
        sra_samples=[]

        # define the URL of the BIGSI instance
        url_front="http://api.bigsi.io/search?seq="
        url_end="&threshold=1.0"

        query_string=nucleotides_string.upper()

        # call the Web API
        r=requests.get(url_front+query_string+url_end)

        # parse the returned data
        result=json.loads(r.text)

        # loop through the samples
        for sra_sample_number in result[query_string]['results']:

            # if no target species is specified, simply count 'em
            if self.species_name is None:
                sra_samples.append(sra_sample_number)

            else:
                # pull out the list of predicted species, and their predicted amounts
                predicted_species_list=result[query_string]['results'][sra_sample_number]["species"].split("; ")

                for line in predicted_species_list:

                    tmp=line.split(" : ")

                    if len(tmp)>1:

                        predicted_species_name=tmp[0]
                        predicted_species_amount=float(tmp[1].split("%")[0])

                        # only keep the sample if it matches the target species and is above the threshold amount
                        if predicted_species_name==self.species_name:
                            if (predicted_species_amount/100.)>=self.species_min_amount:
                                sra_samples.append(sra_sample_number)

        return(len(sra_samples))

    def _create_dataframe(self):
        """Create a Pandas dataframe that lists all the mutations and their occurences, one per line"""

        data_dict={}

        # create a Boolean array of only those positions where sequences have been identified
        positive_elements=self.arrays["number_genomes"]>0

        for key in ['amino_acid_position','original_triplet','new_triplet','number_nucleotide_changes','mutation','number_genomes','number_reference_genomes','original_amino_acid','new_amino_acid','synonymous','non_synonymous']:
            data_dict[key]=(self.arrays[key][positive_elements]).tolist()

        self.df=pandas.DataFrame(data=data_dict)

        self.df["number_nucleotide_changes"]=self.df["number_nucleotide_changes"].astype("int8")
