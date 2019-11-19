#! /usr/bin/env python

import numpy
import json
import requests
import pandas

class NucleotideStretch():

    """ The NucleotideStretch class accepts a nucleotide sequence as a string and then allows you to ask the SRA how many of
    each (amino-acid based) permutation there are using the developmental BIGSI webserver.

    Args:
        nucleotide_sequence (str): nucleotide sequence containing only (a,t,c,g). Can be upper- or lower-case. Must be divisible by three.
        filename (str): allows previously saved NPY files to be loaded, thereby re-creating the instance of the class. Either this on nucleotide_sequence must be specified.
        gene_name (str): name of the gene. Metadata, not used in search. (optional)
        species_name (str) and species_min_amount (float): only results that are predicted by Bracken to have at least this amount (range 0-1) of the
            specified species will be included. Defaults are None and 0.80 (i.e. 80%). A species_name of None means all results are included.
        first_amino_acid_position (int): the number of the first amino acid encoded by nucleotide_sequence (default is 1)

    Important:
        if you wish to find all the sequences with one amino acid difference, you must either (i) provide the nucleotide sequence for amino acids -9 to N+10 or (ii) if you
        only give the nucleotide sequence for the CDS, then you are restricted to examining 11<residues<N-10. This is because the code iteratively examines the variation in a single triplet,
        flanked by 30 bases on either side. If the 30 flanking bases are not present (for example if you try interrogating position 1 but haven't provided the sequence down to the equivalent
        amino acid position of -9), then BIGSI will return false positives. The other implication of using a 63-kmer is that you are implicitly assuming that all mutations are greather
        than 30 bases apart, which may not be the case i.e. unless you explicitly look for them, the approach is blind to double mutations separated by 10 or fewer amino acids.
    """

    def __init__(self,filename=None,nucleotide_sequence=None,gene_name=None,species_name=None,species_min_amount=0.80,first_amino_acid_position=1):

        # insist that either a filename or a nucleotide sequence is specified (and NOT both)
        assert (filename or nucleotide_sequence), "either a nucleotide sequence must be given as a string, or the filename of a .npy containing saved results"
        assert not (filename and nucleotide_sequence), "you cannot specify BOTH a nucleotide sequence and a save file"

        # if the user has given a nucleotide sequence, initialise all the various arrays etc
        if nucleotide_sequence:

            self.mutations=pandas.DataFrame(columns=['ena_accession','mutation','amino_acid_position','new_triplet','original_triplet'])

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

            # insist that the species_min_amount must lie within 0 and 1 (the default is 0.8)
            assert (0 < self.species_min_amount <= 1), "species_min_amount must be > 0 but ≤ 1"

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
            sra_samples=self._interrogate_bigsi(nucleotide_sequence)

            self.number_reference_genomes=len(sra_samples)

            # self.number_reference_genomes=28500
            # instantiate the arrays dictionary that is going to hold all the numpy 2D arrays (this maybe could be done with pandas)
            self.arrays={}

            # make some boring 2D arrays
            self.arrays["original_amino_acid"]=(numpy.broadcast_to(self.amino_acids,(self.number_codons,self.number_amino_acids)))
            self.arrays["amino_acid_position"]=(numpy.broadcast_to(self.amino_acid_position,(self.number_codons,self.number_amino_acids)))
            self.arrays["new_amino_acid"]=(numpy.rot90(numpy.broadcast_to(self.amino_acids_of_codons,(self.number_amino_acids,self.number_codons)),3))
            self.arrays["number_genomes"]=numpy.zeros((self.number_codons,self.number_amino_acids),dtype=int)
            # self.arrays["number_reference_genomes"]=self.number_reference_genomes*numpy.ones((self.number_codons,self.number_amino_acids),dtype=int)

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
            self.arrays["number_nucleotide_changes"]=data[()]["arr_number_nucleotide_changes"]

            # now that all the arrays are loaded, re-create the Pandas dataframe (bit hacky)
            self._create_dataframe()

            # if you thought that was hacky, try this
            stem=filename.split("-summary.npy")[0]

            print(stem)

            self.mutations=pandas.read_csv(stem+"-mutations.csv")


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

    def calculate_metrics(self,threshold=10):

        self.dnds_lookup={}
        for codon in self.codons:
            original_amino_acid=self.triplet_to_amino_acid[codon]
            syn=0
            nonsyn=0
            for position in [1,2,3]:
                for base in ['a','t','c','g']:
                    new_codon=codon[:position-1]+base+codon[position:]
                    if new_codon!=codon:
                        new_amino_acid=self.triplet_to_amino_acid[new_codon]
                        if new_amino_acid==original_amino_acid:
                            syn+=1
                        else:
                            nonsyn+=1
            self.dnds_lookup[codon]=(nonsyn,syn)

        def normalise_number_genomes(row):

            # find out the tuple of (nonsyn,syn) which adds up to nine
            tmp=self.dnds_lookup[row['original_triplet']]

            # if this is a non-synonmous mutation, then normalise
            if row['non_synonymous']:
                return (row['number_genomes']/tmp[0])
            else:
                # if it isn't, and there are no routes to make a synonmous mutation
                if tmp[1]==0:
                    return (row['number_genomes'])
                else:
                    return (row['number_genomes']/tmp[1])

        self.df['normalised_number_genomes']=self.df.apply(normalise_number_genomes,axis=1)

        # create a dataset only containing mutations
        mutations=self.df.loc[self.df.mutation!="-"]

        self.df_mutations_sum=mutations.groupby('amino_acid_position').number_genomes.sum()
        self.df_mutations_num=mutations.groupby('amino_acid_position').number_genomes.count()
        self.df_mutations_max=mutations.groupby('amino_acid_position').number_genomes.max()
        self.df_nucleotide_changes_max=mutations.groupby('amino_acid_position').number_nucleotide_changes.max()
        self.df_nucleotide_changes_num=mutations.groupby('amino_acid_position').number_nucleotide_changes.count()

        synonymous=mutations.loc[mutations.synonymous==True]
        self.df_syn_sum=synonymous.groupby('amino_acid_position').normalised_number_genomes.sum()

        non_synonymous=mutations.loc[mutations.non_synonymous==True]
        self.df_non_syn_sum=non_synonymous.groupby('amino_acid_position').normalised_number_genomes.sum()

        dnds={}
        # only need to consider amino acid positions where there is a non-synonymous mutation, as by definition is there isn't, then dsdn=0
        for i in self.df_non_syn_sum.keys():

            # check if there are also some synonymous mutations at this position
            if i in self.df_syn_sum.keys():

                # now check to see if there are enough non-syn to register a signal!
                if (self.df_non_syn_sum[i]>=threshold):

                    # and enough syn
                    if (self.df_syn_sum[i]>=threshold):

                        dnds[i]=self.df_non_syn_sum[i]/self.df_syn_sum[i]

                    # otherwise there are some syn, but not enough, and enough non-syn
                    else:

                        dnds[i]=150

                # if there aren't enough non-syn, just set to zero as there isn't enough for a signal
                else:
                    dnds[i]=0

            # otherwise there are no synoymous mutations so set to zero as not measurable (inf)
            else:
                if (self.df_non_syn_sum[i]>=threshold):
                    dnds[i]=150.
                else:
                    dnds[i]=0

        else:

            dnds[i]=0.

        self.df_dnds=pandas.Series(data=dnds)

        return(True)

    def __repr__(self):
        """ Overload so that the print() function outputs a summary of the instance.
        """

        amino_acids_string=''.join(self.amino_acids)
        nucleotides_string=''.join(self.nucleotides)

        if self.species_name is None:
            species="All species considered"
        else:
            species=self.species_name

        return("%s\n%s gene\n%s to %s\n%s\n%s\n%i genomes found with this sequence" % (species,self.gene_name,amino_acids_string[0]+str(self.first_amino_acid_position),amino_acids_string[-1]+str(self.first_amino_acid_position+self.number_amino_acids),amino_acids_string,nucleotides_string,self.number_reference_genomes))

    def __len__(self):
        """ Overload so len() simply returns the number of nucleotides stored within the instance of the class.
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
            if triplet_position==1:
                nucleotides_before=nucleotides_string[nucleotide_idx-30:nucleotide_idx]
                nucleotides_after=nucleotides_string[nucleotide_idx+1:nucleotide_idx+33]
            elif triplet_position==2:
                nucleotides_before=nucleotides_string[nucleotide_idx-31:nucleotide_idx]
                nucleotides_after=nucleotides_string[nucleotide_idx+1:nucleotide_idx+32]
            elif triplet_position==3:
                nucleotides_before=nucleotides_string[nucleotide_idx-32:nucleotide_idx]
                nucleotides_after=nucleotides_string[nucleotide_idx+1:nucleotide_idx+31]
            nucleotide_original=nucleotides_string[nucleotide_idx]

            permutations=['a','c','t','g']

            # only consider the original triplet once
            if triplet_position in (1,3):
                 permutations.remove(nucleotide_original)

        else:

            # from the 0-based triplet index, work out the 0-based nuclotide index
            nucleotide_idx=(3*(triplet_idx))

            # split the nucleotide sequence
            nucleotides_before=''.join(self.nucleotides[nucleotide_idx-30:nucleotide_idx])
            tmp=self.nucleotides[nucleotide_idx:nucleotide_idx+3]
            triplet_original=''.join(i for i in tmp)
            nucleotides_after=''.join(self.nucleotides[nucleotide_idx+3:nucleotide_idx+33])

            # codons is a numpy array, so convert to a list so that remove will work
            permutations=self.codons.tolist()

        first_pass=True

        for new_sequence in permutations:

            # create the query sequence
            query_sequence=nucleotides_before+new_sequence+nucleotides_after

            # find out the resulting new triplet
            if triplet_position:
                if triplet_position==1:
                    new_triplet=new_sequence+nucleotides_after[:2]
                elif triplet_position==2:
                    new_triplet=nucleotides_before[-1]+new_sequence+nucleotides_after[0]
                elif triplet_position==3:
                    new_triplet=nucleotides_before[-2:]+new_sequence
            else:
                new_triplet=new_sequence

            # ..and therefore the new amino acid
            new_aminoacid=self.triplet_to_amino_acid[new_triplet]

            # count how many occurences are in the SRA by calling the BIGSI instance
            sra_samples=self._interrogate_bigsi(query_sequence)

            if original_triplet!=new_triplet:
                mut="%s%s%s" % (original_aminoacid,aminoacid_number,new_aminoacid)
                if original_aminoacid==new_aminoacid:
                    synoymous=True
                else:
                    synoymous=False
                for ena_accession in sra_samples:
                    self.mutations=self.mutations.append({'ena_accession':ena_accession, 'mutation':mut,'new_triplet':new_triplet, 'original_triplet':original_triplet, 'amino_acid_position':aminoacid_number, 'original_amino_acid':original_aminoacid, 'new_amino_acid':new_aminoacid,'synoymous':synoymous,'species0_name':sra_samples[ena_accession][0][0],'species0_percentage':sra_samples[ena_accession][0][1],'species1_name':sra_samples[ena_accession][1][0],'species1_percentage':sra_samples[ena_accession][1][1]},ignore_index=True)

            total=len(sra_samples)

            # work out which row this codon belongs to (we already know that the col==triplet_idx)
            row=numpy.where(self.codons==new_triplet)[0][0]

            # if total>0:
            #     print(nucleotides_before,new_sequence,nucleotides_after,aminoacid_number, triplet_idx, nucleotide_idx,triplet_idx+1-10,original_triplet,new_triplet,row,total,len(query_sequence))

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
        sra_samples={}

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

            # pull out the list of predicted species, and their predicted amounts
            predicted_species_list=result[query_string]['results'][sra_sample_number]["species"].split("; ")

            print(predicted_species_list)

            record_sample=False

            if ":" in predicted_species_list[0]:
                predicted_species_name_0 = predicted_species_list[0].split(" : ")[0]
                predicted_species_amount_0 = float(predicted_species_list[0].split(" : ")[1].split("%")[0])

                if len(predicted_species_list)>2:
                    predicted_species_name_1 = predicted_species_list[1].split(" : ")[0]
                    predicted_species_amount_1 = float(predicted_species_list[1].split(" : ")[1].split("%")[0])

                if self.species_name is not None:
                    if predicted_species_name_0==self.species_name:
                        if (predicted_species_amount_0/100.)>=self.species_min_amount:
                            record_sample=True
                    elif len(predicted_species_list)>2:
                        if predicted_species_name_1==self.species_name:
                            if (predicted_species_amount_1/100.)>=self.species_min_amount:
                                record_sample=True
                else:
                    record_sample=True

            if record_sample:
                if sra_sample_number not in sra_samples.keys():
                    if len(predicted_species_list)>2:
                        sra_samples[sra_sample_number]=[(predicted_species_name_0,predicted_species_amount_0),(predicted_species_name_1,predicted_species_amount_1)]
                    else:
                        sra_samples[sra_sample_number]=[(predicted_species_name_0,predicted_species_amount_0),(None,None)]


        return(sra_samples)

    def _create_dataframe(self):
        """Create a Pandas dataframe that lists all the mutations and their occurences, one per line"""

        data_dict={}

        # create a Boolean array of only those positions where sequences have been identified
        positive_elements=self.arrays["number_genomes"]>0

        for key in ['amino_acid_position','original_triplet','new_triplet','number_nucleotide_changes','mutation','number_genomes','original_amino_acid','new_amino_acid','synonymous','non_synonymous']:
            data_dict[key]=(self.arrays[key][positive_elements]).tolist()

        self.df=pandas.DataFrame(data=data_dict)

        self.df["number_nucleotide_changes"]=self.df["number_nucleotide_changes"].astype("int8")
