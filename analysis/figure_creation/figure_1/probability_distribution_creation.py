#This is a python file to create the probability distribution dictionaries for ERVmancer.
#Note: Input and output files are too large for github, contact apatterson@wistar.org for input data.

#Imports--------------------------------------

import pandas as pd
import numpy as np
import pickle
from Bio import SeqIO
from tqdm import tqdm
import random

#Functions------------------------------------
def retrieve_pickled_python_obj(pathname):
    '''
    Gets a pickled python object at the path and name
    :param pathname: the path and name to the matrix, make sure it has the .pickle at the end.
    :return: whatever object is at the path and name location
    '''
    with open(pathname, "rb") as handle: pickled_obj = pickle.load(handle)
    return pickled_obj

def determine_lowest_common_clade_final(virus_list, path_dict):
    '''
    Given a list of viruses and a dictionary with viruses as keys and the path back to the "root" as values
    Determines the lowest common clade of those viruses, or returns that virus if there is only one virus
    This version uses the list method, but with an additional step to always use the shortest list
    '''
    if len(virus_list) == 1:
        return virus_list[0]
    else:

        #get the paths of these viruses
        virus_path_list = [path_dict[i] for i in virus_list]

        inter_set = set.intersection(*map(set,virus_path_list)) #this is unordered, so now we need to find the first one

        v_path_for_lca = virus_path_list[np.argmin([len(v) for v in virus_list])] #get the shortest virus path list, to speed up computation in the next step
        lca = [i for i in v_path_for_lca if i in inter_set][0] #this gets the lowest clade that contains all of the given viruses

    return lca

def reverse_complement(input_seq):
    '''
    Converts the str input_seq into its reverse complement
        '''

    input_seq = input_seq.upper()

    temp_seq = input_seq.replace('A','t')
    temp_seq = temp_seq.replace('T','a')
    temp_seq = temp_seq.replace('G','c')
    temp_seq = temp_seq.replace('C','g')

    temp_seq = temp_seq[::-1].upper()

    return temp_seq

def herv_to_kmer_dict_changed(herv, herv_fasta_pathname, kmer_size):
    '''
    Takes in a herv_name, a path to a herv fasta, and a kmer size and outputs a list of all of the possible 75 BP kmers.
    Imports the herv and chops it up into reads of size kmer
    :param herv: the herv name
    :param herv_fasta_pathname: The path to the herv fasta file
    :param kmer_size: the size of the kmer used.
    '''

    herv_fasta = SeqIO.parse(herv_fasta_pathname,'fasta')

    for record in herv_fasta:
        if  '_'.join(record.id.split('_')[0:2]) == herv:
            herv_seqio = record
            break

    f_name, f_sequence = herv_seqio.id, str(herv_seqio.seq)
    seq_len = len(f_sequence)


    #just in case the sequence is smaller than the window, we want to just return the entire sequence to avoid an error
    num_windows = seq_len-kmer_size
    #iterate over the kmer and extract the sequences (need num_windows + 1 due to how the range function works)
    output_list = []

    for win_start in range(0, num_windows+1):
        forward_kmer = f_sequence[win_start:win_start + kmer_size]

        #we need to also consider the reverse complement, so let's make that
        reverse_kmer = reverse_complement(forward_kmer)
        #add both of these into a list
        output_list.append(forward_kmer)
        output_list.append(reverse_kmer)


    return output_list

def herv_clade_pair_to_probability(herv, clade, kmer_to_herv_dict, herv_path_dict, herv_fasta_pathname, kmer_size = 75):
    '''
    converts every (herv, clade) pair into a probability of the reads from that HERV being assigned to that clade.
    :param herv: the name of the herv
    :param clade: the name of the clade
    :param kmer_to_herv_dict: a dictionary of kmer as keys and originator hervs as values.
    :param herv_path_dict: a dictionary with herv as keys and the path back to the unrooted "root" as values.
    :param herv_fasta_pathname: the path to the fasta containing the herv sequences.
    :param kmer_size: the kmer size for the "reads" that will be generated (not the kmer_to_herv_dict kmers)
    :return: The probability a read from the herv is assigned to the specified clade.
    '''

    all_reads_in_herv = herv_to_kmer_dict_changed(herv, herv_fasta_pathname, kmer_size) #convert that HERV into a list of reads

    reads_in_clade_count = 0

    for read in all_reads_in_herv: #iterate over all of the reads

        if determine_lowest_common_clade_final(kmer_to_herv_dict[read], herv_path_dict) == clade: #if the lca of this read read maps to the current clade

            reads_in_clade_count += 1 #add one to the clade count

    return reads_in_clade_count/len(all_reads_in_herv)

def create_herv_clade_pair_dictionary(herv_path_dict, kmer_to_herv_dict, herv_fasta_pathname, kmer_size = 75):
    '''
    Runs herv_clade_pair_to_probability over all hervs in the herv_path_dict dictionary.
    :param herv_path_dict: a dictionary of hervs as keys and the path back to the root as values.
    :param kmer_to_herv_dict: a dictionary of kmers as keys and the kmer originator hervs as values.
    :param herv_fasta_pathname: the path to the fasta file with the herv sequences.
    :param kmer_size: the kmer size for the "reads" that will be generated (not the kmer_to_herv_dict kmers)
    :return: a dictionary with herv, clade pair as key and the probability a read will be assigned from that herv to that clade as a value.
    '''

    output_herv_clade_pair_probability_dict = {} #The output dictionary, with key as a tuple of herv and clade, and the value as the probability a read from that herv will map to that clade.

    for herv in tqdm(list(herv_path_dict.keys())): #iterate over every herv

        for clade in [herv] + herv_path_dict[herv]: #iterate over every clade in the path back to the root, including the herv itself first.

            tuple_key = (herv, clade) #this is the key for the dictionary

            #add the key into the dictionary, with the value set to the probability of mapping a read from herv to clade

            output_herv_clade_pair_probability_dict[tuple_key] = herv_clade_pair_to_probability(herv, clade, kmer_to_herv_dict, herv_path_dict, herv_fasta_pathname, kmer_size)


    return output_herv_clade_pair_probability_dict

def distribution_dict_to_cumulative_dict(pair_to_dist_dict, herv_path_dict):
    '''
    Converts the pair_to_distribution_dict into a pair_to_cumulative_distribution_dictionary
    :param pair_to_dist_dict: A dictionary with tuple of (HERV, Clade) as key and distribution as value.
    :param herv_path_dict: a dictionary of hervs as keys and the path back to the unrooted "root" as values.
    :return dict: a dictionary of tuple (HERV, Clade) as key and cumulative distribution as value
    '''

    pair_to_cumulative_dict = {}

    for herv in herv_path_dict.keys():

        herv_path_back = herv_path_dict[herv] #this gets the path back from the first lca above the herv to the root

        #sometimes, some hervs do not have any unique kmers. If this happens, set the cumulative probability to 0
        if (herv, herv) in pair_to_dist_dict:

            cumulative_prob = pair_to_dist_dict[(herv, herv)] #The first cumulative probability is the herv, herv tuple from the pair_to_dist_dict

        else: #sometimes, some hervs do not have any unique kmers. If this happens, set the cumulative probability to 0
            cumulative_prob = 0

        pair_to_cumulative_dict[(herv, herv)] = cumulative_prob

        for clade in herv_path_back: #now tracing the clade backwards, we update the cumulative probability.

            if (herv, clade) in pair_to_dist_dict: #only update the cumulative probability if these are present, otherwise the cumulative probability does not update, but we still assign to the previous cumulative probability.

                cumulative_prob += pair_to_dist_dict[(herv, clade)] #we update the cumulative probability using the herv, lca tuple

            pair_to_cumulative_dict[(herv, clade)] = cumulative_prob #Add the tuple key and the cumulative probability to the cumulative dictionary.

    #I am adding in one final step, to take care of the computer numbers by rounding the cumulative probabilities to a very long decimal.

    pair_to_cumulative_dict = {i:np.round(j, 9) for i,j in zip(pair_to_cumulative_dict.keys(), pair_to_cumulative_dict.values())}

    return pair_to_cumulative_dict

if __name__ == "__main__":
    kmer_to_herv_dict = retrieve_pickled_python_obj('kmer_75_herv_no_filter.pickle')
    herv_path_dict = retrieve_pickled_python_obj('cleaned_dictionaries/cleaned_herv_path_dict.pkl')
    pair_distribution_dict = create_herv_clade_pair_dictionary(herv_path_dict, kmer_to_herv_dict, 'active_hervs_fasta.fa')
    cumulative_dictionary = distribution_dict_to_cumulative_dict(pair_distribution_dict, herv_path_dict)
    
    with open("saved_probabilities/cumulative_probability_dictionary.pkl", "wb") as f:
        pickle.dump(cumulative_dictionary, f)

    with open("saved_probabilities/not_cumulative_probability_dictionary.pkl", "wb") as f:
        pickle.dump(pair_distribution_dict, f)