from ervmancer.preprocessing.lca import determine_lowest_common_clade
from typing import Tuple, Optional


def identify_sequence_herv_barcode(sequence, kmer_herv_dict: dict, kmer_size: int = 31) -> Tuple[list, list]:
    """ Takes in a sequence, breaks it into kmers, for each kmer it extracts out the viruses that kmer could belong to.
        then exports the intersection of those kmer virus lists.

    Args:
        sequence (str): The sequence to assign to the viruses.
        kmer_herv_dict (dict): a dictionary of kmers as keys and HERVs that kmer could belong to as a list
        kmer_size (int, optional): _description_. Defaults to 31.

    Returns:
        Tuple[list, list]: kmer list originated from viruses in the potential candidate pool
    """
    # Cut up the sequence into a list of kmers.
    seq_kmers = chop_sequence_into_kmers(str(sequence), kmer_size=kmer_size)

    all_potential_viruses = []
    kmers_used = []

    for sk in seq_kmers:  # for each sequence, we check if it is in the kmer dictionary. If it is, we append that list to the growing list of potential HERVs.

        if sk in kmer_herv_dict:
            all_potential_viruses.append(kmer_herv_dict[sk])
            kmers_used.append(sk)

    # if its empty, we just need to return an empty list
    if not all_potential_viruses:
        return all_potential_viruses, kmers_used

    # This only runs if we actually have anything in the all_potential_viruses
    # We do the intersection of all of the smaller lists first.
    # This is because, if one or more kmers has more specific information than others, we defer to the more specific kmer.
    intersected = list(set.intersection(*map(set, all_potential_viruses)))
    if intersected:  # if we have a non-empty list, then we return the intersection
        return intersected, kmers_used
    else:  # If the list is empty after intersection, BUT we already know there are hervs present, we return the union instead.
        # This will occur if the kmers could have originated from different viruses, but not the same virus. We will still map to the tree, it will just be higher up to indicate the uncertainty.
        return list(set().union(*all_potential_viruses)), kmers_used


def chop_sequence_into_kmers(seq: str, kmer_size: int = 31) -> list:
    """Cuts up the sequence into kmers of the specified length

    Args:
        seq (str): a string of a dna sequence.
        kmer_size (int, optional): the size of the final kmers. Defaults to 31.

    Returns:
        list: a list of overlapping kmers of size kmer_size covering the original sequence.
    """
    seq_len = len(seq)

    output_list = []

    num_windows = seq_len-kmer_size
    # iterate over the kmer and extract the sequences (need num_windows + 1 due to how the range function works)
    for win_start in range(0, num_windows+1):
        output_list.append(seq[win_start:win_start + kmer_size])

    return output_list


def parse_fastq_run_kmer_return_dict(input_sam_pathname: str, kmer_herv_dict: dict, path_dict: dict, kmer_size: int = 31, paired_end: bool = False) -> dict:
    """Goes over every line in a fasta file, for each sequence identifies the potential hervs it could come from.
       Then, exports the list of assignments to a separate file.

    Args:
        input_sam_pathname (str): The path and name of the fastq file to analyze
        kmer_herv_dict (dict): A dictionary of kmers as keys and the hervs which have that kmer as values
        path_dict (dict): A dictionary of the herv as keys and the paths back to the "root" as values
        kmer_size (int, optional): The size of the kmer in the kmer_herv_dict. Defaults to 31.
        paired_end (bool, optional): Flag to determine whether to process it as a paired read or single strand file. Defaults to False.

    Returns:
        out_id_to_assignment_dict: a dictionary with clade as keys and values as the number of reads DIRECTLY assigned to that LCA (use everything_under_clade_dict if you want to see everything assigned to or UNDER/WITHIN an LCA clade)
    """
    # open (text output from previous method)
    in_file = open(input_sam_pathname, 'r')

    # I need another dictionary that will be the read id and the herv or clade assignment of that id.
    out_id_to_assignment_dict = {}

    if paired_end:  # if the data is paired end, need to run this way. This will calculate the kmers for both paired reads and add them together to determine the virus.
        first_line = in_file.readline()
        # need to do to while loops, the first to chew through the @ lines, the second to go through the actual sequences.
        while '@' == first_line[0]:
            # checks for sam file but TODO: removal because its a text file now with the format we want
            first_line = in_file.readline()
        while first_line:
            # the first and second lines are paired end reads.
            second_line = in_file.readline()
            # since this is paired end, we are going to calculate the kmers for both the first and second line
            # the text file has id first then next row is read secondary thats why were reading two lines
            out_herv_assignment = []  # will add the assignment to this list
            kmers_used = []  # will add the kmers used to this list
            out_herv_assignment_list = []
            kmers_used_list = []

            # iterate over the first and second lines, doing the kmer step for each line. Sum up the kmers and viruses for each and then determine the LCA for the pair.
            for line in [first_line, second_line]:
                split_line = line.split('\t')
                # First one is the id. #TODO: check if we need specifically R1 or R2, at this point it does not look like we need to change it.
                id_barcode = split_line[0]
                # tab delineated text file where the first is the ID
                # no longer a sam file, now the read position is the second (1 pythonic)
                seq = split_line[1]
                out_herv_assignment_subset, kmers_used_subset = identify_sequence_herv_barcode(
                    seq, kmer_herv_dict, kmer_size=kmer_size)
                kmers_used_list.append(kmers_used_subset)
                # since we loop, we add the assignment for both read 1 and read 2 to this, for later lca assignment.
                out_herv_assignment_list.append(out_herv_assignment_subset)
                out_herv_assignment = out_herv_assignment + out_herv_assignment_subset
                kmers_used = kmers_used + kmers_used_subset

            # remove duplicates
            out_herv_assignment = list(set(out_herv_assignment))
            kmers_used = list(set(kmers_used))

            if out_herv_assignment:  # if we actually got an assignment, which means any kmer in either read was a viral kmer

                lca_assignment = determine_lowest_common_clade(
                    virus_list=out_herv_assignment, path_dict=path_dict)  # then assign the virus to an lca
                # add the assignment to the output dictionary.

                out_id_to_assignment_dict[id_barcode] = lca_assignment

            # continue the while loop
            first_line = in_file.readline()

    else:  # if the data is not paired end, then we need to run this way.

        for line in in_file:
            if '@' in line:
                pass
            else:
                split_line = line.split('\t')
                id_barcode = split_line[0]  # First one is the id
                # no longer a sam file, now the read position is the second (1 pythonic)
                seq = split_line[1]
                out_herv_assignment, kmers_used = identify_sequence_herv_barcode(
                    seq, kmer_herv_dict, kmer_size=kmer_size)

                if out_herv_assignment:  # if we actually got an assignment, which means any kmer was a viral kmer

                    lca_assignment = determine_lowest_common_clade(
                        virus_list=out_herv_assignment, path_dict=path_dict)  # then assign the virus to an lca
                    out_id_to_assignment_dict[id_barcode] = lca_assignment

    in_file.close()

    return out_id_to_assignment_dict
