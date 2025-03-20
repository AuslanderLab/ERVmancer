from ervmancer.preprocessing.lca import determine_lowest_common_clade


def determine_lowest_common_clade_final(virus_list, path_dict):
    '''
    Given a list of viruses and a dictionary with viruses as keys and the path back to the "root" as values
    Determines the lowest common clade of those viruses, or returns that virus if there is only one virus
    This version uses the list method, but with an additional step to always use the shortest list
    :param virus_list: a list of the viruses a read could have potentially come frome.
    :param path_dict: a dictionary with viruses as keys and the path back to root as values.
    '''
    if len(virus_list) == 1:  # If the virus is just 1 virus, return that virus as the lca
        return virus_list[0]
    else:

        # get the paths of these viruses
        virus_path_list = [path_dict[i] for i in virus_list]

        # this is unordered, so now we need to find the first one
        inter_set = set.intersection(*map(set, virus_path_list))

        # get the shortest virus path list, to speed up computation in the next step
        v_path_for_lca = virus_path_list[np.argmin(
            [len(v) for v in virus_list])]
        # this gets the lowest clade that contains all of the given viruses
        lca = [i for i in v_path_for_lca if i in inter_set][0]

    return lca


def single_multimapped_sam_to_dictionary(input_bed_pathname, path_dict):
    '''
    Takes in a multimapped sam file from bowtie2 and calculates the LCA for each read
    This version has been updated to take in the bed file of HERV mapped reads, and can only make the dictionaries.
    It doesn't matter if it is paired end or single end in this step.
    :param input_bed_pathname: the intersected bed file.
    :param path_dict: a dictionary of hervs as keys and path back to the 'root' as values, for use in the function determine_lowest_common_clade_final.
    '''

    in_sam = open(input_bed_pathname, 'r')

    id_to_herv_list_dict = {}  # instantiate the read id to herv assignment dictionary
    # out_lca_dict = {} #instantiate the output lca dictionary
    # print(input_sam_pathname)
    out_id_to_assignment_dict = {}

    for line in in_sam:
        if line[0] == '@':  # skip the first lines that don't have reads
            pass
        else:  # we need to get the barcode and save to a dictionary with the virus assignments
            # id_barcode = split_line[0] #First one is the id

            # This gets the id_barcode. There are read-1 and read-2, which are signified by a /2 or a /1 at the end, typically.
            # TODO: Check if /2 is a universal divider of read 1 vs read 2.
            id_barcode = line.split('\t')[3].split('/')[0]

            # virus_assignment = split_line[2] #third is the virus assignment from bowtie2

            virus_assignment = '_'.join(line.split(
                '\t')[-1].split(';')[0].split(' ')[1].strip('"').split('_')[0:2])

            # assign these to the dictionary
            # previous versions would export the sequence as well, however I do not have the sequence information now due to using a bed file.
            if id_barcode in id_to_herv_list_dict:  # if we have already encountered this read before
                # add the virus assignment to the previous assignments
                id_to_herv_list_dict[id_barcode] = id_to_herv_list_dict[id_barcode] + [
                    virus_assignment]

            else:  # if we haven't encountered the id before, then make the id the key and the virus assignment in a list as the value
                id_to_herv_list_dict[id_barcode] = [virus_assignment]

    in_sam.close()

    # once we have the dictionary, we need to first take a set of the lists, just in case there are duplicate assignments (shouldn't be, but just in case)
    # there will be duplicate assignments for paired end, so this step reduces those duplicate assignments
    id_to_herv_list_dict = {i: list(set(j)) for i, j in zip(
        id_to_herv_list_dict.keys(), id_to_herv_list_dict.values())}

    # if we have the optional output file, need to open it
    # if optional_seq_assignment_filename: #Option to save the seq assignments to a separate file
    #     out_file = open(optional_seq_assignment_filename, 'w')

    # now, we need to iterate over every single read-herv key value pair in the dict and find out the LCA assignment
    for read in id_to_herv_list_dict.keys():
        lca_assignment = determine_lowest_common_clade_final(
            virus_list=id_to_herv_list_dict[read], path_dict=path_dict)

        # add these to a dictionary of ID to lca assignment.
        out_id_to_assignment_dict[read] = lca_assignment

    return out_id_to_assignment_dict
