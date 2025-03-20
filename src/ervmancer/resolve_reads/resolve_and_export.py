import pandas as pd
import regex as re
import subprocess


def append_dict_to_df(new_dict, df, name_of_sample):
    '''
    Takes a dictionary and a dataframe, and appends the dictionary values as row entries using the keys as columns
    df should already have the column names
    name of sample will be appended to the last entry, for use as index at the end
    '''
    key_list = list(new_dict.keys())
    val_list = list(new_dict.values())

    val_list.append(name_of_sample)
    key_list.append('sample')

    inter_df = pd.DataFrame(val_list).T

    inter_df.columns = key_list
    # out_df = df.merge(inter_df, how = 'outer', on = None) #so this will merge the dataframes together, but will merge on the intersection of the columns.

    return pd.concat([df, inter_df])


def counts_for_everything_in_clade(lca_dict_out, clade_dict, herv_path_dict):
    '''
    For every clade and virus in the phylogenetic tree, calculates the counts in the sample in/under that clade
    :param lca_dict_out: the output of herv_assignment_file_to_clade, which is a dictionary of an LCA or virus as keys and the number of times a read was assigned to that lca as values
    :param clade_dict: the output of create_hervs_and_clades_under_clade_dict, a dictionary of all lca (clades) and everything under those LCA
    :param herv_path_dict: a dictionary of the hervs as keys and the paths to the "root" as a list. Only used to get all of the herv names
    :return: a dictionary of clade as keys and counts of all of the viruses and clades under that clade summed as a value
    :WARNING: FASTQ ID must be uppercase, or else "b" will be removed. This warning I will keep for now, but appears to be before I split out the normalization function separately.
    '''

    # instatiate the dictionary with all hervs and clades, setting counter to 0

    output_dict = {str(i): 0 for i in list(
        clade_dict.keys()) + list(herv_path_dict.keys())}
    total_iterable = list(clade_dict.keys()) + list(herv_path_dict.keys())
    # iterate over all clades in the clade_dict
    for lca in total_iterable:
        # lca = str(lca)
        if lca in herv_path_dict:  # this runs if the lca is just the herv, in that case then we just need to have the direct herv count
            if lca in lca_dict_out:  # but need to check if we actually have this lca in the lca dict before assigment, otherwise it will just remain at zero
                output_dict[lca] = lca_dict_out[lca]
        # extract out all of the subclades and leafs, sum up the counts in the lca_dict_out, only if the clade is in the lca_dict_out (or else would return an error)

        else:  # need to run this if the lca is a clade
            total_count = 0
            # extract out the clades and viruses within this clade
            for subclade in clade_dict[lca]:
                if subclade in lca_dict_out:
                    total_count += lca_dict_out[subclade]
            # have to also add the count for the actual lca
            if lca in lca_dict_out:  # but only if it is actuall in the lca_dict dictionary
                total_count += lca_dict_out[lca]
            output_dict[lca] = total_count

            # output_dict[lca] = np.sum([lca_dict_out[i] for i in clade_dict[lca] if i in lca_dict_out])

    return output_dict


def get_original_fastq_norm_factor(original_fastq_for_normalization):
    '''
    Counts the number of lines in the original_fastq_for_normalization
    :param original_fastq_for_normalization: the combined path and name for the fastq for normalization
    '''
    # Also need normalized. We need to do this on the original file (pre alignment filtered) to get the truly normalized counts.
    normalized_gz_check = r'\.gz+$'  # This will match a .gz only at an ending
    # checks if the .gz is the ending of this pathname
    normalized_check_out = re.search(
        normalized_gz_check, original_fastq_for_normalization)
    if normalized_check_out is not None:  # if we have a gz file, then we need to treat it differently
        # This is from chatgpt
        gunzip_process = subprocess.Popen(
            ['gunzip', '-c', original_fastq_for_normalization], stdout=subprocess.PIPE)
        wc_process = subprocess.Popen(
            ['wc', '-l'], stdin=gunzip_process.stdout, stdout=subprocess.PIPE)
        gunzip_process.stdout.close()
        gunzip_process.stdout.close()
        norm_factor = int(wc_process.communicate()[0])/4
    else:  # If it is not gz, then we can just read the file using wc
        # more chatgpt help
        sub_out = subprocess.run(
            ['wc', '-l', original_fastq_for_normalization], stdout=subprocess.PIPE, text=True)
        # norm_factor = int(sub_out.decode().strip('\n').split(' ')[0])/4
        norm_factor = int(sub_out.stdout.split()[0])/4
    return norm_factor


def convert_final_dict_into_csv(original_fastq_pathname_for_normalization, resolved_dict, clade_dict, herv_path_dict, out_csv):
    '''
    Converts the output of resolve_reads_single_sample_final into a single csv file.
    '''

    # First, need to convert the dictionary into an lca dictionary, IE, clade or herv as key and counts as values.
    count_only_dict = {}

    for val in resolved_dict.values():
        if val in count_only_dict:
            count_only_dict[val] = count_only_dict[val] + 1
        else:
            count_only_dict[val] = 1

    # next, need to create the summed dictionaries as well, and normalize them.

    summed_only_dict = counts_for_everything_in_clade(
        count_only_dict, clade_dict, herv_path_dict)

    # finally, normalize these.
    norm_factor = get_original_fastq_norm_factor(
        original_fastq_for_normalization=original_fastq_pathname_for_normalization)

    # make it counts per million
    summed_normalized_dict = {i: (j/norm_factor)*1000000 for i,
                              j in zip(summed_only_dict.keys(), summed_only_dict.values())}
    count_normalized_dict = {i: (j/norm_factor)*1000000 for i,
                             j in zip(count_only_dict.keys(), count_only_dict.values())}

    # convert these dictionaries into a single csv.

    # instantiate the df
    final_df = pd.DataFrame(columns=list(
        clade_dict.keys()) + list(herv_path_dict.keys()))

    # append the dictionaries to the final df
    final_df = append_dict_to_df(
        summed_normalized_dict, final_df, 'summed_normalized_counts')
    final_df = append_dict_to_df(
        count_normalized_dict, final_df, 'direct_normalized_counts')

    final_df = append_dict_to_df(
        count_only_dict, final_df, 'direct_read_counts')
    final_df = append_dict_to_df(
        summed_only_dict, final_df, 'summed_read_counts')

    final_df = final_df.set_index('sample')
    final_df.to_csv(out_csv)


def resolve_reads_single_sample_final(kmer_dict, multi_dict, clade_dict):
    '''
    takes in all four dictionaries, and resolves the assignments according to the resolution_strategy
    :param resolution_strategy: One of ['keep_kmer_only','keep_multimap_only','keep_highest','keep_lowest']
    :param default_resolution: when the assignments are not a subset of each other and one of ['keep_highest','keep_lowest'], will assign reads to either ['default_multimap','default_kmer']
    '''

    # This setup will always use the keep lowest resolution strategy
    resolved_id_dict = {}  # these will be the keys that we will keep
    # these are used to extract out the reads (again) into a separate file that as the final resolved assignments
    resolved_using_kmer_dict = {}
    # Same thing, this is a dictionary of the multimapped reads resolved
    resolved_using_multimap_dict = {}
    # finally, a dictionary of the reads that were resolved using both at the same level.
    resolved_both_dict = {}

    # for every read id, we need to check the kmer clade assignment and the multimap clade assignment.
    for key in kmer_dict.keys():  # for every read id, we need to check the kmer clade assignment and the multimap clade assignment. This step also filters to only the kmer assigned reads.
        # If they are the same, then add the key and value to resolved_id_dict
        if kmer_dict[key] == multi_dict[key]:
            # since they are the same, just use the kmer assigned dictionary.
            resolved_id_dict[key] = kmer_dict[key]
            # resolved_using_kmer_dict[key] = kmer_dict[key]
            resolved_both_dict[key] = kmer_dict[key]

        # next, consider the cases where kmer is a leaf and multimap is a different leaf. Resolve this case by following the multimapped approach, since it has more information.
        elif ('ERV_' == kmer_dict[key][0:4]) & ('ERV_' == multi_dict[key][0:4]):
            resolved_id_dict[key] = multi_dict[key]
            resolved_using_multimap_dict[key] = multi_dict[key]

        # next, consider the case where the kmer is a leaf and multi is a node with children.
        elif ('ERV_' == kmer_dict[key][0:4]) & ('ERV_' != multi_dict[key][0:4]):
            # if that kmer leaf is in the multimapped subtree, then assign to the herv!
            if kmer_dict[key] in clade_dict[multi_dict[key]]:
                resolved_id_dict[key] = kmer_dict[key]
                resolved_using_kmer_dict[key] = kmer_dict[key]

            # if not, resolve according to the default resolution, which is set to be the multimapped approach. This case will occure when the kmer approach assigned it to something outside of the multimap clade
            else:
                resolved_id_dict[key] = multi_dict[key]
                resolved_using_multimap_dict[key] = multi_dict[key]

        # next, consider the case where the kmer is a node with children and multi is a leaf
        elif ('ERV_' != kmer_dict[key][0:4]) & ('ERV_' == multi_dict[key][0:4]):
            # If the herv is in the clade, then assign to the herv!
            if multi_dict[key] in clade_dict[kmer_dict[key]]:
                resolved_id_dict[key] = multi_dict[key]
                resolved_using_multimap_dict[key] = multi_dict[key]

            # Otherwise, resolve according to the default resolution strategy, which is to default to the multimap.
            else:
                resolved_id_dict[key] = multi_dict[key]
                resolved_using_multimap_dict[key] = multi_dict[key]

        # next, consider the cases where both are clades
        # The first to check is if the kmer clade is under the multi clade, which will resolve to the kmer dictionary clade
        elif kmer_dict[key] in clade_dict[multi_dict[key]]:
            resolved_id_dict[key] = kmer_dict[key]
            resolved_using_kmer_dict[key] = kmer_dict[key]

        # next, see if the multi clade is under the kmer clade. This will resolve to the multimap clade
        elif multi_dict[key] in clade_dict[kmer_dict[key]]:
            resolved_id_dict[key] = multi_dict[key]
            resolved_using_multimap_dict[key] = multi_dict[key]

        # finally, if none of the above worked, just default to whatever the default setting is
        # TODO: Check if this is necessary, maybe with a print statement.
        else:
            resolved_id_dict[key] = multi_dict[key]
            resolved_using_multimap_dict[key] = multi_dict[key]

    # , resolved_using_kmer_dict, resolved_using_multimap_dict, resolved_both_dict
    return resolved_id_dict
