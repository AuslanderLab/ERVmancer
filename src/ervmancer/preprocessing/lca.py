import numpy as np


def determine_lowest_common_clade(virus_list, path_dict):
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
