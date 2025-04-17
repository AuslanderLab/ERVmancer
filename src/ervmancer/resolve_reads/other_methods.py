# This is a python file for the function for the third entry point for ERVmancer
# This is mapping a count CSV file onto the tree.

# imports---------------------------------------------------
import pandas as pd

# Function---------------------------------------------------------------------


def assign_tree_for_other_methods(path_to_other_csv: str, path_for_output_csv: str, clade_under_dict: dict, herv_path_dict: dict):
    """Assigns the other methods to our tree.

    Args:
        path_to_other_df (str): the path to the df of read counts from another method. The format MUST conform to the readme format.
        path_for_output_csv (str): the path where the output csv will be exported.
        clade_under_dict (dict): a dictionary with clade as keys and the clade or hervs under that clade as a list. Will be provided to the end user automatically.
        herv_path_dict (dict): a dictionary with hervs as keys and the path back to the root as a list as values. Will be provided to the end user automatically

    Returns:
        A dataframe of the other method with all clades in addition to the ERV leaves.
    """
    other_df = pd.read_csv(
        path_to_other_csv, index_col=0)  # import the other dataframe to map to the tree.

    # We will use a dictionary of clade as key and a list of values for the final pandas dataframe output.
    out_clade_dict = {}

    for clade in clade_under_dict.keys():  # iterate over all of the clades, identify the hervs under that clade, sum the hervs under that clade, and
        # subset to just the hervs under the clade
        # TODO: Check if final version needs a better ERV check.
        herv_under = ['_'.join(i.split('_')[0:2])
                      for i in clade_under_dict[clade] if i in herv_path_dict]
        # since these methods don't have all of the hervs, I will need an extra step to only look at hervs we consider. Will not affect final calculation.
        herv_under = [i for i in herv_under if i in other_df.columns]

        if len(herv_under) > 1:
            out_clade_dict[clade] = list(other_df[herv_under].sum(axis=1))
        elif len(herv_under) == 1:
            out_clade_dict[clade] = list(other_df[herv_under].values)
        else:
            # print('TEST: herv not found under clade')
            # if no hervs are present under that clade, then just give a list of 0
            out_clade_dict[clade] = [0]*other_df.shape[0]

    # next, give some warnings to the end user if their csv differs from the test csv
    # First, give them a list of hervs not present in their dataframe.
    hashed_csv_columns = {i: 1 for i in list(other_df.columns)}
    hervs_not_present = [i for i in list(
        herv_path_dict.keys()) if i not in hashed_csv_columns]
    if hervs_not_present:
        print('These HERVs in the leaves of the tree were not present in your CSV: \n' +
              str(hervs_not_present) + '\n')
    else:
        print('All HERVs at the leaves of the tree are present in the columns of your CSV')

    # second, give them a list of columns not present in the hervs
    columns_not_hervs = [i for i in list(
        hashed_csv_columns.keys()) if i not in herv_path_dict]
    if columns_not_hervs:
        print('These columns of your CSV were not present in the list of HERVs at the leaves of the tree: \n' +
              str(columns_not_hervs) + '\n')
    else:
        print('No additional non-HERV columns found')

    # convert into the output dataframe

    clade_df = pd.DataFrame(out_clade_dict, index=other_df.index)
    final_df = pd.concat([clade_df, other_df], axis=1)

    # export as a csv
    final_df.to_csv(path_for_output_csv)
