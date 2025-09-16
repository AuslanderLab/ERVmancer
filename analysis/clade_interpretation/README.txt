This folder contains files which assist in interpreting the output of ERVmancer.

convert_clade_to_elements_in_clade.csv
This file is a two column CSV. The first column (clade_or_herv) is the truncated name of a HERVd HERV or a clade in the phylogenetic tree.
The second column (elements_within_herv_or_clade) contains all of the families or subtypes present in that herv or clade.
For clades, the pattern is element-X, where X is the number of times that element appears in that clade.
For leaf HERVs, the full HERV name is returned.

clades_and_elemnents_of_clades.json
This file is a json dictionary. Each key is a clade, and each value is a list of all of the clades and HERVs that are children or subchildren of that clade.