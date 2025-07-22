# ERVmancer Advanced Function Documentation

# Explanation of Advanced Function

Entrypoint 3 of ERVmancer is recommended only for advanced users. The use case is mapping HERV counts generated from another method onto ERVmancer's phylogenetic tree. ERVmancer uses a phylogenetic tree of 3043 HERVs with protein coding potential. In the normal process of ERVmancer, reads are mapped to the lowest common ancestor of the phylonegenic tree, resolving ambiguity in the best possible way given the read information. However, methods which do not account for this ambiguity can also be mapped to the phylogenetic tree. This involves taking in a CSV file that contains columns found in the leaves (HERV counts) of the phylogenetic tree, and then adding the internal nodes of the phylogenetic tree into the CSV. In this case, the counts of the internal nodes are the counts within all of the leaves which are found underneath the internal node in the phylogenetic tree.

# CSV Format
The CSV MUST have columns as genes or HERVs, and rows as samples. The columns can contain a mix of genes and HERVs. The HERVs must match the exact format provided in herv_documentation.csv column "herv_name_truncated". The HERV names follow the HERVdb numeric identity format. Any column that does not exactly match the labels of the 3043 HERVs in the leaves of the phylogenetic tree will be ignored. A list of the columns ignored is provided in the output log. Additionally, any of the 3043 HERVs that are NOT found in the CSV file will be treated as zero in the phylogenetic tree, and a list of these missing HERVs will be provided in the output log file as well.

An example format can be found in the file "counts_ART_10_for_example_added_and_removed_colums.csv". This file contains additional genes and removed HERV columns for illustration purposes.

# Commands

ervmancer --advanced ABSOLUTE_PATH_TO_DATA --output_dir ABSOLUTE_PATH_TO_OUTPUT_DIRECTORY

If the folders "final", "intermediate_files", and "logs" are not present in the output directory, they will be created. The output csv will be in the "final" folder with the same name as the original CSV. The "logs" directory will contain a hashed log which will detail the columns mapped, columns ignored, and HERVs missing.

# Output
The final result will be a copy of the provided CSV with 3042 additional columns. These additional columns are the internal nodes found in the ERVmancer phylogenetic tree. Values of the internal nodes are the summed counts of all of the HERVs that are found beneath that internal node. 

An example result can be found in the file "mapped_ART_10_for_example_added_and_removed_columns.csv"