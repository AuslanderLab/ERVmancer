import argparse
import logging
import os
import sys
import pickle
import subprocess
import pkg_resources
from tqdm import tqdm
import hashlib
import random
import string
from ervmancer.preprocessing.filter_reads import ReadFilter
from ervmancer.preprocessing.subset_reads import extract_out_original_reads_by_subset
from ervmancer.kmer.process_kmer import parse_fastq_run_kmer_return_dict
from ervmancer.multimap.multimapped_sam import single_multimapped_sam_to_dictionary
from ervmancer.resolve_reads.resolve_and_export import resolve_reads_single_sample_final, convert_final_dict_into_csv
from ervmancer.resolve_reads.other_methods import assign_tree_for_other_methods


def retrieve_pickled_python_obj(pathname: str):
    """Gets a pickled python object at the path and name.

    Args:
        pathname (str): the path and name to the matrix, make sure it has the .pickle at the end.

    Returns:
        pickle object at the designated path
    """
    with open(pathname, "rb") as handle:
        pickled_obj = pickle.load(handle)
    return pickled_obj


def generate_random_hash(length: int = 8):
    """Generates a random string to hash using md5 encoding/hashing

    Args:
        length (int, optional): Length of desired hash. Defaults to 8.

    Returns:
        _type_: _description_
    """
    # Create a random string
    random_str = ''.join(random.choices(
        string.ascii_lowercase + string.digits, k=16))
    # Create a hash from it
    hash_obj = hashlib.md5(random_str.encode())
    # Return the first few characters of the hash
    return hash_obj.hexdigest()[:length]


def get_unique_id():
    """Wrapper function to generate a unique ID hash

    Returns:    
        str: random hash string with size @length
    """
    return generate_random_hash()


def run_command(cmd: str):
    """Runs command using subprocess function

    Args:
        cmd (str): desired command to run

    Raises:
        e: An error if the command exits with any failure exit codes
    """
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed: {cmd}")
        raise e


def get_data_path(filename: str, subdir: str = None):
    """Uses pkg_resources to find data paths in src/ervmancer/data

    Args:
        filename (str): filename of desired data object in data
        subdir (_type_, optional): subdirectory for specific path joining if subdir in data. Defaults to None.

    Returns:
        str: absolute path to data module's desired file for usage in main method
    """
    try:
        # pkg_resources for installed package
        data_path = 'data'
        if subdir:
            data_path = os.path.join(data_path, subdir)
        return pkg_resources.resource_filename('ervmancer', os.path.join(data_path, filename))
    except (ImportError, ModuleNotFoundError):
        # dev mode - relative path
        base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        if subdir:
            return os.path.join(base_path, 'data', subdir, filename)
        return os.path.join(base_path, 'data', filename)


def cleanup_intermediate_files(unique_id: str, output_pathname: str, keep_intermediate: bool):
    """Removes intermediate files that belong to the specific process run/command of ervmancer using the unique hash.

    Args:
        unique_id (str): a length n hash that is used to differentiate between files and runs appended to the basename
        output_pathname (str): absolute path for the output base directory
        keep_intermediate (bool): flag to decide whether or not end user wishes to keep intermediate processing files
    """
    # Cleanup files.
    cleanup_dir = os.path.join(output_pathname, 'intermediate_files')
    if not keep_intermediate:
        logging.info("Cleaning up intermediate files...")
        if os.path.exists(cleanup_dir):
            for filename in os.listdir(cleanup_dir):
                file_path = os.path.join(cleanup_dir, filename)
                try:
                    if os.path.isfile(file_path) and unique_id in file_path and '.gitkeep' not in file_path:
                        # if the files that were used in the generation of the final output exist with their hash, delete the file.
                        # this is to prevent array jobs from erasing parallel job intermediate outputs.
                        os.remove(file_path)
                        logging.debug(
                            f"Removed intermediate file: {file_path}")
                except Exception as e:
                    logging.error(f"Error removing {file_path}: {e}")
            logging.info(f"Cleanup of {cleanup_dir} complete")
        else:
            logging.warning(
                f"Cleanup directory {cleanup_dir} does not exist or is not a directory")
    logging.info("Processing complete!")


def main():
    parser = argparse.ArgumentParser(description='ervmancer')
    parser.add_argument('--kd', '--kmer_dict', dest='kmer_dict', required=False,
                        help="Absolute path to the kmer dictionary file (clean_kmer_31_60_percent_cutoff.pkl).")
    parser.add_argument('--b', '--bt2', dest='b', type=str, required=False,
                        help="User provided absolute path to self provided bowtie2 alignment file.")
    parser.add_argument('--advanced', type=str, required=False,
                        help="User provided absolute path to CSV file with user-provided read counts from other methods.")
    parser.add_argument('--r1', required=False,
                        help='Absolute path to paired-end R1 fastq file. A file is also needed for normalization')
    parser.add_argument('--r2', required=False,
                        help='Absolute path to paired-end R2 fastq file.')
    parser.add_argument('--s1', required=False,
                        help='Absolute path to single strand S1 fastq file.')
    parser.add_argument('--keep_files', action='store_true', default=False,
                        help="Keeps intermediate outputs from Samtools and Bedtools.")
    parser.add_argument('--num_cores', type=int, default=8,
                        help="Number of CPU cores used for processing.")
    parser.add_argument('--output_dir', required=True,
                        help='Absolute path to base output directory for CSV output and intermediate files. (Output folders are autogenerated if they do not exist already!)')
    parser.add_argument('--bowtie_index', required=False,
                        help="Absolute path to the Bowtie2 index to be used for processing. If not provided, uses the built-in GRCh38_noalt_as index.")

    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s', force=True)

    # PROVIDED DICTIONARIES
    if args.kmer_dict is None and args.b and not args.advanced:
        # If using b mode (which needs kmer dict) but no dict is provided and not in advanced mode
        logging.error(
            "Kmer dictionary file path (--kmer_dict) is required when using --b without --advanced.")
        sys.exit(1)
    elif args.kmer_dict is not None and not os.path.exists(args.kmer_dict) and (args.b or not args.advanced):
        # If dict is provided but doesn't exist, and we're in a mode that would use it
        logging.error(f"Kmer dictionary file not found at: {args.kmer_dict}")
        sys.exit(1)
    elif args.kmer_dict is not None:
        # Only set the pathname if a dictionary was provided
        kmer_herv_dict_pathname = os.path.abspath(args.kmer_dict)
    else:
        # In advanced mode without kmer_dict, set to None or a default value
        kmer_herv_dict_pathname = None

    clades_under_dict_pathname = get_data_path('cleaned_clades_under_dict.pkl')
    herv_path_dict_pathname = get_data_path('herv_path_dict.pkl')

    logging.info(f"Using kmer dictionary: {kmer_herv_dict_pathname}")
    try:
        if kmer_herv_dict_pathname is not None:
            kmer_herv_dict = retrieve_pickled_python_obj(
                kmer_herv_dict_pathname)
        else:
            # advanced mode does not need kmer dict
            kmer_herv_dict = None

        clades_under_dict = retrieve_pickled_python_obj(
            clades_under_dict_pathname)
        herv_path_dict = retrieve_pickled_python_obj(herv_path_dict_pathname)
    except Exception as e:
        logging.error(f"Error loading dictionary files: {str(e)}")
        sys.exit(1)

    read_filter = ReadFilter(args.output_dir, args.r1,
                             args.r2, args.s1, args.advanced)
    base_name, paired = read_filter.validate_inputs()
    logging.info(f'Base Name: {base_name}')
    # Create output directories if they do not exist
    for subdir in ['intermediate_files', 'final', 'logs']:
        os.makedirs(os.path.join(args.output_dir, subdir), exist_ok=True)
    unique_id = get_unique_id()

    final_csv_out = read_filter.get_path(
        'final', f'{base_name}_{unique_id}_unified_run_final_out', 'csv')

    if args.advanced and not ((args.r1 and args.r2) or args.s1):
        # ENTRYPOINT THREE - USER PROVIDED READS (FROM OTHER METHODS)
        if os.path.isabs(args.advanced) and os.path.exists(args.advanced) and os.path.isfile(args.advanced):
            # if the given absolute path is valid/exists
            logging.info(f"Assigning tree for given CSV file: {args.advanced}")
            assign_tree_for_other_methods(
                args.advanced, final_csv_out, clades_under_dict, herv_path_dict)
            logging.info(
                f"Assigning reads to tree for {args.advanced} and exporting output complete! Terminating program.")
            sys.exit(1)
        else:
            logging.error(
                "Invalid Path - please provide absolute path to CSV file.")

    # ENTRYPOINT ONE - Check if bowtie2 index is provided and we are using the bowtie2 entrypoint
    if not args.b and not args.advanced:
        if args.bowtie_index is None:
            logging.error(
                "Bowtie2 index path (--bowtie_index) is required when not given a provided bowtie2 sam file or advanced mode.")
            sys.exit(1)
        # Check for the existence of all six bowtie2 index files if were using the bowtie2 process entrypoint
        bt2_extensions = [".1.bt2", ".2.bt2", ".3.bt2",
                          ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]
        bt2l_extensions = [".1.bt2l", ".2.bt2l", ".3.bt2l",
                           ".4.bt2l", ".rev.1.bt2l", ".rev.2.bt2l"]

        # Check if regular or large index files exist
        regular_index_exists = all(os.path.exists(
            f"{args.bowtie_index}{ext}") for ext in bt2_extensions)
        large_index_exists = all(os.path.exists(
            f"{args.bowtie_index}{ext}") for ext in bt2l_extensions)

        if not (regular_index_exists or large_index_exists):
            logging.error(
                "Provided Bowtie index is incomplete. All six index files (.1.bt2, .2.bt2, .3.bt2, .4.bt2, .rev.1.bt2, .rev.2.bt2) or their large index equivalent should exist.")
            sys.exit(1)
        logging.info(f"Using bowtie2 index: {args.bowtie_index}")

    try:
        if args.b:
            # ENTRYPOINT TWO - set paired boolean flag for normalization later and absolute path to sam file.
            try:
                # Convert to absolute path
                outsam_pathname = os.path.abspath(args.b)
                # Check if the path exists
                if os.path.exists(outsam_pathname):
                    # Parse given path as an absolute path for consumption
                    logging.info(
                        f"Using user provided alignment file: {outsam_pathname}")
                    if args.r1 and not args.s1:
                        # TODO: Check with Andrew about this 
                        # if they provided a paired end file then set boolean flag to paired
                        logging.info(
                            "User provided Bowtie2 input and paired read file for normalization")
                        paired = True
                    elif args.s1 and not args.r1:
                        logging.info(
                            "User provided Bowtie2 input and single strand file for normalization")
                        paired = False
                    else:
                        logging.error(
                            "Invalid path to either s1 or r1 arguments. User must provide a single strand or paired read file for normalization.")
                else:
                    logging.error(
                        "Path to provided alignment file is invalid. File does not exist. Terminating program.")
                    sys.exit(1)
            except Exception as e:
                logging.error(f"Error processing path: {str(e)}")
                sys.exit(1)
        else:
            outsam_pathname = read_filter.get_path(
                'intermediate_files', f'{base_name}_{unique_id}_bowtie2_out', 'sam')

        outbam_pathname = read_filter.get_path(
            'intermediate_files', f'{base_name}_{unique_id}_all_hervs', 'bam')
        sorted_outbam_pathname = read_filter.get_path(
            'intermediate_files', f'{base_name}_{unique_id}_sorted', 'bam')
        converted_herv_gtf_to_bed = get_data_path('hervs_genomic_coords.bed')
        logging.info(f"Using HERV BED file: {converted_herv_gtf_to_bed}")
        outbed_pathname = read_filter.get_path(
            'intermediate_files', f'{base_name}_{unique_id}_multimap', 'bed')
        # Bedtools intersect with HERV GTF file (Multimap)
        subset_outsam_pathname = read_filter.get_path(
            'intermediate_files', f'{base_name}_{unique_id}_hervs_subset', 'sam')
        # Filter unique read appearances (step for KMER)
        subset_outbam_pathname = read_filter.get_path(
            'intermediate_files', f'{base_name}_{unique_id}_all_hervs.bwa_read', 'bam')

        commands = [
            f"samtools view -bS {outsam_pathname} -o {outbam_pathname}",
            f"samtools sort -@ {args.num_cores} {outbam_pathname} -o {sorted_outbam_pathname}",
            f"samtools index {sorted_outbam_pathname}",
            # Multimap uses this intermediate output below in the final dir
            f"bedtools intersect -abam {sorted_outbam_pathname} -b {converted_herv_gtf_to_bed} -wa -wb -bed > {outbed_pathname}",
            f"bedtools intersect -abam {sorted_outbam_pathname} -b {converted_herv_gtf_to_bed} > {subset_outbam_pathname}",
            # Kmer step uses this intermediate output in the final dir
            f"samtools view -h -o {subset_outsam_pathname} {subset_outbam_pathname}"
        ]
        if not args.b:
            # ENTRYPOINT ONE
            # if bowtie file is not provided by user, add the bowtie 2 commands with unified mode parameters
            if paired:
                bt_cmd = (f"bowtie2 -p {args.num_cores} --very-sensitive --end-to-end -X 1500 "
                          f"--no-mixed --no-discordant --no-dovetail --no-unal --score-min L,-0.1,-0.1 "
                          f"-x {args.bowtie_index} -1 {read_filter.r1_path} -2 {read_filter.r2_path} "
                          f"-k 100 -S {outsam_pathname} "
                          f"2> {read_filter.get_path('logs', base_name, f'bt2_paired_{unique_id}.err')}")
            else:
                bt_cmd = (f"bowtie2 -p {args.num_cores} -N 1 -L 10 --very-sensitive --end-to-end --no-unal "
                          f"-x {args.bowtie_index} --score-min L,-0.1,-0.1 -U {read_filter.s1_path} "
                          f"-k 100 -S {outsam_pathname} "
                          f"2> {read_filter.get_path('logs', base_name, f'bt2_single_{unique_id}.err')}")
            # insert the bowtie 2 command to the front of the pipeline commands queue
            commands.insert(0, bt_cmd)

        logging.info("Starting processing pipeline...")
        logging.info(f"Running commands for file with hash: {unique_id}")
        for cmd in tqdm(commands, desc="Processing commands"):
            logging.info(f"\nExecuting: {cmd}")
            run_command(cmd)

        pathname_extracted_reads = read_filter.get_path(
            'intermediate_files', f'{base_name}_{unique_id}_extracted_reads', 'txt')
        # filter each read so that it appears once for KMER method
        logging.info("Extracting original reads for kmer analysis...")
        if paired:
            extract_out_original_reads_by_subset(
                read_filter.r1_path,
                subset_outsam_pathname,
                pathname_extracted_reads,
                read_filter.r2_path
            )
        else:
            extract_out_original_reads_by_subset(
                read_filter.s1_path,
                subset_outsam_pathname,
                pathname_extracted_reads
            )
        # filter each read so that it appears once for KMER method
        logging.info("Executing kmer method - 31 bp")
        if paired:
            # KMER STEP: Returns a dictionary of reads as keys and kmer assignment as values
            read_to_kmer_assignment_dict = parse_fastq_run_kmer_return_dict(input_sam_pathname=pathname_extracted_reads,
                                                                            kmer_herv_dict=kmer_herv_dict,
                                                                            path_dict=herv_path_dict,
                                                                            kmer_size=31,
                                                                            paired_end=True)
        else:
            read_to_kmer_assignment_dict = parse_fastq_run_kmer_return_dict(input_sam_pathname=pathname_extracted_reads,
                                                                            kmer_herv_dict=kmer_herv_dict,
                                                                            path_dict=herv_path_dict,
                                                                            kmer_size=31,
                                                                            paired_end=False)

        logging.info("Executing multimap method")

        read_to_multimap_assignment_dict = single_multimapped_sam_to_dictionary(input_bed_pathname=outbed_pathname,
                                                                                path_dict=herv_path_dict)
        logging.info(
            "Resolving reads and converting final dictionary into CSV")
        # resolve the differences between the kmer and the multimapping approaches.
        final_resolved_dict = resolve_reads_single_sample_final(kmer_dict=read_to_kmer_assignment_dict,
                                                                multi_dict=read_to_multimap_assignment_dict,
                                                                clade_dict=clades_under_dict)
        if paired:
            convert_final_dict_into_csv(original_fastq_pathname_for_normalization=args.r1,
                                        resolved_dict=final_resolved_dict,
                                        clade_dict=clades_under_dict,
                                        herv_path_dict=herv_path_dict,
                                        output_path=final_csv_out)
        else:
            convert_final_dict_into_csv(original_fastq_pathname_for_normalization=args.s1,
                                        resolved_dict=final_resolved_dict,
                                        clade_dict=clades_under_dict,
                                        herv_path_dict=herv_path_dict,
                                        output_path=final_csv_out)

        cleanup_intermediate_files(unique_id, args.output_dir, args.keep_files)

    except Exception as e:
        logging.error(f"An error occurred: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
