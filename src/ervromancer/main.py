import argparse
import logging
import os
import subprocess
import sys
from tqdm import tqdm
from .preprocessing.filter_reads import ReadFilter


def run_command(cmd):
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed: {cmd}")
        raise e


def main():
    parser = argparse.ArgumentParser(description='ervmancer')
    parser.add_argument('--mode', choices=['strict', 'relaxed'], required=True,
                        help="Alignment bounds (stringent, relaxed)")
    parser.add_argument('--r1', required=False,
                        help='Absolute path to paired-end R1 fastq file.')
    parser.add_argument('--r2', required=False,
                        help='Absolute path to paired-end R2 fastq file.')
    parser.add_argument('--s1', required=False,
                        help='Absolute path to single strand S1 fastq file.')
    parser.add_argument('--keep_files', action='store_true', default=False,
                        help="Keeps intermediate outputs from Samtools and Bedtools.")
    parser.add_argument('--num_cores', type=int, default=8,
                        help="Number of CPU cores used for processing.")
    parser.add_argument('--output_dir', required=True,
                        help='Absolute path to output directory for CSV output and intermediate files.')
    parser.add_argument('--bowtie_index', required=True,
                        help="Absolute path to the Bowtie2 index to be used for processing.")

    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s', force=True)

    try:
        # Create output directories
        for subdir in ['intermediate_files', 'final', 'logs']:
            os.makedirs(os.path.join(args.output_dir, subdir), exist_ok=True)

        read_filter = ReadFilter(args.output_dir, args.r1, args.r2, args.s1)
        base_name, paired = read_filter.validate_inputs()
        logging.info(f'Base Name: {base_name}')

        os.makedirs(os.path.join(args.output_dir,
                    "intermediate_files"), exist_ok=True)

        sam_path = read_filter.get_path(
            'intermediate_files', f'{base_name}_bowtie2_all_hervs.bwa.read', 'sam')
        bam_path = read_filter.get_path(
            'intermediate_files', f'{base_name}_all_hervs.bwa.read', 'bam')
        indexed_bam = read_filter.get_path(
            'intermediate_files', f'{base_name}_all_hervs.bwa.read.indexed', 'bam')
        sorted_bam = read_filter.get_path(
            'intermediate_files', f'{base_name}_all_hervs.bwa.read.sorted', 'bam')
        bed_out = read_filter.get_path(
            'intermediate_files', f'{base_name}_all_hervs.bwa.read', 'bed')
        # Bedtools intersect with HERV GTF file (Multimap)
        intersect_file = read_filter.get_path(
            'final', f'{base_name}_all_hervs_intersected.bwa.read', 'bed')
        # Filter unique read appearances (step for KMER)
        hervs_read = read_filter.get_path(
            'final', f'{base_name}_all_hervs.bwa_read', 'sam')

        logging.info(f"Sam Path: {sam_path}")
        logging.info(f"Sam Path: {bam_path}")

        if args.mode == "strict":
            if paired:
                logging.info("Going into strict mode, paired reads")
                bt_cmd = f"""bowtie2 -p {args.num_cores} --very-sensitive --end-to-end -X 1500 \
                --no-mixed --no-discordant --no-dovetail --no-unal --score-min L,0,0.1 \
                -x {args.bowtie_index} -1 {read_filter.r1_path} -2 {read_filter.r2_path} -k 100 -S {sam_path} \
                2> {read_filter.get_path('logs', base_name, 'bt2_paired_strict.err')}"""
            else:
                logging.info("Going into strict mode, single strand")
                bt_cmd = f"""bowtie2 -p {args.num_cores} -N 1 -L 10 --very-sensitive --end-to-end --no-unal \
                -x {args.bowtie_index} --score-min L,0,0.1 -U {read_filter.s1_path} -k 100 -S {sam_path} \
                2> {read_filter.get_path('logs', base_name, 'bt2_single_strict.err')}"""
        else:  # relaxed mode
            if paired:
                logging.info("Going into relaxed mode, paired reads")
                bt_cmd = f"""bowtie2 -p {args.num_cores} --very-sensitive --end-to-end --no-unal --no-mixed --no-discordant \
                -x {args.bowtie_index} -1 {read_filter.r1_path} -2 {read_filter.r2_path} -k 100 -S {sam_path} \
                2> {read_filter.get_path('logs', base_name, 'bt2_paired_relaxed.err')}"""
            else:
                logging.info("Going into relaxed mode, single strand")
                bt_cmd = f"""bowtie2 -p {args.num_cores} --very-sensitive --end-to-end --no-unal \
                -x {args.bowtie_index} -U {read_filter.s1_path} -k 100 -S {sam_path} \
                2> {read_filter.get_path('logs', base_name, 'bt2_single_relaxed.err')}"""

        commands = [
            bt_cmd,
            f"samtools view -bS {sam_path} -o {bam_path}",
            f"samtools sort -@ {args.num_cores} {bam_path} -o {sorted_bam}",
            f"samtools index {sorted_bam}",
            # Multimap uses this intermediate output below in the final dir
            f"bedtools intersect -abam {indexed_bam} -b {args.herv_bed} -wa -wb -bed > {intersect_file}",
            f"bedtools intersect -abam {sorted_bam} -b {args.herv_bed} > {bed_out}",
            # Kmer step uses this intermediate output in the final dir
            f"samtools view -h -o {hervs_read} {bed_out}"
        ]

        logging.info("Starting processing pipeline...")
        for cmd in tqdm(commands, desc="Processing commands"):
            logging.debug(f"Executing: {cmd}")
            run_command(cmd)

        cleanup_dir = os.path.join(args.output_dir, 'intermediate_files')
        if not args.keep_files:
            logging.info("Cleaning up intermediate files...")
            if os.path.exists(cleanup_dir):
                for filename in os.listdir(cleanup_dir):
                    file_path = os.path.join(cleanup_dir, filename)
                    try:
                        if os.path.isfile(file_path):
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

    except Exception as e:
        logging.error(f"An error occurred: {str(e)}")
        raise


if __name__ == "__main__":
    main()
