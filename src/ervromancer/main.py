import argparse
import logging
import os
import subprocess
from tqdm import tqdm
from ervromancer.preprocessing.filter_reads import ReadFilter


def run_command(cmd):
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed: {cmd}")
        raise e


def main():
    parser = argparse.ArgumentParser(description='Ervromancer')
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
    parser.add_argument('--herv_bed', required=True,
                        help="Absolute path to the HERV BED file for intersection.")

    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s')

    try:
        read_filter = ReadFilter(args.output_dir, args.r1, args.r2, args.s1)
        base_name, paired = read_filter.validate_inputs()

        os.makedirs(os.path.join(args.output_dir,
                    "intermediate_files"), exist_ok=True)

        sam_path = read_filter.get_path(
            'intermediate_files', 'all_hervs.bwa.read', 'sam')
        bam_path = read_filter.get_path(
            'intermediate_files', 'all_hervs.bwa.read', 'bam')
        sorted_bam = read_filter.get_path(
            'intermediate_files', 'all_hervs.bwa.read.sorted', 'bam')
        bed_out = read_filter.get_path(
            'intermediate_files', 'all_hervs.bwa.read', 'bed')
        hervs_read = read_filter.get_path(
            'intermediate_files', base_name, 'bam')

        if args.mode == "strict":
            if paired:
                bt_cmd = f"""bowtie2 -p {args.num_cores} --very-sensitive --end-to-end -X 1500 \
                --no-mixed --no-discordant --no-dovetail --no-unal --score-min L,0,1.6 \
                -x {args.bowtie_index} -1 {read_filter.r1_path} -2 {read_filter.r2_path} -k 100 -S {sam_path} \
                2> {read_filter.get_path('final', base_name, 'bt2_paired_strict.err')}"""
            else:
                bt_cmd = f"""bowtie2 -p {args.num_cores} -N 1 -L 10 --very-sensitive-local --end-to-end --no-unal \
                -x {args.bowtie_index} --score-min L,0,1.6 -U {read_filter.s1_path} -k 100 -S {sam_path} \
                2> {read_filter.get_path('final', base_name, 'bt2_single_strict.err')}"""
        else:  # relaxed mode
            if paired:
                bt_cmd = f"""bowtie2 -p {args.num_cores} --sensitive-local --no-mixed --no-discordant \
                -x {args.bowtie_index} -1 {read_filter.r1_path} -2 {read_filter.r2_path} -k 100 -S {sam_path} \
                2> {read_filter.get_path('final', base_name, 'bt2_paired_relaxed.err')}"""
            else:
                bt_cmd = f"""bowtie2 -p {args.num_cores} --sensitive-local \
                -x {args.bowtie_index} -U {read_filter.s1_path} -k 100 -S {sam_path} \
                2> {read_filter.get_path('final', base_name, 'bt2_single_relaxed.err')}"""

        commands = [
            bt_cmd,
            f"samtools view -bS {sam_path} -o {bam_path}",
            f"samtools sort -@ {args.num_cores} {bam_path} -o {sorted_bam}",
            f"samtools index {sorted_bam}",
            # Multimap uses this intermediate output below
            f"bedtools intersect -abam {sorted_bam} -b {args.herv_bed} -wa -wb -bed > {bed_out}",
            # TODO: refactor this (step 5)
            f"bedtools intersect -abam {sorted_bam} -b {args.herv_bed} > {bed_out}",
            # Kmer step uses this intermediate output
            f"samtools view -h {args.herv_bed} {sorted_bam} -o {hervs_read}"
        ]

        logging.info("Starting processing pipeline...")
        for cmd in tqdm(commands, desc="Processing commands"):
            logging.debug(f"Executing: {cmd}")
            run_command(cmd)

        if not args.keep_files:
            logging.info("Cleaning up intermediate files...")
            for path in [sam_path, bam_path]:
                if os.path.exists(path):
                    os.remove(path)

        logging.info("Processing complete!")

    except Exception as e:
        logging.error(f"An error occurred: {str(e)}")
        raise


if __name__ == "__main__":
    main()
