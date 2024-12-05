# Python Package CLI setup w/ dependencies
from setuptools import setup, find_packages
import argparse
import subprocess
import os
import logging
import time
from tqdm import tqdm

setup(
    name="ervromancer",
    version="0.1.0"
    packages=find_packages(),
    install_requires=[
        "bowtie2",
        "samtools",
        "bedtools",
        "pysam"
    ],
    entry_points={
        "console_scripts": [
            "ervromancer=ervromancer.cli:main",
        ],
    },
)


def run_command(cmd):
    subprocess.run(cmd, shell=True, check=True)


def get_path(args, subdir, filename, ext):
    return os.path.join(args.output_dir, subdir, f"{filename}.{ext}")


def main():
    parser = argparse.ArgumentParser(
        description='Ervromancer')
    # ARGUMENTS
    parser.add_argument(
        '--mode', choices=['strict', 'relaxed'], required=True, help="Alignment bounds (stringent, relaxed)")
    parser.add_argument('--r1', required=False,
                        help='Path to paired-end R1 fastq file')
    parser.add_argument('--r2', required=False,
                        help='Path to paired-end R2 fastq file')
    parser.add_argument('--s1', required=False,
                        help='Path to paired-end R2 fastq file')
    parser.add_argument('--keep_files', default=False,
                        help="Keeps intermediate outputs from Samtools and Bedtools.")
    parser.add_argument('--output-dir', required=True, help='Output directory')

    args = parser.parse_args()
    mode = args.mode
    r1 = args.r1
    r2 = args.r2
    if r1 and r2:
        base_folder = os.path.splitext(os.path.basename(args.r1))[0]

    commands = [
        f"samtools view -bS {sam_path} -o {bam_path}",
        f"samtools sort -@ {args.threads} {bam_path} -o {sorted_bam}",
        f"samtools index {sorted_bam}",

        f"bedtools intersect -abam {sorted_bam} -b {args.herv_bed} -wa -wb -bed > {bed_out}"
    ]

    # for dir_name in ['sam', 'bam', 'sorted_bam', 'bed', 'herv_reads', 'final']:
    #    os.makedirs(os.path.join(args.output_dir, dir_name), exist_ok=True)

    # if mode == "relaxed":
    #    bt_cmd = f"""bowtie2 -p {args.threads} --very-sensitive --end-to-end -X 1500 \
    #        --no-mixed --no-discordant --no-dovetail --no-unal --score-min L,0,1.6 \
    #        -x {args.bowtie_index} -1 {args.r1} -2 {args.r2} -k 100 -S {sam_path}""",
    # else:

    # Pipeline commands

    progress_bar = tqdm(commands)
    for cmd in tqdm(progress_bar):
        time.sleep(1)
        run_command(cmd)


if __name__ == "__main__":
    main()
