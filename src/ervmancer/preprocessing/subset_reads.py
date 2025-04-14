import gzip
from itertools import zip_longest


def open_file(filename):
    '''
    from chatGPT, for opening the original file names
    Specifically, this allows handling either gz or normal files
    '''
    # Check if the file is gzipped based on the extension
    if filename.endswith('.gz'):
        # Open the file with gzip if it is gzipped
        return gzip.open(filename, 'rt')
    else:
        # Open the file normally if it is not gzipped
        return open(filename, 'r')


def process_fastq_chunk(fastq_file, kept_reads):
    """Process FASTQ file in chunks of 4 lines and yield ID and sequence."""
    while True:
        id_line = fastq_file.readline()
        if not id_line:
            break

        seq_line = fastq_file.readline()
        # Skip the next two lines (+ and quality scores)
        fastq_file.readline()
        fastq_file.readline()

        # Extract read ID without @ and /1 or /2 suffix
        read_id = id_line.split('/')[0].strip('@')

        if read_id in kept_reads:
            yield read_id, seq_line


def extract_out_original_reads_by_subset(pathname_r1,
                                         pathname_subset_sam_for_filtering,  # sam file subset from step 5
                                         pathname_extracted_reads,  # output path
                                         pathname_r2=''):
    """
    Extract out the original reads based on the subsetted sam file.
    This is done because we want the actual reads for the kmer step.
    Bowtie2 modifies the reads, including softclipping, etc.
    This returns a tab delineated file which contains the read id and the sequence.

    Parameters:
        pathname_r1 (str): The path and name of the r1 file
        pathname_subset_sam_for_filtering (str): The pathname of the sam file output at the final step of the intersection step
        pathname_extracted_reads (str): The pathname where extracted reads will be written
        pathname_r2 (str, optional): The path and name of the r2 file for paired-end reads. Defaults to ''.
    """
    # Extract read IDs to keep from the sam file
    with open(pathname_subset_sam_for_filtering, 'r') as subset_sam:
        # Skip header lines starting with @
        for line in subset_sam:
            if not line.startswith('@'):
                break

        only_keep_these_reads = {line.split('\t')[0]}

        # Add remaining read IDs
        only_keep_these_reads.update(
            line.split('\t')[0] for line in subset_sam)

    # Process the reads and write to output file
    with open(pathname_extracted_reads, 'w') as out_file:
        # Handle paired-end reads
        if pathname_r2:
            with open_file(pathname_r1) as in_r1, open_file(pathname_r2) as in_r2:
                for (r1_id, r1_seq), (r2_id, r2_seq) in zip_longest(
                        process_fastq_chunk(in_r1, only_keep_these_reads),
                        process_fastq_chunk(in_r2, only_keep_these_reads)):
                    out_file.write(f"{r1_id}\t{r1_seq}")
                    out_file.write(f"{r2_id}\t{r2_seq}")
        # Handle single-end readss
        else:
            with open_file(pathname_r1) as in_r1:
                for read_id, seq in process_fastq_chunk(in_r1, only_keep_these_reads):
                    out_file.write(f"{read_id}\t{seq}")
