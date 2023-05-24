import pysam
import argparse


def count_mutated_reads(bam_file: str, min_mapping_quality: int) -> int:
    """
    Counts the number of unpaired reads or read pairs that have two or more mutations in a BAM file.

    Args:
        bam_file (str): Path to the BAM file.
        min_mapping_quality (int): Minimum mapping quality allowed for each read.

    Returns:
    int: The count of unpaired reads or read pairs with two or more mutations.
    """
    samfile = pysam.AlignmentFile(bam_file, "rb")
    mutated_reads_count = 0
    for read in samfile.fetch():
        if read.mapping_quality < min_mapping_quality:
            continue
        if not read.is_paired:
            # NM tag ( edit distance) >= 2
            if read.get_tag("NM") >= 2:
                mutated_reads_count += 1
        else:
            # For paired-end reads , increment the counter only if this is the first in pair (
            # to avoid double-counting)
            if read.is_read1 and read.get_tag("NM") >= 2:
                mutated_reads_count += 1
    return mutated_reads_count

def main():
    parser = argparse.ArgumentParser(description='This tool parses a bam file and prints the number of reads which have 2 or more mutations. The reads are assumed to be aligned to the same genomic location.')
    parser.add_argument('-v', dest='verbose', action='store_true', help='print verbose info')
    parser.add_argument('-i', dest='infiles', nargs='+', type=str, help='input bam files')
    parser.add_argument('-q', dest='quality', type=int, help='minimum mapping quality')
    
    # parse arguments
    args = parser.parse_args()

    # adding ability to parse and count through multiple bam files
    out = []
    for f in args.infiles:
        # add exception handling in case the function fails to execute on a particular file for any reason
        try:
            x = count_mutated_reads(f, args.quality)
            out.append(x)
        except:
            raise ValueError('Failed to process',f)
    print(sum(out))

main()
