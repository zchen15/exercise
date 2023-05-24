import pysam


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
