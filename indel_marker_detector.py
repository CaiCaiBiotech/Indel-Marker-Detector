import pysam
import sys
import re
import csv
import argparse
from itertools import groupby

def read_indel_file(file_path, min_indel_length=10):
    indels = set()
    vcf_in = pysam.VariantFile(file_path, "r")
    for rec in vcf_in.fetch():
        for alt in rec.alts:
            indel_length = abs(len(alt) - len(rec.ref))
            if indel_length >= min_indel_length:
                indel_seq = alt if len(alt) > len(rec.ref) else rec.ref
                indels.add((rec.chrom, rec.pos, indel_seq))
    return indels

def get_unique_indel_list(strain1_indels, strain2_indels, genome_file, marker_length_min, marker_length_max):
    unique_indels = strain1_indels.symmetric_difference(strain2_indels)
    genome = pysam.FastaFile(genome_file)

    marker_list = []
    for indel_chrom, indel_pos, indel_seq in unique_indels:
        # Skip if indel sequence is a palindrome or has a repeated subsequence occupying more than 50% of the sequence
        if is_palindrome(indel_seq) or is_simple_repeat(indel_seq):
            continue

        found_marker = False
        for marker_length in range(marker_length_min, marker_length_max + 1):
            start = indel_pos - marker_length // 2
            end = indel_pos + marker_length // 2
            
            # Ensure start is not negative
            start = max(0, start)

            seq = genome.fetch(region=indel_chrom, start=start, end=end)

            if 'N' in seq:  # Skip sequences with 'N'
                continue

            if not is_palindrome(seq) and not is_simple_repeat(seq):
                strain_with_indel = "strain1" if (indel_chrom, indel_pos, indel_seq) in strain1_indels else "strain2"
                marker_list.append((indel_chrom, start, end, seq, strain_with_indel, indel_seq))
                found_marker = True
                break
        if not found_marker:
            print(f"No suitable marker found for indel at {indel_chrom}:{indel_pos}")
    return marker_list

def is_palindrome(seq):
    seq_len = len(seq)
    for i in range(seq_len // 2):
        if seq[i] != seq[seq_len - 1 - i]:
            return False
    return True

def is_simple_repeat(seq):
    seq_len = len(seq)
    min_subseq_len = 3  # Minimum length of a subsequence to be considered
    max_same_char_stretch = max(len(list(g)) for _, g in groupby(seq))  # Length of the longest stretch of the same character
    if max_same_char_stretch > seq_len * 0.5:  # If the longest stretch of the same character occupies more than 50% of the sequence
        return True

    for repeat_len in range(min_subseq_len, seq_len // 2 + 1):
        # Check for repeated subsequences
        for i in range(seq_len - repeat_len + 1):
            subseq = seq[i:i+repeat_len]
            if seq.count(subseq) * len(subseq) > seq_len * 0.5:  # If the repeated subsequence occupies more than 50% of the sequence
                return True
    return False

def main(indel_files1, indel_files2, genome_file):
    strain1_indels = set.intersection(*(read_indel_file(file) for file in indel_files1))
    strain2_indels = set.intersection(*(read_indel_file(file) for file in indel_files2))
    marker_list = get_unique_indel_list(strain1_indels, strain2_indels, genome_file, 300, 400)

    strain1_markers = []
    strain2_markers = []
    for marker in marker_list:
        chrom, start, end, seq, strain_with_indel, indel_seq = marker  # Corrected here
        indel_pos = (start + end) // 2
        marker_info = (chrom, start, end, indel_pos, seq, indel_seq)  # Corrected here
        if strain_with_indel == "strain1":
            strain1_markers.append(marker_info)
        else:
            strain2_markers.append(marker_info)

    # Save strain1 markers to a CSV file
    with open('strain1_markers.csv', 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(['Chromosome', 'Marker Start', 'Marker End', 'Indel Position', 'Marker Sequence', 'Indel Sequence'])
        for marker in strain1_markers:
            csvwriter.writerow(marker)

    # Save strain2 markers to a CSV file
    with open('strain2_markers.csv', 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(['Chromosome', 'Marker Start', 'Marker End', 'Indel Position', 'Marker Sequence', 'Indel Sequence'])
        for marker in strain2_markers:
            csvwriter.writerow(marker)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='This script detects unique indel markers from multiple VCF files for two different strains.')
    parser.add_argument("--strain1", nargs='+', required=True, help="Specify one or more VCF files for strain 1.")
    parser.add_argument("--strain2", nargs='+', required=True, help="Specify one or more VCF files for strain 2.")
    parser.add_argument("--genome", required=True, help="Specify the genome file.")
    args = parser.parse_args()

    main(args.strain1, args.strain2, args.genome)