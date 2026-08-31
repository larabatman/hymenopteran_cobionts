#!/usr/bin/env python3
"""
This script selects one representative protein sequence for a given Wolbachia orthogroup. 
The input is an aligned fasta, one file for each OG. The protein sequences contain gap characters and may go through severale lines because of the alignment
The output is a fasta with one chosen protein for each OG, and the header >OG000043, >OG0000044, etc

Usage: 
python3 1_extract_of_representatives.py alignment_directory output_fasta
"""

import gzip
import re
import sys
from pathlib import Path
from Bio import SeqIO

# Kepts sequences are at least 10 aa long when the gaps are removed
MIN_LEN = 10

# Function to go through the alignment files in the deposit
def find_alignment_files(aln_dir):
    """
    Returns the OG alignment files in aln_dir, in sorted order.
    """
    return sorted(aln_dir.glob("*.aln.fa.gz"))

# Function to chose the least gapped protein in the alignment
# Each alignment file holds the same protein from many Wolbachia strains, all with - characters so that equivalent positions line up in columns
# MKTLLIGAGGSA---GFTQLAAKHV
# MKT--IGAGG-----GFTQ-----
# Decide which protein is going to represent the OG for the blastp
def least_gapped(sequences):
    """
    Return the sequence with the smallest fraction of gap characters.
    """
    # best_seq and best_gap_frac are the best seen so far
    # Find minimum: go through list, and whenever something is better than best seen this far, it becomes the new best
    best_seq = None
    # Starting somewhere that guarantees that the first seen is better than nothing
    best_gap_frac = 2.0
    for seq in sequences:
        if not seq:
            # Skip empty strings, division below
            continue 
        # Count the - characters and divide by total length to get the fraction
        gap_frac = seq.count("-") / len(seq)
        # If tie, first seen wins with the strictly less than
        if gap_frac < best_gap_frac:
            best_gap_frac = gap_frac
            best_seq = seq
    return best_seq

def main():
    aln_dir = Path(sys.argv[1])
    out_fasta = Path(sys.argv[2])

    aln_files = find_alignment_files(aln_dir)
    # Counter
    n_written = 0
    with open(out_fasta, "w") as out:
        for aln_file in aln_files:
            # Get OG id r"....": OG, then one or more digits. re.search to scan anywhere in the string for that pattern
            match = re.search(r"(OG\d+)", aln_file.name) 
            if not match:
                # Move to next file
                continue

            og_id = match.group(1)
            # Decompress the alignment fasta, read text mode 
            with gzip.open(aln_file, "rt") as handle: 
                # SeqIO.parse gives one SeqRecord for each fasta entry, one for each Wolbachia strain in that orthogroup
                # The sequence lines between two > headers get joined into one string
                # sequences: list of strings with every strain copy of that one protein
                sequences = [str(record.seq) for record in SeqIO.parse(handle, "fasta")]
            # Get the least gapped one of those sequences
            best_seq = least_gapped(sequences)
            # Skip it if there is none
            if best_seq is None:
                continue
            # Remove the - characters
            ungapped = best_seq.replace("-", "")
            # Make sure it is at least 10 aa for blastp
            if len(ungapped) < MIN_LEN:
                continue
            out.write(f">{og_id}\n{ungapped}\n")
            n_written += 1

    # Exit non-zero if nothing was written. An empty output file would otherwise pass silently to the blastp step
    if n_written == 0:
        sys.exit("[ERROR] No representative sequences written")

if __name__ == "__main__":
    main()