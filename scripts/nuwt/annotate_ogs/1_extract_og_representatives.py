#!/usr/bin/env python3
"""
Extract one representative protein sequence per Wolbachia orthogroup.

Called by download_og_alignments.sh after the Zenodo archive has been unpacked

Usage:
python3 extract_og_representatives.py <alignment_dir> <output_fasta>

Input : aligned FASTA, protein (amino acid), one file per OG, optionally gzipped. Sequences contain '-' gap characters and may be wrapped across several lines.
Output: unaligned FASTA, protein, one record per OG, header ">OG#######".

The only transformation is the removal of alignment gap characters, which are an artefact of the multiple sequence alignment and not part of the protein.

Output extension note: 
The pipeline calls the output og_representatives.fna, but it holds protein sequence and by convention should be .faa. The name is kept for consistency with the rest of the pipeline
"""

import gzip
import re
import sys
from pathlib import Path

# Minimum ungapped length, in amino acids, for a sequence to be kept.
MIN_LEN = 10

# Alignment file extensions found in the Zenodo archive
# This defines which files we are going to loop over for the main iteration
# For each alignement file, pick a representative
def find_alignment_files(aln_dir):
    """Return the OG alignment files in aln_dir, in sorted order.
    All 3,786 files in the Zenodo archive are named OG#######.aln.fa.gz, so a single pattern is enough. Sorting makes the output file byte-identical across runs, since directory listing order is not guaranteed.
    """
    return sorted(aln_dir.glob("*.aln.fa.gz")) #aln_dir.glob() asks the filesystem for every file in that dir matching the pattern, returning 3786 PATH objects, one per OG. Sorted to have them in alphabetical order

# FASTA: a record is a header line starting with > with the sequence spread across as many lines as it takes
#>wAlb_00421
#MKTLLIGAGGSA----GFTQLAAKHNVELV
#LIDRNRDR----IPYFTGD
# The example record contains one protein, but might look like two separate things
# We need to stitch them back together to have the complete record
# Using yield instead of return makes it a generator that hands out records one at a time as they are asked for, rather than building a list of all of them first
# This matters to read and emit as the function goes rather than hold everything in memory
def parse_fasta(handle):
    """
    Yield (header, sequence) pairs from an open FASTA file handle.
    Sequence lines between two headers are joined.
    """
    # header and parts are the accumulator: header holds the name of the record currently being buit, and parts collects its sequence lines.
    header = None
    parts = []
    for line in handle: 
        line = line.strip() # for each line,strip() removes the trailing newline and any stray whitespace
        if line.startswith(">"): # case 1: the line starts with >: new record begins. Hand out the preivous record
            if header is not None: # is not None check ecists for the first header in the file, where there is no previous record to emit
                yield header, "".join(parts) # Glue the collected sequence lines into one string
            header = line[1:] # store the name minus the > 
            parts = [] # Empty the accumulator for the new record
        elif line: # case 2: the line is anything else, true for any non-empty line: it is sequence and gets appened to parts. Empyt lines are silently skipped
            parts.append(line)
    if header is not None:          # the final record has no header after it
        yield header, "".join(parts)

# Each alignment file holds the same protein from many Wolbachia strains, all padded with - characters so that equivalent positions line up in columns
#MKTLLIGAGGSA---GFTQLAAKHV
#MKT--IGAGG-----GFTQ-----
# the first has fewer gaps than the second. We need a function that decides which becomes the blastp query for that OG 
def least_gapped(sequences):
    """Return the sequence with the smallest fraction of gap characters.

    A low gap fraction means this taxon is represented across most alignment columns, so its ungapped sequence is the longest and most complete version of the protein available in the OG and the best for query.
    Returns None if no non-empty sequence is present.
    """
    # best_seq and best_seq_frac are the running winner. Find the minimum pattern: walk the list, and whenever something beats the current best, it becomes the new best
    best_seq = None
    best_gap_frac = 2.0 # sentinel: a fraction can never exceed 1.0, which would be a sequence of nothing but dashes. Starting at 2.0 guarantees the first real sequence always beats it and gets installed as the winner
    for seq in sequences:
        if not seq:
            continue # skip empty strings to abandon iteration because of the division by len(seq) below
        gap_frac = seq.count("-") / len(seq) # counts dash characters and divides by total length. A sequence that has 30- out of 100 characters scores 0.3
        if gap_frac < best_gap_frac: # strictly less than: on a tie, the earlier sequence keeps the title 
            best_gap_frac = gap_frac
            best_seq = seq
    return best_seq # hands the winning string, still gapped 


def main():
    aln_dir = Path(sys.argv[1])
    out_fasta = Path(sys.argv[2])

    aln_files = find_alignment_files(aln_dir)
    # Define file cases to see if everything adds up
    n_written = 0
    n_skipped_short = 0  # best ungapped sequence below MIN_LEN


    with open(out_fasta, "w") as out:
        for aln_file in aln_files:
            # Get the OG id from the filename with a regex rather than by splitting on '.'
            match = re.search(r"(OG\d+)", aln_file.name) #r"....": the literal letters OG, then one or more digits. re.search to scan anywhere in the string for that pattern
            if not match: #re.search retuns None when nothing matches, continue abandons the file and moves to the next one
                print(f"[WARN] No OG id in filename: {aln_file.name}", file=sys.stderr)
                n_skipped_parse += 1
                continue
            og_id = match.group(1) #text captured by the first parenthsis, wrapping the entire pattern

            with gzip.open(aln_file, "rt") as handle: # open the compressed file and decompress it on the fly, 
                sequences = [seq for _header, seq in parse_fasta(handle)] # list comprehension: parse_fasta(handle) yiled (header, sequence) pairs one at a time and for _header, seq in unpacks each pair into its two parts, seq is kept
                # result: plain list of sequence strings for that one file, every strain's copy of that protein expcted by least_gapped()

            best_seq = least_gapped(sequences) # get the single winner as the row eith the smallest gap fraction: the longest real protein in the OG

            ungapped = best_seq.replace("-", "") # removes every - from the string by replacing it with nothing, the moment an alignment stop being an ailgnment and becomes a protein sequence
            if len(ungapped) < MIN_LEN: # length measured after stripping
                n_skipped_short += 1
                continue

            out.write(f">{og_id}\n{ungapped}\n")
            n_written += 1

    print(f"[INFO] Written: {n_written} sequences", file=sys.stderr)
    print(f"[INFO] Skipped (short) : {n_skipped_short} "
          f"(ungapped length < {MIN_LEN} aa)", file=sys.stderr)
    print(f"[INFO] Output: {out_fasta}", file=sys.stderr)

    # Exit non-zero if nothing was written. An empty output file would otherwise pass silently to the blastp step
    if n_written == 0:
        sys.exit("[ERROR] No representative sequences written")


if __name__ == "__main__":
    main()