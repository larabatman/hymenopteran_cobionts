#!/usr/bin/env python3
"""
build_og_gene_annotations.py

Step 3 of 3 in the OG annotation database build.
Run after download_og_alignments.sh and build_og_lookup_table.sh.

Annotates each Wolbachia orthogroup with gene product names derived from the Prokka annotations of the 1444 Wolbachia genomes used to build the HMMs.
The result is og_gene_annotations.tsv, which is joined alongside og_lookup_table.tsv to every NUWT hit during the annotation step.

The two annotation sources catch different things and have zero overlap in which OGs they flag as transposable-element related:
SwissProt blastp catches proteins where the Wolbachia protein by comparison to a eukaryotic database
With the Prokka annotations, we can catch proteins whose genes were potentially bacterial transposases or IS elements directly within the Wolbachia reference genomes.
Using both sources together gives more complete TE classification.


Requires:
nuwt_scan/hmm_database/NUWTs_Release_v3.0/Data/1_WolbachiaAnnotation/ Prokka FAA files (one per Wolbachia genome, gzipped)
nuwt_scan/hmm_database/NUWTs_Release_v3.0/Data/2_WolbachiaOrthoFinder/ Orthogroups.tsv.gz (OrthoFinder output mapping OG IDs to proteins)

Output columns:
og_id Orthogroup identifier (OG0000001)
n_members: Total protein members across all Wolbachia genomes
n_annotated: Members with a non-hypothetical product name
n_hypothetical: Members annotated as "hypothetical protein"
consensus_name: Most common non-hypothetical product name
consensus_pct: Fraction of annotated members with the consensus name
top_names: Top 3 product names with counts (semicolon-separated)
is_te_by_name: TRUE if consensus name matches TE-related keywords

Usage:
python3 -u scripts/nuwt/database/build_og_gene_annotations.py \\
    nuwt_scan/hmm_database/NUWTs_Release_v3.0/Data/1_WolbachiaAnnotation \\
    nuwt_scan/hmm_database/NUWTs_Release_v3.0/Data/2_WolbachiaOrthoFinder/Orthogroups.tsv.gz \\
    nuwt_scan/hmm_database/og_gene_annotations.tsv
"""

import sys
import os
import gzip
import re
from collections import defaultdict, Counter


def main():
    annot_dir = sys.argv[1]
    orthogroups_gz = sys.argv[2]
    output_tsv = sys.argv[3]
    # Build a dictionary of protein_id and product_name: Parse all Prokka FAA files from 1_WolbachiaAnnotation/. where each file represents one Wolbachia genome. 
    # FASTA header format: >PROTEIN_ID product name >GCA_902644825_00001 Glutamate--tRNA ligase 1
    # Map every protein ID to its product name across all  genomes. This is later used to look up names for each protein ID found in Orthogroups.tsv.
    # The goal is to have a dictionary of protein ID and their product name for every protein in every genome, so that when we look at the actual OGs from the Orthogroups table, we can look each one up there
    protein_names = {} # protein_id to product_name
    n_faa = 0

    # Only the Prokka FAA files, sorted for deterministic order.
    faa_files = sorted(f for f in os.listdir(annot_dir) if f.endswith(".faa.gz"))

    for faa_file in faa_files: #loop opens each gzpped FAA and walks line by line
        faa_path = os.path.join(annot_dir, faa_file)

        with gzip.open(faa_path, "rt") as f: # "rt" is read-text: gzip.open defaults to binary
            for line in f:
                if not line.startswith(">"): # only want the headers, skip the sequences
                    continue # sequence line, only headers needed
                parts = line[1:].rstrip().split(None, 1) # line[1:] drops the >, .rstrip() removes the newline and .split(None, 1) splits on whitespace once such that the ID and the entire remaining name are kept as ["Ace_000025", "ankyrin repeat protein"]
                if len(parts) < 2: # handle headers with IDs and nothing else, putting them in a hypothetical protein bucket
                    protein_names[parts[0]] = "hypothetical protein"
                else:
                    protein_names[parts[0]] = parts[1]
        n_faa += 1

    # Parse Orthogroups.tsv.gz
    # Format: tab-separated
    # Column 0: OG ID (e.g. OG0000001)
    # Columns 1..N: one column per Wolbachia genome, containing comma-separated protein IDs for that genome. Empty cell = that genome has no member in this OG
    # Row 0: header with genome prefix names (skipped)

    og_proteins = defaultdict(list)   # defaultdict(list) means accessing a new key creates an empty list: og_id to [protein_id, ...]
    with gzip.open(orthogroups_gz, "rt") as f:
        f.readline() # Skip the header row: the first row is genome names, not data
        for line in f:
            parts = line.rstrip("\n").split("\t") # strip only the newline, not all whitespace
            og_id = parts[0]
            if not og_id.startswith("OG"): # Skip any malformed or unexpected rows
                continue
            for cell in parts[1:]: # Walk every genome column, an empty cell means that genome contributes nothing to the OG, but a non-empty cell holds that genom's member proteins. Collect protein IDs from all genome columns (everything after the OG ID in column 0)
                cell = cell.strip()
                if not cell:
                    continue # Empty cell = this genome has no member in this OG
                for pid in cell.split(","): # Each cell is a comma-separated list of protein IDs "Ace_00025, Ace_00029, Ace_00041"
                    pid = pid.strip()
                    if pid:
                        og_proteins[og_id].append(pid)

    # Build consensus annotation
    # For each OG, collect all product names for its member proteins, separate informative names from "hypothetical protein", and take the majority-vote
    # TE keyword pattern: flag OGs whose consensus name contains keywords indicating a transposable element or mobile element.
    te_pattern = re.compile(
        r"transpos" # transposon, transposase, retrotransposon
        r"|reverse.transcriptase"
        r"|gag.pol" # gag-pol polyprotein (retrovirus/retrotransposon)
        r"|pol.polyprotein"
        r"|retrovir" # retrovirus, retroviral
        r"|\bLTR\b" # word boundary prevents matching e.g. "filter"
        r"|\bSINE\b"
        r"|\bLINE\b"
        r"|Ty[0-9]" # Ty1, Ty3: yeast/insect LTR retrotransposons
        r"|\bcopia\b"
        r"|\bgypsy\b"
        r"|mariner"
        r"|helitron"
        r"|piggyBac"
        r"|\bIS[0-9]\b" # IS1, IS3, IS4 etc: bacterial insertion sequences
        r"|integrase"
        r"|resolvase", # site-specific recombinases used by IS elements
        re.IGNORECASE
    )

    # Pattern for uninformative product names, excluded before computing the consensus
    hypo_pattern = re.compile(
        r"hypothetical"
        r"|putative protein$"
        r"|unknown function"
        r"|uncharacterized protein",
        re.IGNORECASE
    )

    rows = []
    for og_id, proteins in sorted(og_proteins.items()):
        n_members = len(proteins)
        names = [] # Look up each protein's product name in the dictionary we built from the FAA files
        for pid in proteins:
            name = protein_names.get(pid)
            names.append(name) # name is a list with one entry per member protein, so for a family of 900 members, 900 strings
        # Split the OG member names
        n_hypothetical = sum(1 for n in names if hypo_pattern.search(n)) # produces number 1 for each name that matches, summed across
        informative = [n for n in names if not hypo_pattern.search(n)] # same tested but inverted
        n_annotated = len(informative)

        if n_annotated == 0: # case 1: all members are hypothetical
            consensus_name = "hypothetical protein"
            consensus_pct  = 0.0
            top_names_str  = "hypothetical protein"
        else: # case 2: majority vote, run when at least one member had an infromative name
            name_counts = Counter(informative) # Counter is a dictionnary subclass that, from a list, returns each distinct value mapped to how many times it appeared: {"IS110 family transposase": 240, "ecombinase": 45} and so on
            top3 = name_counts.most_common(3) # Returns the top three as a list of (name, count) pairs, ordered by count descending [("ISS110 family transposase", 240), ("recombinase", 45), ("integrase", 15)]
            consensus_name = top3[0][0] #the first element of the first pair is the most frequenc informative name and becomes the consensus
            consensus_pct  = top3[0][1] / n_annotated # consensus_pct: fraction of informative (non-hypothetical) members that have the consensus name. A value of 1.0 means all annotated members agree on this name.

            # Summary of top names with their counts
            top_names_str  = "; ".join(
                f"{name} ({count})" for name, count in top3
            )

        # Flag as TE if the consensus name matches TE keywords.
        is_te_by_name = bool(te_pattern.search(consensus_name))
        rows.append({
            "og_id": og_id,
            "n_members": n_members,
            "n_annotated": n_annotated,
            "n_hypothetical": n_hypothetical,
            "consensus_name": consensus_name,
            "consensus_pct": round(consensus_pct, 3),
            "top_names": top_names_str,
            "is_te_by_name": is_te_by_name,
        })

    #TSV output
    cols = ["og_id", "n_members", "n_annotated", "n_hypothetical","consensus_name", "consensus_pct", "top_names", "is_te_by_name"]
    with open(output_tsv, "w") as out:
        out.write("\t".join(cols) + "\n")
        for row in rows:
            out.write("\t".join(str(row[c]) for c in cols) + "\n")