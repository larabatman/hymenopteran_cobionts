#!/usr/bin/env python3
"""
extract_nuwt_seqs.py

Reads all per-species filtered dfam tables and extracts the corresponding nuleotide sequences from each host genome assembly using samtools faidx.
Writes one FASTA per orthogroup, pooling sequences from every species that had hits to that OG. These per-OG FASTAs are the input for the placement step.

The dfam tables are organised by species: placement needs one tree per OG, so this step is the reorganisation from one to the other.

This is an adaptation of FeatchRegionSeq.py from Vancaester et al., that performs the same task: read hit coordinate, extract the corresponding sequence from the host genome, write FASTA with a header naming the species, scaffol and orthogroup.
Here the same header format is used with an added INS_ prefix so that we can reuse their Reroot_rename.py script
Changes that were made: 
Extraction method: instead of SeqIO.parse, we already had the .fai indexes from the filter step so samtools faidx is used directly
Input: reading the contamination-filters --dfamtblout table directly
Coordinates: here we use env_from/env_to instead of ali_start/ali_stop which is a wider estimate of the region that include weak edges. Everything is 1-based inclusive, like samtools expects
Output organisation: instead of one FASTA per species, we do one FASTA per orthogroup pooled across all species as each OG is placed in its own tree later

Coordinates: env_from/env_to (the envelope, wider than the alignment) are used, and they are 1-based inclusive, the same convention samtools faidx expects

Header format:
>INS_<Species>_<scaffold>_<OG>:<start>-<stop>
The INS_ prefix lets Reroot_rename_tree.py tell insect NUWT sequences from Wolbachia reference sequences in the tree.

A fragment table is written alongside the FASTAs: one row per extracted fragment,
carrying the same fields as the header but in separate columns. The header has to
pack species, scaffold, OG and coordinates into one string, and both species and
scaffold names contain "_", so taking that string apart again downstream needs
guesswork. The table removes that: later steps join on `name` instead of parsing it.

Usage:
    python3 extract_nuwt_seqs.py 
"""

import os
import glob
import subprocess
from collections import defaultdict

COBIONTS_ROOT = "/data/users/lland/cobionts"
NUWT_DIR = os.path.join(COBIONTS_ROOT, "nuwt_scan")
RESULTS_DIR = os.path.join(NUWT_DIR, "results")
ASSEMBLY_DIR = os.path.join(NUWT_DIR, "host_assemblies")
OG_SEQ_DIR = os.path.join(NUWT_DIR, "og_sequences")

os.makedirs(OG_SEQ_DIR, exist_ok=True)

# Define the positions in the dfam line
# Dfam column indices (0-based), matching annotate_nuwt_hits.R: 0:og_id 1:acc 2:scaffold 3:score 4:evalue 5:bias 6:hmm_from 7:hmm_to 8:strand 9:ali_from 10:ali_to 11:env_from 12:env_to 13:modlen 14:description
COL_OG = 0
COL_SCAFFOLD = 2
COL_STRAND = 8
COL_ENV_FROM = 11
COL_ENV_TO = 12

# find any .fna and drop anything ending with .fai, sort, return the first or None
def find_assembly(species):
    """Return the assembly path for a species, or None if there isn't one.

    Assemblies are named <Species>_<Accession>.fna. The .fai index sits beside them and matches the same glob, so it is filtered out.
    """
    pattern = os.path.join(ASSEMBLY_DIR, f"{species}_*.fna")
    matches = sorted(f for f in glob.glob(pattern) if not f.endswith(".fai"))
    return matches[0] if matches else None

# Read dfam line by line, skip # headers and blanks, line.split() on whitespace and pull the five fileds on a dictionary
def parse_dfam(filepath):
    """Return one dict per hit from a dfam table, skipping "#" header lines."""
    hits = []
    with open(filepath) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.split()
            hits.append({
                "og": cols[COL_OG],
                "scaffold": cols[COL_SCAFFOLD],
                "strand": cols[COL_STRAND],
                "env_from": int(cols[COL_ENV_FROM]), # int() needed for maths below
                "env_to": int(cols[COL_ENV_TO]),
            })
    return hits

# parse_dfam turns dfam line into: {"og": "OG0001339", "scaffold": "OZ125706.1", "strand": "+", "env_from": 15060639, "env_to": 15061970}
# we know which OG, which scaffold, where, on which strand, but we don't have the DNA
# this needs to be converted into an address that samtools can pick up: regions.append() to produce a string "OZ125706.1:15060639-15061970"
# and we can store in hit_meta.append the OG information: ("OG0001339", "OZ125706.1", 15060639, 15061970)
# to ask for the DNA: cmd with the argument list, samtools run and out with the text: >OZ125706.1:15060639-15061970 ATGGATAGCATAGATAGCATATACAG
# and re-attach the OG by walking the text. when record 0 finished, it calls save(0, seq_lines) to look up hit_meta[0] and rebuild the proper: >INS_Abia_candens_OZ125706.1_OG0001339:15060639-15061970
def extract_sequences(species, assembly, hits):
    """Return (og, header, sequence, scaffold, strand, start, stop) for every hit in one species.

    Two samtools calls, one per strand: -i reverse-complements everything a call returns, so plus and minus hits cannot share one. 
    Each call passes all of that strand's regions at once rather than one subprocess per hit.
    The coordinate fields are returned as well as packed into the header, so the
    fragment table can be written without parsing the header back apart.
    """
    results = [] # Accumulates the return value across both strands. Each query is a tuple 3: (og, header, sequence)

    for strand, rc_args in [("+", []), ("-", ["-i"])]: # two element list unpacked into two variables. Iteration 1: strand="+", rc_args=[] and iteration 2: strand"-", rc_args=["-i"]
        hit_group = [h for h in hits if h["strand"] == strand] #iteration subset: out of 500 hits, 260 are + and 240 are -
        if not hit_group:
            continue
        regions, hit_meta = [], [] #Two empty lists. regions is what samtools sees, and hit_meta is what we need: samtools output header will say OZ125706.1:15060639-15061970 but the OG needs to be remembered by position
        for h in hit_group:
            if strand == "+": # start and stop: one hit at a time. HMMER writes minus-strand hits descending, like 5070033 - 5065624 and samtools requires acsending. -i reverse-complements what comes back
                start, stop = h["env_from"], h["env_to"]
            else:
                start, stop = h["env_to"], h["env_from"]
            regions.append(f"{h['scaffold']}:{start}-{stop}") #entry i must fo into both lists so they are kept aligned
            hit_meta.append((h["og"], h["scaffold"], start, stop))

        cmd = ["samtools", "faidx"] + rc_args + [assembly] + regions #one command for all regions: argument vector ["samtools", "faidx", "-i", "/path/genome.fna", "scafA:1-100", ...]
        out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout # out: what samtools printed with text=True to have str instead of bytes
        idx, seq_lines = -1, [] # the parser state: idx is the region we are currently reading, starting at -1 so the first > changes it to 0 and seq_lines collects the wrapped sequence lines of the current record

        def save(idx, seq_lines): # to recover the OG: hit_meta[idx] if the lookup to recover OG and "".join(seq_lines) glues the wrapped lines into one sequence
            og, scaffold, start, stop = hit_meta[idx]
            header = f"INS_{species}_{scaffold}_{og}:{start}-{stop}"
            results.append((og, header, "".join(seq_lines).upper(),
                            scaffold, strand, start, stop))

        for line in out.splitlines():
            if line.startswith(">"): # header line: finish the previous record first, then advance and reset
                if idx >= 0 and seq_lines:
                    save(idx, seq_lines)
                idx += 1
                seq_lines = []
            else:
                seq_lines.append(line.strip())
        if idx >= 0 and seq_lines:# the last record has no > after it. Save the last record
            save(idx, seq_lines)
    return results


# Main
dfam_files = sorted(glob.glob(os.path.join(RESULTS_DIR, "*", "filter", "*_nuwt_hits_dfam_filtered.tbl")))

n_done = n_seqs_total = 0
og_counts = defaultdict(int)

# Keep a fragment table for merging after placement: one row per extracted fragment. Opened once and appended to as the species loop turns. Name is the FASTA header
frag_table = os.path.join(NUWT_DIR, "nuwt_fragments.tsv")
frag_fh = open(frag_table, "w")
frag_fh.write("name\tspecies\tscaffold\tog\tstrand\tstart\tend\tlength\n")

for dfam_file in dfam_files:
    species = os.path.basename(os.path.dirname(os.path.dirname(dfam_file))) # .../results/<Species>/filter/<Species>_nuwt_hits_dfam_filtered.tbl to recover the species name from the directory structure rather than by parsing the filename
    assembly = find_assembly(species) # find DNA source
    hits = parse_dfam(dfam_file) # find hit source
    if not hits: # case where there are not hits in the filtered dfam
        n_done += 1
        continue
    seqs = extract_sequences(species, assembly, hits) # samttools for the coordinates, returning (og, header, seq) tuples
    for og, header, seq, scaffold, strand, start, stop in seqs: # the destination file is chosen for each OG rather than for each species: one species sequences scattoer across different OG files, and each OG file received contributions from many species as the outer loop turns
        with open(os.path.join(OG_SEQ_DIR, f"{og}.fa"), "a") as out_fh:
            out_fh.write(f">{header}\n{seq}\n")
        og_counts[og] += 1
        # same fields as the header, one per column. start/stop are already
        # ascending here, so length is stop - start + 1 (1-based inclusive).
        frag_fh.write(f"{header}\t{species}\t{scaffold}\t{og}\t{strand}\t"
                      f"{start}\t{stop}\t{stop - start + 1}\n")
    n_seqs_total += len(seqs)
    n_done += 1

frag_fh.close()

# Build manifest
# The manifest is OG id + sequence count, one per line. filter_hit_ogs.sh reads it as its OG list, and its line count sizes the placement job array.
manifest = os.path.join(NUWT_DIR, "og_hit_manifest.txt")
with open(manifest, "w") as mf:
    for og in sorted(og_counts): # iterating a dict gives its keys: sorted(og_counts) is the OG IDs in order
        mf.write(f"{og}\t{og_counts[og]}\n")