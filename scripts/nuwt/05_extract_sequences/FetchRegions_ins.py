#!/usr/bin/env python3
"""
Adapted from Vancaester et al. 2025 FetchRegionSeq.py
Here: input is the contamination filtered dfam tableand the nhmmscan already applies the cutoff on E-value
We also used the env_from env_to coordinates instead of the ali_start ali_stop ones
The header prefix is INS.
The output here is appended to one fasta for each orthogroup: across species, each contributes its fragment that got a hit to a given orthogroup to place in its tree

Usage:
python3 FetchRegions_ins.py filtered_dfam_table genome_fasta OG_fasta_dir final_fasta
"""
import argparse
import os
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str, dest="input", help="filtered dfam table")
parser.add_argument("-f", type=str, dest="fasta", help="host genome fasta")
parser.add_argument("-d", type=str, dest="outdir", help="OG fasta directory")
parser.add_argument("-t", type=str, dest="table", help="fragment table")
results = parser.parse_args()

os.makedirs(results.outdir, exist_ok=True)
# dfam cols: 0 og, 2 scaffold, 8 strand, 11 env_from, 12 env_to
COL_OG, COL_SCAFFOLD, COL_STRAND, COL_ENV_FROM, COL_ENV_TO = 0, 2, 8, 11, 12
# Take species from the filename
speciesname = os.path.basename(results.input).split("_nuwt_hits")[0]

coords = {}
m = open(results.input, "r")
for record in m:
    if record.startswith("#") or not record.strip():
        continue
    f = record.split()
    readname = f[COL_SCAFFOLD]
    OG = f[COL_OG]
    start = f[COL_ENV_FROM]
    stop = f[COL_ENV_TO]
    if readname not in coords:
        coords[readname] = {}
        i = 1
    coords[readname][i] = {}
    coords[readname][i]["coords"] = (start, stop)
    coords[readname][i]["name"] = OG
    i = i + 1
m.close()

# Extract the coordinates and write one fasta for each orthogroup
tab = open(results.table, "a")

for seq_record in SeqIO.parse(results.fasta, "fasta"):
    if seq_record.id not in coords:
        continue
    for elem in coords[seq_record.id]:
        a = int(coords[seq_record.id][elem]["coords"][0])
        b = int(coords[seq_record.id][elem]["coords"][1])
        OG = coords[seq_record.id][elem]["name"]

        if a > b:
            start, stop, strand = b, a, "-"
            seq = seq_record.seq[start - 1:stop].reverse_complement()
        else:
            start, stop, strand = a, b, "+"
            seq = seq_record.seq[start - 1:stop]                      

        header = "INS_" + speciesname + "_" + str(seq_record.id) + "_" \
                 + OG + ":" + str(start) + "-" + str(stop)              

        with open(os.path.join(results.outdir, OG + ".fa"), "a") as f:
            f.write(">" + header + "\n" + str(seq).upper() + "\n")

        tab.write("\t".join([header, speciesname, str(seq_record.id), OG,
                             strand, str(start), str(stop),
                             str(stop - start + 1)]) + "\n")

tab.close()