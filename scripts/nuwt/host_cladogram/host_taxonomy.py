#!/usr/bin/env python3
"""
Build a host cladogram
This is an NCBI taxonomy tree, not a molecular phylogeny

Here we use ete3's NCBITaxa, whih uses a local taxonomy databse built at ~/.etetoolkit/taxa.sqlit
Input files:
metadata.tsv with header that contains the host genome metadata with species, taxid and suborder, family  kept for colouring

Output files:
.nkw: the cladogram, with tip labels as species
_tips.ts: one row per tip that actually got in the tree for the ggtree colouring 

Usage: 
/data/projects/p2025-0083_mining_cobionts/.conda_envs/ete3/bin/python build_host_tree.py <metadata_tsv> <out_prefix>
"""

import sys
import csv
from ete3 import NCBITaxa

def main():
    meta_tsv, out_prefix = sys.argv[1], sys.argv[2]
    ncbi = NCBITaxa()
    rows = []
    with open(meta_tsv, newline="") as fh:# COlumns are fixed in the metadata sheet 
        for r in csv.DictReader(fh, delimiter="\t"):
            rows.append({
                "species": r["species"],
                "taxid": int(r["taxid"]),
                "suborder": r["suborder"],
                "family": r["family"],
            })
    # Build the tree
    # get_topology returns the smallest tree connecting these taxids
    tree = ncbi.get_topology([r["taxid"] for r in rows])

    # Rename the tips from taxid to species name
    # Map taxid and species
    tx2sp = {r["taxid"]: r["species"] for r in rows}
    for leaf in tree.get_leaves():
        leaf.name = tx2sp[int(leaf.name)]

    # Write to Newick format:
    # format=9 writes leaf names only
    nwk_path = out_prefix + ".nwk"
    tree.write(outfile=nwk_path, format=9)

    # Tip annotation TSV for ggtree: only species that ended up as tips in the tree.
    tip_species = {l.name for l in tree.get_leaves()}
    tips_path = out_prefix + "_tips.tsv"
    with open(tips_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["species", "taxid", "suborder", "family"])
        for r in rows:
            if r["species"] in tip_species:
                w.writerow([r["species"], r["taxid"], r["suborder"], r["family"]])

if __name__ == "__main__":
    main()