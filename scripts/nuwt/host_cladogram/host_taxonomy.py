#!/usr/bin/env python3
"""
Build a host cladogram with NCBI taxonomy tree using ETE3 NCBITaxa and a local taxonomy datavase at ~/.etetoolkit/taxa.sqlit

Input: metadata.tsv with headers species, taxid, suborder, family 

Usage:
/data/projects/p2025-0083_mining_cobionts/.conda_envs/ete3/bin/python host_taxonomy metadata_tsv out_prefix
"""

import sys
import csv
from ete3 import NCBITaxa

def main():
    meta_tsv, out_prefix = sys.argv[1], sys.argv[2]
    ncbi = NCBITaxa()
    rows = []
    # Read the file into a dictionary
    with open(meta_tsv, newline="") as fh:
        for species_entry in csv.DictReader(fh, delimiter="\t"):
            rows.append({
                "species": species_entry["species"],
                "taxid": int(species_entry["taxid"]),
                "suborder": species_entry["suborder"],
                "family": species_entry["family"],
            })
    # Build the host tree: get_topology returns the smallest tree that connects the taxids
    tree = ncbi.get_topology([species_entry["taxid"] for r in rows])
    # Rename the tips from te taxid to the species name
    taxaid_sp = {species_entry["taxid"]: species_entry["species"] for species_entry in rows}
    for leaf in tree.get_leaves():
        leaf.name = taxaid_sp[int(leaf.name)]

    # Newick format: format=9 to have the leaf names only
    nwk_path = out_prefix + ".nwk"
    tree.write(outfile=nwk_path, format=9)
    # Get the species that are tips in the tree
    tip_species = {l.name for l in tree.get_leaves()}
    tips_path = out_prefix + "_tips.tsv"
    with open(tips_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["species", "taxid", "suborder", "family"])
        for species_entry in rows:
            if species_entry["species"] in tip_species:
                w.writerow([species_entry["species"], species_entry["taxid"], species_entry["suborder"], species_entry["family"]])

if __name__ == "__main__":
    main()