#!/usr/bin/env python3
"""
This script was adapted from Vancaester et al. Reroot_rename_tree.py
In here, the nuwt fragments carry the INS_ header prefix.
Our accession to supergroup table contains two columns: GCA_accession   supergroup
Usage:
python3 Reroot_rename_tree_INS.py -t treefile -a alignment.fasta -og og_id -c accession_to_supergorup.tsv -m skipped_ogs
"""

import sys
import argparse
from ete3 import Tree

parser = argparse.ArgumentParser()
parser.add_argument("-t",  dest="tree", help="IQ-TREE .treefile")
parser.add_argument("-a",  dest="aln",  help="alignment FASTA")
parser.add_argument("-og", dest="og",   help="OG id")
parser.add_argument("-c",  dest="conv", help="accession to supergroup table")
parser.add_argument("-m",  dest="man",  help="OG to skip")
results = parser.parse_args()

# Define outgroups
OUTGROUP = {x.lower() for x in ["Ace", "Ama", "Ech", "Eru", "mAli", "mBlo", "mPat"]}

def collapsed_leaf(node):
    if len(node2labels[node]) == 1:
        return True
    else:
        return False

def node_name_internal(node):
    if len(node2labels[node]) == 1:
        name = node2labels[node]
        return True, name
    else: 
        name = node.name
        return False, name

# Define conversion table: accession tab supergroup
species2supergroup = {}
with open(results.conv) as fh:
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2:
            species2supergroup[p[0]] = p[1]

# Define the OGs to skip
manual = []
with open(results.man) as fh:
    for line in fh:
        manual.append(line.strip())

# Read alignment header to find outgroups or references
# Dictionary of header to supergroup for the references
conversiontable = {}
# List of insect fragments
genome_hits = []
outgroupseqs = []
inseqs = []
with open(results.aln) as fh:
    for line in fh:
        line = line.strip()
        if not line.startswith(">"):
            continue
        header = line[1:]
        if header.split("_")[0].lower() in OUTGROUP:
            outgroupseqs.append(header)
        else:
            if ":" in header:
                header = header.replace(":", "_")
            inseqs.append(header)
            speciesname = "_".join(header.split("_")[:-1])
            if speciesname in species2supergroup:
                conversiontable[header] = species2supergroup[speciesname]
            else:
                genome_hits.append(header)

t = Tree(results.tree)

for leaf in t.iter_leaves():
    if leaf.name in conversiontable:
        leaf.add_features(supergroup=conversiontable[leaf.name])
    else:
        leaf.add_features(supergroup=leaf.name)

node2labels = t.get_cached_content(store_attr="supergroup")

for node in t.traverse("postorder"):
    nm = str(node_name_internal(node)[1])
    if "{" in nm:
        node.name = nm.split("{'")[1].split("'}")[0]
    else:
        node.name = nm

t2 = Tree(t.write(is_leaf_fn=collapsed_leaf))
t3 = t2

# Root to outgroup or midpoint if no outgroup
if len(outgroupseqs) > 1:
    for node in t2.traverse("postorder"):
        t3.set_outgroup(node)
        ancestor = t3.get_common_ancestor(outgroupseqs)
        if ancestor is not t2:
            t3.set_outgroup(ancestor)
            break
elif len(outgroupseqs) == 1:
    ancestor = t3.search_nodes(name=outgroupseqs[0])[0]
    t3.set_outgroup(ancestor)
else:
    ancestor = t3.get_midpoint_outgroup()
    t3.set_outgroup(ancestor)

print(t3, file=sys.stderr)

# Assign nuwt 
hitinfo = {}
# Loop through one nuwt hit at a time
for hit in genome_hits:
    # Skip it is inside outgroups
    if not ancestor.search_nodes(name=hit):
        print(hit, file=sys.stderr)
        for node in t3.traverse("postorder"):
            if node.name == "":
                leafnames = []
                leafnames_all = []
                for leaf in node.iter_descendants():
                    if leaf.name != "":
                        leafnames_all.append(leaf.name)
                        # Changed condition: our genomes start with INS
                        if not leaf.name.startswith("INS_"):
                            leafnames.append(leaf.name)
                if hit in leafnames_all:
                    print(node, file=sys.stderr)
                    if len(set(leafnames)) == 1:
                        hitinfo[hit] = leafnames[0]
                        break

for hit in hitinfo:
    if results.og not in manual:
        print(results.og+'\t'+hit.split('_')[0]+'\t'+hitinfo[hit]+'\t'+hit)