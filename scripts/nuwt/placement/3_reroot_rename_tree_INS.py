#!/usr/bin/env python3

# Adapted from Vancaester et al. Reroot_rename_tree.py.
# Edits to make it work iwith our setup:
# 1. NUWTs are now recognized by INS_ header prefix
# 2. Outgroup detection is now case insensitive
# 3. The conversion table matches our 2 columns one with accession tab supergroup
# 4. We keep in stdout the results row only with OG tab INS tab supergroup tab nuwt_name
# 
# Labels each reference leaf with its supergroup, collapse same-supergroup clades, root on the outgroup MRCA (or midpoint if no outgroup is present in this OG), then call each NUWT by the supergroup of the smallest clade whose reference leaves are all one supergroup.
#
# Usage:
# python Reroot_rename_tree_INS.py -t <treefile> -a <aln.fasta> -og <OG_ID> -c <conversion.tsv> -m <manual_skip.txt>

import sys
import argparse
from ete3 import Tree

parser = argparse.ArgumentParser()
parser.add_argument("-t",  dest="tree", help="IQ-TREE .treefile")
parser.add_argument("-a",  dest="aln",  help="alignment FASTA (headers are read)")
parser.add_argument("-og", dest="og",   help="orthogroup id")
parser.add_argument("-c",  dest="conv", help="conversion table: accession<TAB>supergroup")
parser.add_argument("-m",  dest="man",  help="OGs to skip (one per line); may be empty")
results = parser.parse_args()

# Outgroups compared case insensitive
OUTGROUP = {x.lower() for x in ["Ace", "Ama", "Ech", "Eru", "mAli", "mBlo", "mPat"]}


def collapsed_leaf(node):
    return len(node2labels[node]) == 1


def node_name_internal(node):
    if len(node2labels[node]) == 1:
        return True, node2labels[node]
    return False, node.name


# Conversion table: 2 columns
species2supergroup = {}
with open(results.conv) as fh:
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2:
            species2supergroup[p[0]] = p[1]

# Manual skip list
manual = []
with open(results.man) as fh:
    for line in fh:
        manual.append(line.strip())

# Read alignment headers: outgroup or reference or NUWT 
conversiontable = {} # full header to supergroup for references only
genome_hits = [] # NUWTs
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

# Label leaves, collapse uniform clades
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

# Root: outgroup MRCA if present, else midpoint
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

# Assign each NUWT by its uniform reference sister clade
hitinfo = {}
for hit in genome_hits:
    if not ancestor.search_nodes(name=hit):
        print(hit, file=sys.stderr)
        for node in t3.traverse("postorder"):
            if node.name == "":
                leafnames = []
                leafnames_all = []
                for leaf in node.iter_descendants():
                    if leaf.name != "":
                        leafnames_all.append(leaf.name)
                        if not leaf.name.startswith("INS_"):  # changed here to recognize INS_
                            leafnames.append(leaf.name)
                if hit in leafnames_all:
                    print(node, file=sys.stderr)
                    if len(set(leafnames)) == 1:
                        hitinfo[hit] = leafnames[0]
                        break

# Results row: 
for hit in hitinfo:
    if results.og not in manual:
        sys.stdout.write("%s\t%s\t%s\t%s\n" %
                         (results.og, hit.split("_")[0], hitinfo[hit], hit))