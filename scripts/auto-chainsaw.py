from chainsaw import chainsaw, mutual_info, unroot
import sys
import csv
from Bio import Phylo
import argparse

description = """
This script is used to generate the data required to produce Figures
2A and 3A.  The inputs were trees reconstructed using FastTree2.
Results are written to stdout in CSV format.
"""

parser = argparse.ArgumentParser(description)
parser.add_argument("infile", type=str, help="Input tree")
parser.add_argument("gene", choices=["HA", "NA"], help="Which gene/segment?")
args = parser.parse_args()

if args.gene == "HA":
    cutoffs = [0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15,
               0.16, 0.17, 0.175, 0.1775, 0.18, 0.1825, 0.185, 0.19, 0.2, 0.25,
               0.26, 0.265, 0.2675, 0.27, 0.275, 0.28, 0.29, 0.3]
else:
    cutoffs = [0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.2, 0.21, 0.215,
               0.22, 0.222, 0.225, 0.23, 0.24, 0.25, 0.3, 0.4, 0.405, 0.41, 0.425, 0.43,
               0.4325, 0.435, 0.44, 0.45, 0.46, 0.48, 0.49, 0.5]

writer = csv.writer(sys.stdout)
writer.writerow(["cutoff", "nsubtrees", "mutual.inf", "normalized"])
for cutoff in cutoffs:
    phy = Phylo.read(args.infile, "newick")
    unroot(phy)
    subtrees = chainsaw(phy, cutoff)
    minfo, norm_mi = mutual_info(subtrees, hema=(args.gene == 'HA'))
    writer.writerow([cutoff, len(subtrees), minfo, norm_mi])
    sys.stdout.flush()
