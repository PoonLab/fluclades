"""
initialize list with original tree
for each subtree in list
    locate longest internal branch
    if branch length exceeds threshold
        cut the branch
        replace the subtree with the two new subtrees
    continue until no internal branches are longer than threshold
"""

from Bio import Phylo
from Bio.Phylo.BaseTree import Tree
import argparse
import bisect
import sys
from io import StringIO


parser = argparse.ArgumentParser(
    description="Partition tree by cutting on internal branches with length "
    "exceeding threshold."
)
parser.add_argument("infile", type=argparse.FileType('r'), 
                    help="<input> Path to Newick file to process.")
parser.add_argument("outfile", type=argparse.FileType('w'), nargs="?", default=sys.stdout,
                    help="<output, optional> File to write output, defaults to stdout.")
parser.add_argument("--cutoff", type=float,
                    help="<input> Maximum internal branch length. "
                    "If none given, display summary of lengths.")
parser.add_argument("--nbin", type=int, default=20, 
                    help="<option> number of bins for branch length summary.")
parser.add_argument("-f", "--format", choices=["summary", "labels", "trees"], default="summary",
                    help="<option> Format to write output. Defaults to 'summary'. "
                    "'labels' writes all tip labels for each subtree index. "
                    "'trees' writes a set of Newick tree strings.")
args = parser.parse_args()

phy = Phylo.read(args.infile, 'newick')

if args.cutoff is None:
    bl = [node.branch_length for node in phy.get_nonterminals() if node is not phy.root]
    bl.sort()  # ascending order
    blmax = bl[-1]
    blstep = blmax/args.nbin
    cutval = 0
    left = 0
    for i in range(args.nbin):
        cutval += blstep
        right = bisect.bisect_left(bl, cutval)
        print(cutval, right-left)
        left = right
    sys.exit()


def get_parents(phy):
    # Clades do not store references to parents
    parents = {}
    for parent in phy.get_nonterminals():
        for child in parent.clades:
            parents.update({child: parent})
    return parents


def unroot(phy):
    """ 
    Shorthand, designates a trifurcating Clade as root.
    Note this sets phy.rooted=True
    """
    if len(phy.root.clades) == 2:
        phy.root_with_outgroup(phy.root.clades[0])


def cuttree(phy, clade):
    """
    :param phy: BaseTree.Tree object
    :param clade: BaseTree Clade object, corresponds to internal branch to cut
    :returns: new BaseTree.Tree object rooted on clade; input `phy` is modified 
              in-place with clade removed.
    """
    if clade.is_terminal():
        sys.stderr.write("Warning: cannot cut tree on terminal branch. Use prune method. "
                         "Returning unaltered tree.\n")
        return phy
    
    parents = get_parents(phy)
    parent = parents.get(clade, None)
    if parent is None:
        sys.stderr.write("Clade has no parent in tree\n")
        return phy
    
    clade.branch_length = 0
    subtree1 = Tree(clade)
    unroot(subtree1)

    parent.clades.remove(clade)
    nchild = len(parent.clades)
    
    grandpar = parents.get(parent, None)
    if grandpar is None:
        # parent is current root    
        if nchild == 0:
            sys.stderr.write("cuttree: parent had only one child\n")
            sys.exit()
        elif nchild == 1:
            # parent is stem of one subtree, drop it
            subtree2 = Tree(parent.clades[0])
            del parent
        else:
            subtree2 = Tree(parent)
    else:
        if nchild < 2:
            for child in parent.clades:
                child.branch_length += parent.branch_length
                grandpar.clades.append(child)        
            grandpar.clades.remove(parent)
            del parent
        phy.root_with_outgroup(grandpar)
        subtree2 = Tree(phy.root)
    
    subtree2.root.branch_length = 0
    unroot(subtree2)
    return subtree1, subtree2


def longest(phy):
    """ Get longest internal branch """
    nodes = phy.get_nonterminals()
    if len(nodes) < 2:
        return None  # only terminal descendants left
    intermed = [(node.branch_length, ni) for ni, node in enumerate(nodes)
                if node.branch_length]
    intermed.sort(reverse=True)
    return nodes[intermed[0][1]]

"""
s = StringIO("(((t1:0.06559279552,t2:0.06559279552):0.7043641405,(t6:0.5937900807,"
             "((t4:0.0129620425,t10:0.0129620425):0.005046938947,t8:0.01800898145):"
             "0.5757810992):0.1761668554):1.108476543,((t5:0.1602943302,t9:0.1602943302)"
             ":0.2262072535,(t7:0.09450417782,t3:0.09450417782):0.2919974059):1.491931895);")
phy = Phylo.read(s, "newick")
"""

unroot(phy)
subtrees = [phy]
cutting = True  # enter loop
while cutting:
    cutting = False
    newtrees = []
    for subtree in subtrees:
        node = longest(subtree)
        if node and node.branch_length > args.cutoff:
            cutting = True
            st1, st2 = cuttree(subtree, node)
            newtrees.extend([st1, st2])
        else:
            newtrees.append(subtree)
    subtrees = newtrees

# write output
if args.format == "labels":
    args.outfile.write("subtree,tip.label\n")
    for idx, subtree in enumerate(subtrees):
        for tip in subtree.get_terminals():
            args.outfile.write(f"{idx},{tip.name}\n")
elif args.format == "summary":
    args.outfile.write("subtree,ntips,tip1\n")
    for idx, subtree in enumerate(subtrees):
        tips = subtree.get_terminals()
        args.outfile.write(f"{idx},{len(tips)},{tips[0].name}\n")
else:
    for subtree in subtrees:
        Phylo.write(subtree, args.outfile, "newick")

if args.outfile is not sys.stdout:
    args.outfile.close()
