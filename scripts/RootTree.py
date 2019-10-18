import Bio.Phylo
import argparse
import sys
sys.setrecursionlimit(2000)

parser = argparse.ArgumentParser("Root tree based on outgroup.")
parser.add_argument("input",help="Input tree in newick format.")
parser.add_argument("outgroup",help="Name of outgroup to use for rooting.")
parser.add_argument("output",help="Output tree in newick format.")

args=parser.parse_args()

# Import tree.
tree=Bio.Phylo.read(args.input,"newick")

# Identify outgroup and root tree.
for clade in tree.find_clades():
	if clade.name==args.outgroup:
		tree.root_with_outgroup(clade)

# Export tree.
Bio.Phylo.write(tree,args.output,"newick")