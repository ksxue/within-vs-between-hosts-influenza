import Bio.Phylo
import argparse

parser = argparse.ArgumentParser("Parse tree with ancestral reconstruction, get mutation frequencies.")
parser.add_argument("tree",help="Input tree in nexus format, as produced by treetime ancestral.")
parser.add_argument("sequences",help="Ancestral sequences in FASTA format, as produced by treetime ancestral.")
parser.add_argument("outgroup",help="Name of outgroup sequence to use for rooting.")
parser.add_argument("output",help="Output table of mutations and their frequencies.")

args=parser.parse_args()

import sys, subprocess, glob, os, shutil, re, importlib, csv, json
from subprocess import call
import collections
from collections import Counter
from Bio import SeqIO
from Bio import Seq
import Bio.Phylo
import re
import math

sys.setrecursionlimit(3000)

# Given a comment on the annotated nexus tree produced by treetime ancestral,
# parse the comment and return a flattened list of mutations.
def comment_to_mutation(comment):
	mutations=comment.split("\"")[1]
	mutations=mutations.split(",")
	return(mutations)
	
# Given a mutation in the format produced by treetime ancestral,
# parse the mutation and return the ancestral identity, site, and derived identity.
def parse_mutation(mutation):
	ancestral=mutation[0]
	site=re.findall("\d+",mutation)[0]
	derived=mutation[-1]
	return(ancestral,site,derived)
	
# Given a DNA triplet, translate it to a one-letter amino-acid code.
def translate(triplet):
	table = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W', 
    }
	AA='X'
	if len(triplet)==3 and triplet in table:
		AA=table[triplet]
	return AA

# Given a mutation in the form ancestral nt, site (1-indexed nt site), and derived nt,
# the name of the sequence that the mutation creates,
# and a list of sequences annotated by name,
# return the 1-indexed codon site, ancestral codon, derived codon,
# ancestral amino acid, and derived amino acid.
def translate_mutation(ancestral, site, derived, seqname, seqs):
	# Determine the codon site at which the mutation occurs.
	# This gives the 1-indexed codon site.
	codonsite=math.ceil(int(site)/3)
	# Verify that the derived nt matches the site in the given sequence.
	if derived!=str(seqs[seqname].seq[site-1]):
		print("Error: mutation does not result in derived sequence.")
	# Extract the derived codon.
	derivedcodon=str(seqs[seqname].seq[math.floor((site-1)/3)*3:math.floor((site-1)/3)*3+3])
	# Determine the ancestral codon.
	ancestralcodon=list(derivedcodon)
	ancestralcodon[(site%3-1)] = ancestral
	ancestralcodon=''.join(ancestralcodon)
	return(codonsite, ancestralcodon, derivedcodon, translate(ancestralcodon), translate(derivedcodon))

# Print extra output when DEBUG is set to True.
DEBUG=False
NODENAME="NODE_0000298"
	
# Import tree.
tree=Bio.Phylo.read(args.tree,"nexus")

# Identify outgroup and root tree.
for clade in tree.find_clades():
	if clade.name==args.outgroup:
		tree.root_with_outgroup(clade)
		break

# Calculate the total number of sequences in the tree.
# Omit the outgroup from this tally.
numtips=tree.count_terminals()-1
		
# Import the ancestral sequences as a SeqIO object.
ancestralseqs=SeqIO.parse(args.sequences,"fasta")
# Convert the ancestral sequences to a dictionary format.
ancestralseqs=SeqIO.to_dict(ancestralseqs)

	
# Iterate through all of the clades in a tree.
# If mutations are present, then print the mutations present.
# Calculate the number of descendents of that clade that carry
# the mutation that occurred to form that clade.
# Subtract descendents in sub-clades that have had a reversion
# or second mutation at the same site.
file=open(args.output,"w")
# Print file header.
file.write("Seq NtAnc NtSite NtDer CodonSite CodonAnc CodonDer AAAnc AADer NumDescendants TotalSeqs\n")
for clade in tree.find_clades(order="level"):
	# Do not analyze mutations that are present on the outgroup sequence.
	if clade.name==args.outgroup:
		continue
	# Determine the name of the clade,
	# which is listed under the "name" attribute for tips
	# and the "confidence" attribute for nodes.
	name=clade.name if str(clade.name)!="None" else clade.confidence
	# Insert a DEBUG setting that focuses on a particular node.
	if DEBUG:
		if name!=NODENAME:
			continue
	# Identify clades that are annotated with mutations.
	if str(clade.comment)!="None":
		# Parse the mutations that are encountered.
		mutations=comment_to_mutation(clade.comment)
		# For each mutation, determine the number of descendents.
		# that carry that mutation.
		# Calculate the total number of descendents in that clade.
		# Then subtract the number of descendents that are present in clades
		# that have had a second mutation at the same site (could be reversion or not).
		for mutation in mutations:
		
			# Determine at which site the mutation occurred.
			mutation=parse_mutation(mutation)
			ancestral=mutation[0]
			site=int(mutation[1])
			derived=mutation[2]
			# Determine the one-indexed codon site at which the mutation occurred.
			codonsite=math.ceil(int(site)/3)
			# Calculate the total number of descendents of this node.
			alldescendents=clade.count_terminals()
			if DEBUG:
				print(name, mutation, alldescendents)
			
			# Identify sub-clades that have second-hit mutations.
			secondhitclades=[]
			for subclade in clade.find_clades(order="level"):
				# Determine the name of the clade,
				# which is listed under the "name" attribute for tips
				# and the "confidence" attribute for nodes.
				subclade_name=subclade.name if str(subclade.name)!="None" else subclade.confidence
				# Parse mutations that occur in the subclade.
				# Exclude clades that are the same as the parent clade
				# or that are descendants of clades that already have a
				# second-hit mutation.
				if str(subclade.comment)!="None" and subclade!=clade and subclade not in secondhitclades:
					# Ensure that the subclade is not descended from a subclade
					# that is already known to have a second-hit mutation.
					descendant=False
					for secondhitclade in secondhitclades:
						if secondhitclade.is_parent_of(subclade):
							descendant=True
							if DEBUG:
								print("third hit clade")
					if not descendant:
						# Parse mutations associated with subclade.
						subclade_mutations=comment_to_mutation(subclade.comment)
						# Iterate through all of the mutations in the subclade
						# and determine if any mutations occur at the same site
						# as the mutation in the original clade.
						for subclade_mutation in subclade_mutations:
							# Determine at which site the mutation occurred in the subclade.
							subclade_site=parse_mutation(subclade_mutation)[1]
							# Determine the codon site at which the mutation occurred in the subclade.
							subclade_codonsite=math.ceil(int(subclade_site)/3)
							# If the subclade has a mutation at the same codon site as
							# the original mutation, then subtract the number of descendents
							# in the subclade.
							# Add the subclade to the list of second-hit clades.
							if subclade_codonsite==codonsite:
								alldescendents=alldescendents-subclade.count_terminals()
								secondhitclades.append(subclade)
								if DEBUG:
									print("second hit clade", subclade_name, subclade_mutation, subclade.count_terminals())
								# If a mutation is found that constitutes a second hit,
								# exit the loop and stop analyzing the remaining mutations.
								# This break statement accounts for the uncommon case in which
								# a clade has two mutations at the same codon site.
								break
			# Extract the attributes of the mutation for output.
			codon=translate_mutation(ancestral,site,derived,name,ancestralseqs)
			file.write("%s %s %d %s %d %s %s %s %s %d %d\n" % (name, ancestral, site, derived, 
			codon[0], codon[1], codon[2], codon[3], codon[4], alldescendents, numtips))
file.close()

