#workflow to extract protein-coding sequences from a set of mapcalled genomes. This script exploits the fact that mapped reads share coordinates with the reference genome. It outputs only those sequences shared between a given isolate and a reference genome to which all isolates have been mapped.

from Bio import SeqIO, AlignIO
from Align import MultipleSeqAlignment
from sys import argv

seq_file, seq_list, anno = argv

def get_annotation(anno):
	ref = SeqIO.read(anno, "genbank")
	return ref

#load sequences and return alignment filtered against list of ids
def load_seqs(seqs, list):
	aln = AlignIO.read(seqs, "fasta")
	ids = [line.strip() for line in open(list)]
	filtered_aln = MultipleSeqAlignment([])
	for i in aln:
		for j in ids:
			if i.id.find(j) != -1:
				filtered_aln.append(i)
	return filtered_aln
		

#slice sequences to get "alignments" for each gene
def slice_seqs(f, aln)
		gene = []
		for i in aln:
			gene.append(f.extract(i))
		return gene	
		

#step through alignment for each gene in a pairwise manner
def pairwise_comp

#for each pair of genomes count dN and ds using Nei & Gojobori method as implmented in paml yn00)
def count_diffs()
def count_sites()
def calc_subs()

#actual run
ref_genome = get_annotation(anno)

filtered_aln = load_seqs(seq_file, seq_list)

for f in ref_genome:
	if f.type = "CDS":
		gene = slice_seqs(f, filtered_aln)

