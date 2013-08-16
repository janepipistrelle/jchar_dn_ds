'''workflow to extract protein-coding sequences from a set of mapcalled genomes. This script exploits the fact that mapped reads share coordinates with the reference genome. It outputs only those sequences shared between a given isolate and a reference genome to which all isolates have been mapped.
'''
import Bio
import sys
import gzip
from Bio import SeqIO, AlignIO, SeqFeature
from Bio.Align import MultipleSeqAlignment
from Bio.Data import CodonTable
from sys import argv

#load genetic code table
genetic_code= CodonTable.unambiguous_dna_by_name["Standard"].forward_table

def fail(msg):
	'''Print error message and exit'''
	print >> sys.stderr, msg
	sys.exit(1)

def get_annotation(ref):
	'''read in a Genbank annotation file'''	
	ref_genome = SeqIO.read(ref, "genbank")
	return ref_genome


def load_seqs(data_dir, ref_code, comid_list,):
	''' load a given list of comids from the gorm v3 directory structure- need to add testing for duplicates/missing data'''
	aln = MultipleSeqAlignment([])
	ids = [line.strip() for line in open(comid_list)]
	for i in ids:
		search_path = i.lstrip("C,0")
		if len(search_path) == 3:
			folder = "00000000"
		elif len(search_path) == 4:
			folder = "000%s0000" % search_path[0]
		elif len(search_path) == 5:
			folder = "000%s000" % search_path[0:2]
		path = "%s/%s/%s_%s.mapcall.fasta.gz" %(data_dir, folder, i, ref_code)
		try:
			seq_handle = gzip.open(path, 'r')
			seq = SeqIO.read(seq_handle, "fasta")
			aln.append(seq)
		except IOError:
   			print 'File does not exist: %s' % i
	return aln


def slice_seqs(f, aln):
	'''For a given CDS in the annotation, get that slice of the alignment'''
	gene = []
	for i in aln:
		gene_record = (f.extract(i))
		gene_record.id = i.id
		gene.append(gene_record)
	return gene

def clean_seqs(gene):
	'''clean up sequences to remove N & - characters'''
	clean_gene = MultipleSeqAlignment([])
	for genome in gene:
		if genome.seq.count("N") + genome.seq.count("-") <  0.1*(len(genome.seq)):
			clean_gene.append(genome)
	return clean_gene

			
def pairwise_counts(sample1, sample2):
	'''tally numbers of synonymous and nonsynonymous differences for a pair of sequences'''
	Dn = 0
	Ds = 0
	conserved = 0
	ambiguous = 0
	stop = 0

	codons1= get_codons(sample1)
	codons2= get_codons(sample2)

	for i in range(0,len(codons1)):
			if ambiguous_codon(codons1[i], codons2[i]) == 1:
				ambiguous += 1
			elif stop_codon(codons1[i], codons2[i]) == 1:
					stop += 1
			elif codons1[i] == codons2[i]:
				conserved += 1
			else: 
				codon_score = codon_compare(codons1[i], codons2[i])
				if codon_score == 0:
					Ds += 1
				elif codon_score == 1:
					Dn +=1
	return [Dn, Ds, conserved, stop, ambiguous]

def get_codons(sequence):
        '''split a sequence into codon. check that the specified sequence is a CDS.'''
        stop_codons = set(['TAA', 'TAG', 'TGA'])
        codons = []
        seq = str(sequence.seq)
	assert len(seq) % 3 == 0, "Sequence is not a CDS"
        for c in range(0, len(seq), 3):
                codons.append(seq[c:c+3])
        del codons[-1] #trim terminal stop codon
	return codons   

def ambiguous_codon(codon1, codon2):
	'''checks for ambiguous codons'''
	chars = ['N', '-']
	for c in chars:
		if codon1.find(c) != -1 or codon2.find(c) != -1:
			return 1

def stop_codon(codon1, codon2):
	'''checks if a codon is a stop codon'''
	stop_codons=['TAA','TAG','TGA']
	for c in stop_codons:
		if codon1.find(c) != -1 or codon2.find(c) != -1:
			return 1
	 
def codon_compare(codon_1, codon_2):
	'''compare codons from two sequences. return 0 for synonymous and 1 for nonsyonymous sub '''
	if genetic_code[codon_1] == genetic_code[codon_2]:
		return 0
	else:
		return 1

def p_prior(Dn_total, Ds_total):
	'''calculate p-prior parameter from counts summed across all pairwise comparisons'''
	p_prior = round(Dn_total/float(Dn_total + Ds_total),4)
	return p_prior

def beta_params(p_prior, n_prior):
	'''calculate a and b parameters for beta distribution'''
	a = n_prior*p_prior
	b = n_prior*(1-p_prior)
	return [a,b]

def dn_ds_bayes(a,b, gene_counts):
	'''calculates posterior estimate of dn/ds, adjusted for to account for prior; outputs dictionary of values for all pairwise comparisons'''
	dn_ds_sum = 0
	for i in gene_counts:
		Dn=float(i[0])
		Ds=float(i[1]) 
		dn_ds_sum += round(((Dn + a)/(Ds + b))/3, 4)
	avg_dn_ds = dn_ds_sum/len(gene_counts)
	return avg_dn_ds
		
#Main

if __name__ == '__main__':

#get command line params
	args = sys.argv[1:]
	try:
		data_dir = str(args[0]) #path to data directory
		ref_code = str(args[1]) #code for reference genome
		seq_list = str(args[2]) #path to list of comids
		anno = str(args[3])	#path to ref genome
		nprior =float(args[4]) #set weighting for beta prior

	except IndexError:
		fail('Expected 4 arguments, got %d' % len(args))
	except ValueError:
		fail('Expected string arguments got %d' % int(args))

#get reference genome
	ref_genome = get_annotation(anno)

#load sequences
	aln = load_seqs(data_dir, ref_code, seq_list)

#open file for writing dn/ds estimates to
	outstring = "%s_dn_ds" % seq_list
	dn_ds_genome = open(outstring, 'w')

#get alignments and id string for each CDS in genome
	for f in ref_genome.features:
		if f.type == "CDS":
			gene = clean_seqs(slice_seqs(f, aln))
			gene_id= f.qualifiers['locus_tag'][0]
#step through pairs of genomes in sample; store counts a file
			counts_file = open(gene_id, 'w')
			counts_file.write("Dn\tDs\tconserved\tstop\tambiguous")
	
			gene_counts = []
			Dn_total = 0
			Ds_total = 0
	
			if len(gene) >= 2:
				for i in range(0,len(gene)):		
					for j in range(i, len(gene)):
                                        	if i != j:
							counts = pairwise_counts(gene[i],gene[j])
							counts_file.write('\t'.join(map(str,counts)))
							counts_file.write('\n')
							Dn_total += counts[0]
                                               		Ds_total += counts[1]
                                                	gene_counts.append(counts)
		
#calculate Bayesian estimate of dn/ds; averaged across all pairs
				if Dn_total == 0 or Ds_total == 0:
					a,b = 3,1
				else:
					a,b = beta_params(p_prior(Dn_total, Ds_total), nprior)
				gene_dn_ds = dn_ds_bayes(a, b, gene_counts)			
				dn_ds_genome.write('\t'.join(map(str,[gene_id, gene_dn_ds, a,b, "\n"])))
				print('\t'.join(map(str,[gene_id, gene_dn_ds, a,b])))
