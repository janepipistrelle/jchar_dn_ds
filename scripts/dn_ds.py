# coding: utf-8
seq1
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment
from cogent import LoadSeqs,DNA
ref = SeqIO.read("sequence.gb", "genbank")
aln = AlignIO.read("all_carriage_isolates.fasta, "fasta")
aln = AlignIO.read("all_carriage_isolates.fasta", "fasta")
ids = [line.strip() for line in open("scripts/carriage_isolate_ids.txt")]
filtered_aln = MultipleSeqAlignment([])
get_ipython().magic(u'cpaste')
get_ipython().magic(u'cpaste')
get_ipython().magic(u'cpaste')
get_ipython().magic(u'cpaste')
for f in ref.features:
    if f.type == "CDS":
        gene = clean_seqs(slice_seqs(f, filtered_aln))
        
print gene[0]
seq1 = DNA.makeSequence(str(gene[0].seq)
)
seq1
codons1 = seq1.getInMotifSize(3)
codons1
codons1[0]
seq2 = DNA.makeSequence(str(gene[1].seq)
)
codons1[]
codons1 = seq2.getInMotifSize(3)
codons1[0]==codons2[0]
codons1 = seq1.getInMotifSize(3)
codons2 = seq2.getInMotifSize(3)
codons1[0]==codons2[0]
codons1[1]==codons2[1]
codons1[1]==codons2[1]
seq1
seq_1 = str(gene[0].seq)
seq_1 = str(gene[1].seq)
seq_1 = str(gene[0].seq)
seq_2= str(gene[1].seq)
seq_1
codon_1 = seq_1[0:2]
codon_1
codon_1 = seq_1[0:3]
codon_1
seq_1[0]
seq_1[1]
seq_1[2]
codon_1 = seq_1[4:6]
codon_1
codon_1 = seq_1[3:6]
codon_1
codon_1 = seq_1[6:9]
codon_1
for c in range(10,3):
    print c
    
for c in range(10):
    print c
    
for c in range(10,3):
    print c
    
for c in range(0,10,3):
    print c
    
for c in range(0, len(str_1), 3):
    codons.append(str_1[c:c+3]
    
    
    )
    
for c in range(0, len(seq_1), 3):
             codons.append(seq_1[c:c+3])
    
codons1=[]
for c in range(0, len(seq_1), 3):
             codons1.append(seq_1[c:c+3])
    
codons1
from Bio.Data import CodonTable
print CodonTable
standard_table = CodonTable.unambiguous_dna_by_id[1]
dir(CodonTable)
CodonTable.list_ambiguous_codons
CodonTable.NCBICodonTable
print CodonTable.NCBICodonTable
print standard_table
print standard_table[1]
print standard_table.codons
print standard_table.NCBIcodons
dir(CodonTable)
dir(CodonTable\])
dir(CodonTable)
print dir(NCBICodonTable)
print dir(CodonTable.NCBICodonTable)
type(standard_table)
standard_table = CodonTable
print standard_table
standard_table
standard_table = CodonTable.NCBICodonTable
print standard_table
standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
print standard_table
standard_table.protein_alphabet
standard_table.back_table
standard_table.forward_table
genetic_code= CodonTable.unambiguous_dna_by_name["Standard"].forward_table
genetic_code
genetic_code["TTT"]
genetic_code["AAA"]
genetic_code[codons1[0]]== genetic_code[codons2[0]]
for i in codons1:
    if codons1[i] != codons2[i]:
        if genetic_code[codons1[i]] == genetic_code[codons2[i]]:
            ds += 1
            