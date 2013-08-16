import Bio
import sys
from Bio import SeqIO, AlignIO, SeqFeature
from sys import argv

def get_annotation(ref):
        '''read in a Genbank annotation file'''
        ref_genome = SeqIO.read(ref, "genbank")
        return ref_genome

if __name__ == '__main__':

#get command line params
        args = sys.argv[1:]
        try:
                anno = str(args[0])     #path to ref genome

        except ValueError:
                fail('Expected string arguments got %d' % int(args))

#get reference genome
        ref_genome = get_annotation(anno)

#open file for writing gene list to
	outfile = open("MRSA_252_list.txt", 'w')

#get each gene using Biopython's parser
        for f in ref_genome.features:
                if f.type == "CDS":
			loc = str(f.location).rsplit('(')
			ref_range = loc[0].strip('[]').split(':')
			strand = loc[1].rstrip(')')
			if ref_range[0].isdigit() and ref_range[1].isdigit():
				if strand == "+":
					start = str(int(ref_range[0]) + 1)
					end = str(int(ref_range[1]) -1)
				elif strand == "-":
					start = str(int(ref_range[0]) +1)
					end = str(int(ref_range[1])- 1) 
				strand = loc[1].rstrip(')')
				id = str(f.qualifiers['locus_tag']).strip("[]'")
				outline = "%s\t%s\t%s-%s" %(id, strand, start, end)
				outfile.write(outline)
				outfile.write("\n")
			
