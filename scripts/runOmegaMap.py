import Bio
from Bio import SeqIO,AlignIO
from Bio.Align import MultipleSeqAlignment
import os
import subprocess
import sys

wdir = "~/Work/jchar_dn_ds/denovo_dn_ds/"

def strip_stop_codons(f):
  '''opens alignment file and strips terminal stop codons from each sequence. writes stripped alignments to original file'''
  if os.stat(f).st_size > 0: #check file size in bytes
    aln = AlignIO.read(f, 'fasta')
    for i in aln:
      if str(i.seq).endswith(('TAA','TAG','TGA')):
        i.seq = i.seq[:-3]
    AlignIO.write(f, aln, 'fasta')

def write_ini(f):
  '''writes .ini file for omegaMapML'''
  orders = subprocess.check_output(["/usr/local/bin/omegaMap-mle-dist/order", str(len(aln)), "1"])
  ini_file = "omml.ini"
  ini_out = open(inifile, 'w')
  ini_template= '''
    FASTA = {a}
    pi = .016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393
    norders = 1
    orders = {b}
    outfile = {c}omega_map_output/{a}.out
    '''
  ini=ini_template.format(a = file,  b = orders, c= wdir)
  ini_out.write(ini)
  ini_out.close

f = "glnA+2.mafft.fa"
strip_stop_codons(f)
write_ini(sys.argv[1])
#run omegaMapML
subprocess.call(["/usr/local/bin/omegaMap-mle-dist/omegaMapML", "omml.ini"])




