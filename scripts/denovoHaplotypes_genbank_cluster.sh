#version of this code to run using set of 16 Staph reference genomes from Genbank, as found: /dipro/mmm/gorm/v2/ana/Saur/xcarriage/mapping/v2b/staph_reference_isolates
#!/bin/sh
#$ -t 1-2646
WORKDIR=/gpfs1/well/bag/dipro/mmm/ana/Saur/jchar_denovo
INFILE=$WORKDIR/data.$SGE_TASK_ID
OUTFILE=$WORKDIR/tmp/$SGE_TASK_ID
# See comment below about paths to R
PATHTOR=/usr/local/bin/R
if [ -e $OUTFILE ]
then
rm -f $OUTFILE
fi
# Below, the phrase "EOF" marks the beginning and end of the HERE document.
# Basically, what’s going on is that we’re running R, and suppressing all of
# it’s output to STDOUT, and then redirecting whatever’s between the EOF words
# as an R script, and using variable substitution to act on the desired files.
$PATHTOR --quiet --no-save > /dev/null <<EOF
setwd("$WORKDIR/staph_reference_genomes")

ref_file = "$WORKDIR/R00000022.fa"

gene=read.table("$INFILE", header=FALSE, sep="\t")
gene_id=gene\$V1
strand=gene\$V2
ref_range=gene\$V3

genbank_genomes = read.table("$WORKDIR/ref_genomes.txt",h=F,skip=0,as.is=T,sep="\t")
genome = paste("$WORKDIR/staph_ref_genomes/", genbank_genomes[,1]", file, sep="")
cmd = paste("bash -c \"blastn -query <(cat ",genome,") -subject ",ref_file," -subject_loc ",ref_range," -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qseq'\"",sep="");

res = read.table(pipe(cmd[1]),h=FALSE,as.is=TRUE)
for(i in 2:length(cmd)) {
	res = rbind(res,read.table(pipe(cmd[i]),h=FALSE,as.is=TRUE))
	cat("Done",i,"\n")
}

# Where necessary, take the reverse complement
need.rev = res[,9]>res[,10]
revc = Vectorize(function(s) {
    v = unlist(strsplit(s,""))
    rcv = rev(c("A"="T","C"="G","G"="C","T"="A","-"="-")[v])
    paste(rcv,collapse="")
})
res[need.rev,11] = revc(res[need.rev,11])
#reverse complement all sequences for genes that are on negative strand in the reference
if (strand == "-") res[,11]=revc(res[,11])

#remove isolates that contain premature stop codons
stop = c("TAA", "TGA", "TAG")
codons = function(s){
    c = substring(s, seq(1, nchar(s)-2, 3), seq(3, nchar(s), 3))
    return(sum(as.numeric(c==stop[1] | c==stop[2] | c==stop[3])))
}
stop_counts = as.numeric(sapply(res[,11],codons))
cat(stop_counts)
res=res[as.numeric(sapply(res[,11],codons)) == 1,]

# Warn about contigs split by the BLAST search
if(any(duplicated(res[,1]))) cat(paste(res[duplicated(res[,1]),1])), file = paste("$OUTFILE.", gene_id,".split_contigs.txt", sep="")


# Construct FASTA representation
lab = paste(">",res[,1],sep="")
fa = paste(apply(matrix(c(lab,res[,ncol(res)]),ncol=2,byrow=FALSE),1,paste,collapse="\n"),collapse="\n")
# Output to temporary location
cat(fa,file="$OUTFILE.unaligned.fa")

# Run aligner using default EBI settings for nucleotide sequences (http://www.ebi.ac.uk/Tools/msa/mafft/)

mafft_cmd=paste("mafft --thread 8 --anysymbol --bl 62 --op 1.53 --ep 0.123 --reorder --retree 1 --treeout --maxiterate 0 --localpair $OUTFILE.unaligned.fa >", "$OUTFILE.",gene_id, ".aligned.mafft.fa",sep="")
system(mafft_cmd)

# Read FASTA file
#source("~/R/myutils/fasta.R")
#fa = read.fasta("/home/wilson/temp/aligned.mafft.fa")

EOF
