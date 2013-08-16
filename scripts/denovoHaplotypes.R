setwd("/dipro/mmm/gorm/v2/dat/Saur/velvet")
comid = c(1075:1166,1274:1284)
velvet_file = paste("/dipro/mmm/gorm/v2/dat/Saur/velvet/C0000",comid,"_n1_contigs.fa.gz",sep="")
#velvet_file = paste("/dipro/mmm/gorm/v3/dat/Saur/velvet/",formatC(1000*round(comid/1000),width=8,flag="0"),"/C0000",comid,"_contigs.fa.gz",sep="")
ref_file = "/dipro/mmm/gorm/v2/dat/Saur/refgen/R00000022.fa"

# Note that ref_range is 0-based

#ref_range = "531000-533000" # ksgA
#ref_range = "346000-351000" # SAR0304
#ref_range = "6000-12000" # gyrA
#ref_range = "126000-130000" # sirABC
#ref_range = "674000-683000" # antiporter ion resistance
# Note that ref_range is 0-based
#ref_range = "223635-226424" # SAR0196 (hsdR)
#ref_range = "1349005-1350345" # glnA
#ref_range = "2381938-2383602" # SAR2297
#ref_range = "2176993-2178608" # groEL
#ref_range = "371783-373129"#SAR0326, test frameshift handling 
ref_range = "1502887-1535127" #ebh/SAR1447
strand = "-"

cmd = paste("bash -c \"blastn -query <(zcat ",velvet_file,") -subject ",ref_file," -subject_loc ",ref_range," -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qseq'\"",sep="");

res = read.table(pipe(cmd[1]),h=FALSE,as.is=TRUE)
for(i in 2:length(cmd)) {
	res = rbind(res,read.table(pipe(cmd[i]),h=FALSE,as.is=TRUE))
	#cat("Done",i,"\n")
}

# Where necessary, take the reverse complement
need.rev = res[,9]>res[,10]
revc = Vectorize(function(s) {
    v = unlist(strsplit(s,""))
    rcv = rev(c("A"="T","C"="G","G"="C","T"="A","-"="-")[v])
    paste(rcv,collapse="")
})
res[need.rev,11] = revc(res[need.rev,11])
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
cat(res[,1])
# Warn about contigs split by the BLAST search
if(any(duplicated(res[,1]))) cat(paste("The following contigs contained discontiguous BLAST matches:",paste(res[duplicated(res[,1]),1])), file = "/home/clme1585/temp/split_contigs.txt")

# Construct FASTA representation
lab = paste(">",res[,1],sep="")
fa = paste(apply(matrix(c(lab,res[,ncol(res)]),ncol=2,byrow=FALSE),1,paste,collapse="\n"),collapse="\n")
# Output to temporary location
cat(fa,file="/home/clme1585/temp/unaligned.fa")

# Run aligner using default EBI settings for nucleotide sequences (http://www.ebi.ac.uk/Tools/msa/mafft/)
system("mafft --thread 8 --anysymbol --bl 62 --op 1.53 --ep 0.123 --reorder --retree 1 --treeout --maxiterate 0 --localpair /home/clme1585/temp/unaligned.fa > /home/clme1585/temp/aligned.mafft.fa")

# Read FASTA file
#source("~/R/myutils/fasta.R")
#fa = read.fasta("/home/clme1585/temp/aligned.mafft.fa")
