#code for blast comparison of Nanopore 16S and Illumina 16S

#concatenate all fastq.gz files 
cat *.fastq.gz > Nanopore_all_seqs.fastq

#start interactive session 
salloc --mem=50G -c 20 -N 1 -n 1  -t 04:00:00

#convert to fasta format 
conda activate bbtools 
reformat.sh in=Nanopore_all_seqs.fastq out=Nanopore_all_seqs.fasta

#load blast
conda activate blast

#make a blast db of concatenated files 
makeblastdb -in Nanopore_all_seqs.fasta -out Nanopore_DB -dbtype nucl


#run blast of Illumina ASVs against Nanopore database
blastn -query ../Illumina/seaquencing_ASVseqs.fasta -db Nanopore_DB -outfmt 6 -out blast_results.tbl -max_target_seqs 1

 blastn -query ../Illumina/ASVseqs.fasta -db Nanopore_DB -outfmt 6 -out blast_results_1target_97.tbl -max_target_seqs 1 -perc_identity 97

#increase max number of target seqs to see how many Nanopore sequences each Illumina ASV hits
blastn -query ../Illumina/ASVseqs.fasta -db Nanopore_DB -outfmt 6 -out blast_results_1target_100_1000seqs.tbl -max_target_seqs 1000 -perc_identity 100

#address: -max_target_seqs 1 simply returns the first good hit found in the database, not the best hit as one would assume. Also depends on db order
#potentially can increase threshold in order to make sure best hit is found 


#remove duplicates from fasta file 
conda activate seqkit

seqkit rmdup -s < Nanopore_all_seqs.fasta > Nanopore_all_seqs_rmdup.fasta
#no sequences removed anyway - all unique 



#add short seq names to Nanopore sequences
awk '{for(x=1;x<=NF;x++)if($x~/>/){sub(/>/,">SEQ"++i" ")}}1' Nanopore_all_seqs.fasta > Nanopore_all_seqs_numbers.fasta

#redo blast with updated sequence names
makeblastdb -in Nanopore_all_seqs_numbers.fasta -out Nanopore_DB -dbtype nucl

#combine all taxonomic classifications into one large file 
cat /work/ebg_lab/gm/gapp/jzorz/Atlantic_Condor/Seaquencing/Nanopore/tax_barcode_taxonomy_files/*.csv > tax_barcode_all.csv

##example to search for Nanopore Sequence, get sequence, and then use that to search for taxonomy in tax_bacode_all.csv
 grep "SEQ633675 " Nanopore_all_seqs_numbers_singleline.fasta -A1 | tail -n 1 | grep -f - tax_barcode_all.csv



###
#In R... 
#count number of ASVs with greater >97% identity over 230 bp. 
#Need to first remove any ASV IDs that are not from Bacteria classification 
x = read.csv("blast_97_asv_compare.csv", header = TRUE)
new <- ifelse(x$Illumina_blast %in% x$Illumina_No_Arch, x$Illumina_blast, NA)
new2 = na.omit(new)
length(new2)
length(new)

#Percentage without blast hits 
 x = read.csv("blast_97_asv_compare.csv", header = TRUE)
 new <- ifelse(x$Illumina_No_Arch %in% x$Illumina_blast, NA, x$Illumina_No_Arch)
 new2 = na.omit(new)
 x2 <- x[x$Illumina_No_Arch %in% new2, -1]
 #write.csv(x2, "Blast_no_matches_Illumina.csv")
 
 #summarize ASVs without blast hits by phylum
  x3 = x2 %>% group_by(Phylum) %>% summarize(n = n())
 





