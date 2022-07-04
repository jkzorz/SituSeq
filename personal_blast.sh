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


