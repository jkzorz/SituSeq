###SituSeq Stream2: Search Nanopore 16S rRNA sequences against a custom database of sequences of interest

#Start with a fasta file of sequences of interest: e.g. seqs.fasta
#open your command line program (e.g. wsl in windows, Terminal on a mac, or command prompt on a linux system) 
#move into the directory containing your seqs.fasta file with the cd command 

#make a Blast database out of the fasta sequences found in seqs.fasta (change name accordingly) 

makeblastdb -in seqs.fasta -out seqs_DB -dbtype nucl

#Combine and concatenate Nanopore 16S rRNA sequences. Best to start from filter and trimmed sequences (Stream1).
#move into directory containing filter and trimmed sequences e.g. cd filtered
#convert all fastq files to fasta files
for i in *.fastq; do sed -n '1~4s/^@/>/p;2~4p' $i > $(basename $i fastq)fasta; done

#append sample name to fasta header. Fasta files with sample name in header are saved as "header_samplename"
for i in *.fasta; do sed 's/>/>"$(basename $i .fasta)"/g' $i > 'header_'$(basename $i)

#concatenate all Nanopore fasta files with changed headers together
cat header_* > NanoporeSeqs.fasta

#run blast of Nanopore 16S rRNA sequences against the custom database 

blastn -query Nanoporeseqs.fasta -db seqs_DB -outfmt 6 -out blast_results.tbl -max_target_seqs 1 -perc_identity 97




