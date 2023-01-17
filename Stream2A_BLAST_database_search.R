###SituSeq Stream2: Search Nanopore 16S rRNA sequences against a custom database of sequences of interest
#Note that the Stream 2A BLAST search can not be run from a path (list of directories leading to your working directory) that contains any space characters
#This code is meant to work on the combined, filtered and trimmed sequences generated from the Preproccessing.R code. 
#The filtered and trimmed sequences should be found in the directory "filtered" 
#The end result of this code is a table with BLAST search results 
#There is an option for absolute number of hits or percent abundance using the number of reads per sample
#There is also an option to sum the hits of all the sequences of interest, or to show the abundance of each sequence separately in each sample
#Choose the same working directory as Preprocessing script (Folder containing the "filtered" folder from Preprocessing)

######
#parameters to set before running:
path_to_working_directory = "." #leave as a "." if you want to set your working directory manually in RStudio "Session"--> "Set Directory" --> "Choose Directory"
path_to_seqs4_BLAST_db = "BLAST_DB.fasta" #change to file name and location of sequences to use for BLAST database
######

#load packages
library('rBLAST')
library(ShortRead)
library(tidyverse)

#set working directory
setwd(paste0(path_to_working_directory,"/filtered"))

#create fasta files from filtered fastq files, and change headers to contain sample names
k = list.files(pattern = "_filt.fastq", full.names = TRUE)
for (fq in k) {
  fq2 = basename(fq)
  fa = gsub("fastq", "fasta",fq2)
  sample = gsub("_combined_filt.fastq", "", fq2)
  print(fq2)
  writeFasta(readFastq(fq), fa)
  seq <- readDNAStringSet(fa)
  seq@ranges@NAMES = paste0(sample, " ", seq@ranges@NAMES)
  #overwrites fasta file with new headers
  writeFasta(seq, fa)
}

#concatenate fasta files together into "All_combined_seqs.fasta"
m <- list.files(pattern = "_filt.fasta" )
for (fa in m) {
  print(fa) 
  fout = file.path("../All_combined_seqs.fasta")
  writeFasta(readFasta(fa), fout, mode = "a")
}

#start blast steps
setwd("../")
seq <- readDNAStringSet("All_combined_seqs.fasta")

#Make BLAST db and perform BLAST search 
makeblastdb(path_to_seqs4_BLAST_db, dbtype = "nucl")
dbb <- blast(db=path_to_seqs4_BLAST_db)
#change parameters here as required
results = predict(dbb, seq, BLAST_args= c("-max_target_seqs 1"))
write.csv(results, "blast_results.csv", row.names = FALSE)

