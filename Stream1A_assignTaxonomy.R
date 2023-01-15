##This code is meant to be run after completing the Preproccessing.R code on your raw fastq files. After the Preprocessing.R code you will have filtered and trimmed fastq files for each sample in your analysis.
##This code can be copied and pasted directly into R to assign taxonomy to full length 16S rRNA reads sequenced using the Nanopore MinION platform and 16S barcode kit.
##The working directory folder should be set to the folder containing subfolders of fastq files from each barcoded sample. E.g From the Nanopore default output, the "fastq_pass" folder would be the working directory (same as for the Preprocessing step)
##The end result of this code will be a series of csv files containing the taxonomic information for each Nanopore sequence
#Use Stream 1B to visualize this data (bubble plot with abundance of all taxa, and a bar plot with the abundance of the top 10 taxa per sample) 
##You can change the parameters in the following section before running the code 

######
#parameters to set before running
subsample_depth = 1000 #each sample will be randomly subsampled to this number of reads, prior to taxonomic assignment (after filtering and trimming). For no subsampling see Nanopore_no_rarefaction.R under "backups" 
taxonomic_level = "Phylum" #choose from "Phylum" "Class" "Order" "Family" "Genus" 
sample_number = 12
path_to_taxonomy_database = "silva_nr99_v138.1_train_set.fa.gz" #change to location of taxonomy database in relation to working directory (easiest to copy taxonomy database to working directory)
path_to_working_directory = "." #leave as a "." if you want to set your working directory manually in RStudio "Session"--> "Set Directory" --> "Choose Directory"
######

#in R 
#set working directory 
setwd(path_to_working_directory)

#load packages in this order to avoid masking issues
library(ShortRead)
library(dada2)
library(tidyverse)

#save path to object
path = getwd()

#fastq filenames have format: 
#barcode01_combined.fastq 
fnFs = sort(list.files(path, pattern="_combined.fastq", full.names = TRUE))

#extract sample names, assuming filenames have format: #samplename_XXX.fastq
sample.names = sapply(strsplit(basename(fnFs), "\\."), `[`, 1)

#path for filtered and trimmed reads 
filtFs <- file.path(path, "filtered", paste0(sample.names, "_filt.fastq"))
names(filtFs) = sample.names

#import sequences and assign taxonomy - with subsetting to subsampling depth
#this will create a csv file for each sample with the sequence and its assigned taxonomy
for (fastq in filtFs) {
print(fastq)
seqs = getSequences(fastq)
sub = sample(1:length(seqs), subsample_depth, replace=FALSE) 
seq2 = seqs[sub]
tax_rc = assignTaxonomy(seq2, path_to_taxonomy_database, multithread=TRUE, tryRC = TRUE)
base = basename(fastq)
samples = gsub("_filt.fastq", "", base)
write.csv(tax_rc, paste('tax', samples, 'csv', sep = '.' ))
}

