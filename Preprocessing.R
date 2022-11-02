##This code can be copied and pasted directly into R to concatenate, filter, and trim Nanopore 16S rRNA gene sequences. 
##The working directory folder should be set to the folder containing the subfolders of fastq files from each barcoded sample. E.g From the Nanopore default output, the "fastq_pass" folder would be the working directory
##Each subdirectory should start with "barcode". If you would like to add extra identifying text, please add it to the subdirectory with an underscore after "barcode". E.g. "barcode01_sample1" 
##The end result of this code will be filter and trimmed sequences and a csv file containing information on the reads in and out
##You can change the parameters in the following section before running the code 


######
#parameters to set before running:
path_to_working_directory = "." #leave as a "." if you want to set your working directory manually in RStudio "Session"--> "Set Directory" --> "Choose Directory"
#parameters for the filterAndTrim command:
minLength = 1200 #Removes reads shorter than this length. Minimum length is enforced AFTER trimming
maxLength = 1800 #Removes reads longer than this length. Maximum length is enforced BEFORE trimming 
trimLeft = 100 #The number of nucleotides to remove from the start of each read. Must cover the primer length
trimRight = 100 #The number of nucleotides to remove from the end of each read. Must cover the primer length

######
#in R
#set working directory if you haven't set it manually
setwd(path_to_working_directory)

#load packages (in this order to avoid masking issues)
library(ShortRead)
library(dada2)
library(tidyverse)

#Concatenate nanopore files in R
folders <- list.files(pattern = "barcode" )

for (directory in folders) {
print(directory) 
files = list.files(path = paste(directory, "/", sep = ""), pattern = ".fastq*") 
print(files)
fout = file.path(paste(directory, "combined.fastq.gz", sep = "_"))
    for (fl in files) {
    fq = readFastq(paste(directory,"/",fl, sep = ""))
        writeFastq(fq, fout, mode="a")
        }
}


#save path to object
path = getwd()

#Forward and fastq filenames have format: 
#barcode01_combined.fastq 
fnFs = sort(list.files(path, pattern="_combined.fastq", full.names = TRUE))

#Extract sample names, assuming filenames have format: #samplename_XXX.fastq
sample.names = sapply(strsplit(basename(fnFs), "\\."), `[`, 1)

#filter and trim reads- create new paths for new files 
filtFs <- file.path(path, "filtered", paste0(sample.names, "_filt.fastq"))
names(filtFs) = sample.names

#filter and trim command - will create new fastq files with filtered and trimmed sequences
out = filterAndTrim(fnFs, filtFs, trimLeft = trimLeft, trimRight = trimRight, maxLen = maxLength, minLen = minLength,  truncQ = 0, compress = FALSE)

#see how many reads were lost and write to csv file
head(out,12) 
write.csv(out, "Filtered_sequence_summary.csv")
