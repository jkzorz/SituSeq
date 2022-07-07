#Use previously sequenced (Winand et al. 2020) mock community Nanopore and Illumina data to compare methods 
setwd("~/University of Calgary/PostDoc/Atlantic Condor 2021/UofC_Analysis/Seaquencing/Mock_community_Winand_2020")


library(ShortRead)
library(dada2)
library(tidyverse)


#filter and trim nanopore sequences
out = filterAndTrim("Winand_2020_Long_read_SRR10391201.fastq.gz", "Winand_2020_Long_read_filt.fastq.gz", trimLeft = 100, trimRight = 100, maxLen = 1800, minLen = 1200,  truncQ = 0, compress = FALSE)
head(out,12) 

write.csv(out, "Filtered_sequence_summary.csv")

#assign taxonomy
#read in filtered fastq
seqs = getSequences("Winand_2020_Long_read_filt.fastq.gz")
#subset filtered fastq file to 5000 reads
sub = sample(1:length(seqs), 5000, replace=FALSE)
seq2 = seqs[sub]
rm(seqs)

#assign taxonomy
tax_rc = assignTaxonomy(seq2, "../silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE, tryRC = TRUE)
write.csv(tax_rc, "Mock_community_nanopore_taxonomy.csv")

#try adding species
taxa <- addSpecies(tax_rc, "../silva_species_assignment_v138.1.fa.gz")
write.csv(taxa, "Mock_community_nanopore_taxonomy_with_species.csv")


###Illumina sequences
library(dada2); packageVersion("dada2")
library(seqinr)


