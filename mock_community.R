#Use previously sequenced (Winand et al. 2020) mock community Nanopore and Illumina data to compare methods 
setwd("~/University of Calgary/PostDoc/Atlantic Condor 2021/UofC_Analysis/Seaquencing/Mock_community_Winand_2020")


library(ShortRead)
library(dada2)
library(tidyverse)


#filter and trim nanopore sequences
out = filterAndTrim("Winand_2020_Long_read_SRR10391201.fastq.gz", "Winand_2020_Long_read_filt.fastq.gz", trimLeft = 100, trimRight = 100, maxLen = 1800, minLen = 1200,  truncQ = 0, compress = FALSE)

