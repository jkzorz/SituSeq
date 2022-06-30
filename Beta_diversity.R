##Nanopore alone

setwd("C:/Users/jacqu/Documents/University of Calgary/PostDoc/Atlantic Condor 2021/UofC_Analysis/Seaquencing/16S_Nanopore/Seaquencing_all_fastq_pass/Seaquences-no-rarefaction")

library(tidyverse)
library(vegan)

#load in phylum summary
tax_df = read.csv("Phylum_summary.csv", header = TRUE)

tax_df[is.na(tax_df)] <- 0
tax_dft = t(tax_df[,-1])
colnames(tax_dft) = tax_df[,1]
m_com = as.matrix(tax_dft[,2:ncol(tax_dft)])

set.seed(123)
nmds = metaMDS(m_com, distance = "bray")
nmds
#Stress:     0.07474386 

data.scores = as.data.frame(scores(nmds)$sites)
data.scores$Sample = colnames(tax_df[,2:ncol(tax_df)])

data.scores$Sample = gsub("_04_", "_00_",data.scores$Sample)
data.scores$Sample = gsub("barcode[0-9][0-9]_", "",data.scores$Sample)
data.scores$Sample = gsub("_combined", "",data.scores$Sample)

data.scores2 = data.scores %>% separate(Sample, into = c("Site", "Subsite", "Depth"), sep = "_")
data.scores2$Depth = gsub("0408", "4", data.scores2$Depth)
data.scores2$Depth = gsub("0812", "8", data.scores2$Depth)
data.scores2$Depth = gsub("1216", "12", data.scores2$Depth)
data.scores2$Depth = gsub("1620", "16", data.scores2$Depth)
data.scores2$Depth = gsub("2024", "20", data.scores2$Depth)
data.scores2$Depth = gsub("2428", "24", data.scores2$Depth)
data.scores2$Depth = gsub("2832", "28", data.scores2$Depth)
data.scores2$Depth = gsub("3236", "32", data.scores2$Depth)
data.scores2$Depth = gsub("3640", "36", data.scores2$Depth)
data.scores2$Depth = gsub("2430", "24", data.scores2$Depth)





