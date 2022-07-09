#set wd
setwd("~/University of Calgary/PostDoc/Atlantic Condor 2021/UofC_Analysis/Seaquencing/PhyloFlash")

library(tidyverse)
library(vegan)



pf = read.csv("PhyloFlash_Relative_abundance.csv", header = TRUE)
pf2 = pf %>% group_by(Phylum,Sample) %>% summarize(Abund = sum(Abund))
pf3 = pf2 %>% pivot_wider(names_from = Sample, values_from = Abund)
colnames(pf3)[2:ncol(pf3)] = paste("PF",  colnames(pf3)[2:ncol(pf3)], sep = "-")


##########################################################
####Illumina data
setwd("~/University of Calgary/PostDoc/Atlantic Condor 2021/UofC_Analysis/Seaquencing/16S_Nanopore/Seaquencing_all_fastq_pass")
asv = read.csv("../../Illumina_reads/run_separation/seaquencing_ASVseq_taxa_4analysis.csv")

library(reshape2)
library(tidyverse)

asv = asv %>% filter(Kingdom == "Bacteria")


abund = asv[,2:41]
sums = colSums(abund)
abund2 = t(t(abund)/sums*100)
abund3 = data.frame(Phylum = asv$Phylum, abund2) 
abund4 = abund3 %>% group_by(Phylum) %>% summarize_all(list(sum))
abund4$Phylum[is.na(abund4$Phylum)] <- "Unknown"



tax_df = abund4
colnames(tax_df) = gsub("JZO.2021.Condor.", "", colnames(tax_df))
colnames(tax_df) = gsub(".Univ.20220222", "", colnames(tax_df))
colnames(tax_df) = gsub("TSU.2021.Condor.", "", colnames(tax_df))
colnames(tax_df) = gsub(".20211213", "", colnames(tax_df))

colnames(tax_df)[2:ncol(tax_df)] = paste("illumina",  colnames(tax_df)[2:ncol(tax_df)], sep = "-")

#write.csv(tax_df, "../../Illumina_reads/run_separation/Illumina_Phylum_summary.csv")




###########################################


#nanopore data
nano = read.csv("Seaquences-no-rarefaction/Phylum_summary.csv",header = TRUE)
colnames(nano) = gsub("barcode[0-9][0-9]_", "", colnames(nano))
colnames(nano) = gsub("_combined", "", colnames(nano))

colnames(nano)[2:ncol(nano)] = paste("nano",  colnames(nano)[2:ncol(nano)], sep = "-")
nano$Phylum[is.na(nano$Phylum)] <- "Unknown"

#join illumina and nanopore tables

tax_comb = full_join(tax_df, nano, by = "Phylum")
tax_comb0 = full_join(tax_comb, pf3, by = "Phylum")
tax_comb2 = tax_comb0 %>% pivot_longer(!Phylum, names_to = "Sample", values_to = "Abundance") %>% separate(Sample, into = c("Tech", "Sample"), sep = "-")


#change underscores to periods in sample names 
tax_comb2$Sample = gsub("_", ".", tax_comb2$Sample)
tax_comb2$Sample = gsub("2A2.D52", "2A2", tax_comb2$Sample)
tax_comb2$Sample = gsub("2AT.E46", "2AT", tax_comb2$Sample)
tax_comb2$Sample = gsub("2B1.B22", "2B1", tax_comb2$Sample)
tax_comb2$Sample = gsub("2B1.C18", "2B1", tax_comb2$Sample)
tax_comb2$Sample = gsub("51.E21", "51", tax_comb2$Sample)
tax_comb2$Sample = gsub("Transect","", tax_comb2$Sample)
tax_comb2$Sample = gsub("0408", "4.8", tax_comb2$Sample)
tax_comb2$Sample = gsub("04", "0.4", tax_comb2$Sample)
tax_comb2$Sample = gsub("0812", "8.12", tax_comb2$Sample)
tax_comb2$Sample = gsub("1216", "12.16", tax_comb2$Sample)
tax_comb2$Sample = gsub("1620", "16.20", tax_comb2$Sample)
tax_comb2$Sample = gsub("2024", "20.24", tax_comb2$Sample)
tax_comb2$Sample = gsub("2428", "24.28", tax_comb2$Sample)
tax_comb2$Sample = gsub("2832", "28.32", tax_comb2$Sample)
tax_comb2$Sample = gsub("3236", "32.36", tax_comb2$Sample)
tax_comb2$Sample = gsub("3640", "36.40", tax_comb2$Sample)
tax_comb2$Sample = gsub("2430", "24.30", tax_comb2$Sample)
tax_comb2$Sample = gsub("51.Kilo.1.4", "51.Kilo.0.4", tax_comb2$Sample)

tax_comb3 = tax_comb2 %>% pivot_wider(names_from = Tech, values_from = Abundance)
tax_comb3[,3:5][is.na(tax_comb3[,3:5])] <- 0


tax_comb3$Phylum2 = ifelse((tax_comb3$illumina+tax_comb3$nano) > 7, tax_comb3$Phylum, "other")

tax_comb4 = tax_comb3 %>% filter(Sample == "2A2.PurpleHaze.0.4" | Sample == "2A2.PurpleHaze.12.16" | Sample == "2A2.PurpleHaze.24.28" | Sample == "2AT.175NW.0.4" | Sample == "2AT.175NW.24.28" | Sample == "2B1.TinyBubbles.0.4" | Sample == "2B1.TinyBubbles.12.16" | Sample == "2B1.TinyBubbles.20.24" | Sample == "2B1.TinyBubbles.24.30" ) 


cor(tax_comb4$PF, tax_comb4$illumina, method = "pearson")

cor(tax_comb4$PF, tax_comb4$nano, method = "pearson")

cor(tax_comb4$nano, tax_comb4$illumina, method = "pearson")






