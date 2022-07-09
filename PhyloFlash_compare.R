#set wd
setwd("~/University of Calgary/PostDoc/Atlantic Condor 2021/UofC_Analysis/Seaquencing/PhyloFlash")

library(tidyverse)
library(vegan)



pf = read.csv("PhyloFlash_Relative_abundance.csv", header = TRUE)
pf2 = pf %>% group_by(Phylum,Sample) %>% summarize(Abund = sum(Abund))
pf3 = pf2 %>% pivot_wider(names_from = Sample, values_from = Abund)
colnames(pf3)[2:ncol(pf3)] = paste("PF",  colnames(pf3)[2:ncol(pf3)], sep = "-")

tax_comb = full_join(tax_df, nano, by = "Phylum")
tax_comb0 = full_join(tax_comb, pf3, by = "Phylum")
tax_comb2 = tax_comb0 %>% pivot_longer(!Phylum, names_to = "Sample", values_to = "Abundance") %>% separate(Sample, into = c("Tech", "Sample"), sep = "-")

tax_comb3[,3:5][is.na(tax_comb3[,3:5])] <- 0


tax_comb4 = tax_comb3 %>% filter(Sample == "2A2.PurpleHaze.0.4" | Sample == "2A2.PurpleHaze.12.16" | Sample == "2A2.PurpleHaze.24.28" | Sample == "2AT.175NW.0.4" | Sample == "2AT.175NW.24.28" | Sample == "2B1.TinyBubbles.0.4" | Sample == "2B1.TinyBubbles.12.16" | Sample == "2B1.TinyBubbles.20.24" | Sample == "2B1.TinyBubbles.24.30" ) 
