#Command for subsampling Nanopore 16S sequences from taxonomy files 


setwd("~/University of Calgary/PostDoc/Atlantic Condor 2021/UofC_Analysis/Seaquencing/16S_Nanopore/Seaquencing_all_fastq_pass/Seaquences-no-rarefaction")
library(tidyverse)

temp = list.files(pattern="tax.*.csv")
temp_list = list()
for (i in 1:length(temp)) {
  try({
  sample = gsub(".csv", "", temp[[i]])
  sample2 = gsub("tax.","",sample)
  new = read.csv(temp[i]) %>% filter(Kingdom == "Bacteria")
  sub = sample(1:nrow(new), 5000, replace=FALSE)
  seq2 = new[sub,]
  seq3 = seq2 %>% group_by(Phylum) %>% summarise(n = n()) %>% mutate(abund = n/(colSums(as.matrix(n)))*100) %>% select(-n)
  colnames(seq3) = c("Phylum", sample2)
  temp_list[[length(temp_list) + 1]] <- seq3
  }) 
  }


#merge all data frames in list
tax_df = temp_list %>% reduce(full_join, by='Phylum')


#write phylum summary 
write.csv(tax_df, "Phylum_summary_5000.csv")
