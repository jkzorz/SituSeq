library(tidyverse)


setwd("~/University of Calgary/PostDoc/Atlantic Condor 2021/UofC_Analysis/Seaquencing/Blast")

#read in blast results table
blast = read.csv("blast_gas_associated_results_97_1000seqs.tbl", sep = "\t", header = FALSE)


blast$V1 = gsub("_SEQ.*", "", blast$V1)
blast_hq = blast %>% filter(V11 <= 0)
blast_sum = blast_hq %>% group_by(V1) %>% summarize(n = n())

#read in file with information on number of sequences per sample
filt = read.csv("Filtered_sequence_summary.csv")

filt$X = gsub("_combined.fastq.*", "", filt$X)
filt2 = data.frame(V1 = filt$X, num_seqs = filt$reads.out)



 blast_sum2 = left_join(blast_sum, filt2, "V1")
 blast_sum2$rel_abund = (blast_sum2$n/blast_sum2$num_seqs)*100
 blast_sum2$V1 = gsub("barcode[0-9][0-9]_", "", blast_sum2$V1)
  


blast_sum2$copy = blast_sum2$V1

blast_sum3 = blast_sum2 %>% separate(copy, into = c("Site", "Subsite", "Depth"), sep = "_")
blast_sum3$Depth = gsub("2430", "2428", blast_sum3$Depth)
blast_sum3$Subsite = gsub("Transect", "", blast_sum3$Subsite)

gg = ggplot(blast_sum3, aes(y = Depth, x = rel_abund))+ facet_grid(.~Subsite) + geom_bar(stat = "identity", aes(fill = Subsite), colour = "black") + labs(y = "Depth below surface (cm)", x = "Relative abundance (%) of sequences of interest", fill = "") + theme(panel.background = element_blank(), strip.background = element_rect(fill = "grey90", colour = "black"), panel.border = element_rect(fill =NA, colour = "black"), axis.text.x = element_text( colour = "black"), axis.text.y = element_text(colour = "black"), legend.position = "none") + scale_x_continuous(expand = c(0,0), limits = c(0, max(blast_sum2$rel_abund + max(blast_sum2$rel_abund/10)))) + scale_y_discrete(limits=rev) + scale_fill_manual(values = c("#ED254E", "#023778", "#A9E2A2", "#94A2B3", "#F9DC5C"))


#ggsave("Blast_hits_per_sample_percent.png", height = 6, width = 6)
#ggsave("Blast_hits_per_sample_percent.svg", height = 6, width = 6)

