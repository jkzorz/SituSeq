
###SituSeq Stream 2B: Summarize and visualize BLAST search results from Stream 2A
#This code is meant to work on the BLAST results table generated from the Stream 2A code. 
#The end result of this code is a figure with the number of hits or percent abundance of the sequences of interest per sample 
#There is an option for absolute number of hits or percent abundance using the number of reads per sample
#There is also an option to sum the hits of all the sequences of interest, or to show the abundance of each sequence separately in each sample
#Choose the same working directory as Preprocessing script (Folder containing the "filtered" folder from Preprocessing) and Stream 2A

######
#parameters to set before running:
path_to_working_directory = "." #leave as a "." if you want to set your working directory manually in RStudio "Session"--> "Set Directory" --> "Choose Directory"
sample_number = 12 #number of samples 
alignment_length = 400 #minimum length of alignment between Nanopore sequence and database sequence to keep 
e_value = 0.01 #maximum e-value of BLAST match to keep
#######


#set working directory
setwd(path_to_working_directory)


#load packages
library(tidyverse)


##Processing BLAST results##

#load blast results table 
blast = read.csv("blast_results.csv",  header = FALSE)

#filter out hits that are smaller than expected size
blast_hq = blast %>% filter(V4 > alignment_length)

#or filter out hits less than a certain e-value 
blast_hq = blast %>% filter(V11 <= e_value)

#count the number of high quality hits per sample 
blast_sum = blast_hq %>% group_by(V1) %>% summarize(n = n())

#plot the number of high quality hits per sample as bar plot
#colour scheme, change according to number of samples and desired colours
colours = colorRampPalette(c("#2F4858", "#33658A", "#86BBD8", "#830689", "#F5A614", "#F26419", "#BB3551",  "#C1D7AE", "#68AC5D", "#EBDDAD"))(sample_number)

#bar plot
gg = ggplot(blast_sum, aes(y = V1, x = n)) + geom_bar(stat = "identity", aes(fill = V1), colour = "black") + labs(y = "", x = "Number of hits to sequences of interest", fill = "") + theme(panel.background = element_blank(), panel.border = element_rect(fill =NA, colour = "black"), axis.text.x = element_text( colour = "black"), axis.text.y = element_text(colour = "black"), legend.position = "none") + scale_x_continuous(expand = c(0,0), limits = c(0, max(blast_sum$n + max(blast_sum$n/10)))) + scale_fill_manual(values = colours)
gg

#save image
ggsave("Blast_hits_per_sample.png", height = 6, width = 6)

#If also interested in which of the sequences of interest had hits use the following code to plot. 
#If there are many sequences of interest in the database (e.g. >30) the figure will get crowded
blast_sum_seqs = blast_hq %>% group_by(V1,V2) %>% summarize(n = n())

#plot as a bar plot
gg = ggplot(blast_sum_seqs, aes(y = V2, x = n))+ facet_grid(.~V1) + geom_bar(stat = "identity", aes(fill = V1), colour = "black") + labs(y = "", x = "Number of hits to sequences of interest", fill = "") + theme(panel.background = element_blank(), panel.border = element_rect(fill =NA, colour = "black"), axis.text.x = element_text( colour = "black"), strip.background = element_rect(colour ="black"), axis.text.y = element_text(colour = "black"), legend.position = "none") + scale_x_continuous(expand = c(0,0), limits = c(0, max(blast_sum_seqs$n + max(blast_sum_seqs$n/10)))) + scale_fill_manual(values = colours)
gg

#save image
ggsave("Blast_hits_per_sample_and sequence.png", height = 8, width = 7)

#######################################################################
#This code is for calculating relative abundance of hits
#You will need the Filtered_sequence_summary.csv file from the Preprocessing.R script filterAndTrim step

#open R and set your working directory
#e.g. setwd("~/My_working_directory")

#load tidyverse library
library(tidyverse)

#load blast results table 
blast = read.csv("blast_results.csv", header = FALSE)

#filter out hits that are smaller than expected size
blast_hq = blast %>% filter(V4 > 400)

#or filter out hits less than a certain e-value 
blast_hq = blast %>% filter(V11 <= 0)

#Count the number of high quality hits per sample 
blast_sum = blast_hq %>% group_by(V1) %>% summarize(n = n())

#load in Filtered_sequence_summary.csv file from the Stream 1 Filter and Trim step
filt = read.csv("Filtered_sequence_summary.csv")

#remove file extension from first column 
filt$X = gsub("_combined.fastq.*", "", filt$X)

#select the columns with the sample name and with the final sequence count for each sample
filt2 = data.frame(V1 = filt$X, num_seqs = filt$reads.out)

#join the blast_sum and filt2 data frames
blast_sum2 = full_join(blast_sum, filt2, "V1")

#calculate relative abundance of sequences of interest
blast_sum2$rel_abund = (blast_sum2$n/blast_sum2$num_seqs)*100

#plot relative abundance of sequences of interest per sample as a bar plot
#colour scheme, change according to number of samples and desired colours
colours = c("#2F4858","#68AC5D", "#830689", "#C1D7AE",   "#86BBD8",  "#F5A614","#33658A",  "#BB3551",  "#F26419",  "#EBDDAD")

gg = ggplot(blast_sum2, aes(y = V1, x = rel_abund)) + geom_bar(stat = "identity", aes(fill = V1), colour = "black") + labs(y = "", x = "Relative abundance (%) of sequences of interest", fill = "") + theme(panel.background = element_blank(), panel.border = element_rect(fill =NA, colour = "black"), axis.text.x = element_text( colour = "black"), axis.text.y = element_text(colour = "black"), legend.position = "none") + scale_x_continuous(expand = c(0,0), limits = c(0, max(blast_sum2$rel_abund + max(blast_sum2$rel_abund/10)))) + scale_fill_manual(values = colours)
gg

#save image
ggsave("Blast_hits_per_sample_percent.png", height = 6, width = 6)


#If also interested in which of the sequences of interest had hits use the following code to plot. 
#If there are many sequences of interest in the database (e.g. >30) the figure will get crowded
blast_sum_seqs = blast_hq %>% group_by(V1,V2) %>% summarize(n = n())

#join the blast_sum_seqs and filt2 data frames
blast_sum_seqs2 = full_join(blast_sum_seqs, filt2, "V1")

#calculate relative abundance of sequences of interest
blast_sum_seqs2$rel_abund = (blast_sum_seqs2$n/blast_sum_seqs2$num_seqs)*100

#plot as a bar plot
gg = ggplot(blast_sum_seqs2, aes(y = V2, x = rel_abund))+ facet_grid(.~V1) + geom_bar(stat = "identity", aes(fill = V1), colour = "black") + labs(y = "", x = "Relative abundance (%)", fill = "") + theme(panel.background = element_blank(), panel.border = element_rect(fill =NA, colour = "black"), axis.text.x = element_text( colour = "black"), strip.background = element_rect(colour ="black"), axis.text.y = element_text(colour = "black"), legend.position = "none") + scale_x_continuous(expand = c(0,0), limits = c(0, max(blast_sum_seqs2$rel_abund + max(blast_sum_seqs2$rel_abund/10)))) + scale_fill_manual(values = colours)
gg

#save image
ggsave("Blast_hits_per_sample_and_sequence_percent.png", height = 8, width = 7)
