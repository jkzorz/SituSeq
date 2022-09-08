###SituSeq Stream2: Search Nanopore 16S rRNA sequences against a custom database of sequences of interest
#This code is meant to work on the combined, filtered and trimmed sequences generated from the first part of Stream 1. 
#If using Stream 1, they should be found in the directory "filtered" 


#open your command line program (e.g. wsl in windows, Terminal on a mac or a linux system) 

#move into directory containing filter and trimmed sequences from the first part of Stream 1 e.g. cd filtered
#convert all fastq files to fasta files
for i in *.fastq; do sed -n '1~4s/^@/>/p;2~4p' $i > $(basename $i fastq)fasta; done

#append sample name to fasta header. Fasta files with sample name in header are saved as "header_samplename"
for i in *.fasta; do name="$(basename $i _combined_filt.fasta)"; sed "s/>/>${name} /g" $i > 'header_'$(basename $i); done

#concatenate all Nanopore fasta files with changed headers together
cat header_* > NanoporeSeqs.fasta

#Now you need a fasta file containing the sequences of interest: e.g. seqs.fasta
#move into the directory containing your seqs.fasta file with the cd command 
#make a Blast database out of the fasta sequences found in seqs.fasta (change name accordingly) 

makeblastdb -in seqs.fasta -out seqs_DB -dbtype nucl

#run blast of Nanopore 16S rRNA sequences against the custom database. The cutoff here is 97% identity, but can be changed accordingly 
#note that 1 max_target_seqs is not recommended when the user is interested in the "best" hit, as the blast algorithm only reports the first hit that matches the requirements 
#for this analysis, we are just interested in the samples with the most hits above 97% to our sequences of interest, and not necessarily the "best" hit.
#if using this code for a different purpose keep this in mind and change the -max_target_seqs parameter accordingly

blastn -query filtered/NanoporeSeqs.fasta -db seqs_DB -outfmt 6 -out blast_results.tbl -max_target_seqs 1 -perc_identity 97

################################################################################################################################################
### The rest of the analysis will be done in R
################################################################################################################################################
#This first section is the analysis without normalizing for number of sequences per sample 
#The following section normalizes for number of sequences per sample. Skip to the following section if you would rather calculate relative abundance of hits
#than absolute number of hits

#open R and set your working directory
#e.g. setwd("~/My_working_directory")

#load tidyverse library
library(tidyverse)

#load blast results table 
blast = read.csv("blast_results.tbl", sep = "\t", header = FALSE)

#filter out hits that are smaller than expected size
blast_hq = blast %>% filter(V4 > 400)

#or filter out hits less than a certain e-value 
blast_hq = blast %>% filter(V11 <= 0)

#Count the number of high quality hits per sample 
blast_sum = blast_hq %>% group_by(V1) %>% summarize(n = n())

#Plot the number of high quality hits per sample as bar plot
#colour scheme, change according to number of samples and desired colours
colours = c("#2F4858","#68AC5D", "#830689", "#C1D7AE",   "#86BBD8",  "#F5A614","#33658A",  "#BB3551",  "#F26419",  "#EBDDAD")

gg = ggplot(blast_sum, aes(y = V1, x = n)) + geom_bar(stat = "identity", aes(fill = V1), colour = "black") + labs(y = "", x = "Number of hits to sequences of interest", fill = "") + theme(panel.background = element_blank(), panel.border = element_rect(fill =NA, colour = "black"), axis.text.x = element_text( colour = "black"), axis.text.y = element_text(colour = "black"), legend.position = "none") + scale_x_continuous(expand = c(0,0), limits = c(0, max(blast_sum$n + max(blast_sum$n/10)))) + scale_fill_manual(values = colours)
gg

#save image
ggsave("Blast_hits_per_sample.png", height = 6, width = 6)

#If also interested in which of the sequences of interest had hits use the following code to plot. 
#If there are many sequences of interest in the database (e.g. >40) the figure will get crowded
blast_sum_seqs = blast_hq %>% group_by(V1,V2) %>% summarize(n = n())

#plot as a bar plot
gg = ggplot(blast_sum_seqs, aes(y = V2, x = n))+ facet_grid(.~V1) + geom_bar(stat = "identity", aes(fill = V1), colour = "black") + labs(y = "", x = "Number of hits to sequences of interest", fill = "") + theme(panel.background = element_blank(), panel.border = element_rect(fill =NA, colour = "black"), axis.text.x = element_text( colour = "black"), strip.background = element_rect(colour ="black"), axis.text.y = element_text(colour = "black"), legend.position = "none") + scale_x_continuous(expand = c(0,0), limits = c(0, max(blast_sum_seqs$n + max(blast_sum_seqs$n/10)))) + scale_fill_manual(values = colours)
gg

#save image
ggsave("Blast_hits_per_sample_and sequence.png", height = 8, width = 7)

#######################################################################
#This code is for calculating relative abundance of hits
#You will need the Filtered_sequence_summary.csv file from the Stream 1 Filter and Trim step

#open R and set your working directory
#e.g. setwd("~/My_working_directory")

#load tidyverse library
library(tidyverse)

#load blast results table 
blast = read.csv("blast_results.tbl", sep = "\t", header = FALSE)

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
blast_sum2 = left_join(blast_sum, filt2, "V1")

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
#If there are many sequences of interest in the database (e.g. >40) the figure will get crowded
blast_sum_seqs = blast_hq %>% group_by(V1,V2) %>% summarize(n = n())

#join the blast_sum_seqs and filt2 data frames
blast_sum_seqs2 = left_join(blast_sum_seqs, filt2, "V1")

#calculate relative abundance of sequences of interest
blast_sum_seqs2$rel_abund = (blast_sum_seqs2$n/blast_sum_seqs2$num_seqs)*100

#plot as a bar plot
 gg = ggplot(blast_sum_seqs2, aes(y = V2, x = rel_abund))+ facet_grid(.~V1) + geom_bar(stat = "identity", aes(fill = V1), colour = "black") + labs(y = "", x = "Relative abundance (%)", fill = "") + theme(panel.background = element_blank(), panel.border = element_rect(fill =NA, colour = "black"), axis.text.x = element_text( colour = "black"), strip.background = element_rect(colour ="black"), axis.text.y = element_text(colour = "black"), legend.position = "none") + scale_x_continuous(expand = c(0,0), limits = c(0, max(blast_sum_seqs2$rel_abund + max(blast_sum_seqs2$rel_abund/10)))) + scale_fill_manual(values = colours)
gg

#save image
ggsave("Blast_hits_per_sample_and sequence_percent.png", height = 8, width = 7)


