####can do everything in R now... 

##########################################
##parameters to set before running
rarefaction_depth = 1000 #for no rarefaction substitute number for: length(seqs)
taxonomic_level = "Phylum" #choose from "Phylum" "Class" "Order" "Family" "Genus" 
sample_number = 40
path_to_taxonomy_database = "../../silva_nr99_v138.1_train_set.fa.gz"
path_to_working_directory = "."


########################################

#in R - load packages in this order to avoid masking issues
setwd(path_to_working_directory)
library(ShortRead)
library(dada2)
library(tidyverse)


#try concatenating nanopore files in R
folders <- list.files(pattern = "barcode" )

for (directory in folders) {
print(directory) 
files = list.files(path = paste(directory, "/", sep = ""), pattern = ".fastq") 
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
fnFs = sort(list.files(path, pattern="_combined.fastq.gz", full.names = TRUE))

#Extract sample names, assuming filenames have format: #samplename_XXX.fastq.gz
sample.names = sapply(strsplit(basename(fnFs), "\\."), `[`, 1)


#filter and trim reads- create new paths for new files 
filtFs <- file.path(path, "filtered", paste0(sample.names, "_filt.fastq.gz"))
names(filtFs) = sample.names


#filter and trim command - compress true or false?, originally 1350 and 1850, try 1200 and 1800?
out = filterAndTrim(fnFs, filtFs, trimLeft = 100, trimRight = 100, maxLen = 1800, minLen = 1200,  truncQ = 0, compress = FALSE)
#see how many reads were lost
head(out,12) 
write.csv(out, "Filtered_sequence_summary.csv")



#for loop for getting sequences and assigning taxonomy - with subsetting to rarefaction depth
for (fastq in filtFs) {
print(fastq)
seqs = getSequences(fastq)
sub = sample(1:length(seqs), rarefaction_depth, replace=FALSE) 
seq2 = seqs[sub]
tax_rc = assignTaxonomy(seq2, path_to_taxonomy_database, multithread=TRUE, tryRC = TRUE)
base = basename(fastq)
samples = gsub("_filt.fastq.gz", "", base)
write.csv(tax_rc, paste('tax', samples, 'csv', sep = '.' ))
}


#read in newly made csv files 
temp = list.files(pattern="tax.*.csv")
temp_list = list()
for (i in 1:length(temp)) {
    sample = gsub(".csv", "", temp[[i]])
    sample2 = gsub("tax.","",sample)
    new = read.csv(temp[i], header = TRUE) 
    new2 = new %>% filter(Kingdom == "Bacteria") %>% select(all_of(taxonomic_level)) %>% group_by_all() %>% summarise(n = n()) %>% mutate(abund = n/(colSums(as.matrix(n)))*100) %>% select(-n)
    colnames(new2) = c(taxonomic_level, sample2)
    temp_list[[length(temp_list) + 1]] <- new2 }


#merge all data frames in list
tax_df = temp_list %>% reduce(full_join, by=taxonomic_level)

#write phylum summary csv
write.csv(tax_df, paste0(taxonomic_level,"_summary.csv"))

#long format
tax_df_long = tax_df %>% pivot_longer(!taxonomic_level, names_to = "Sample", values_to = "Abundance")

#colour scheme 
colours = colorRampPalette(c('brown', 'red',"orange", 'gold',  'forestgreen', 'turquoise', 'dodgerblue', 'navy', 'purple', 'pink',  'black'))(sample_number)

#remove "_combined" from sample name
tax_df_long$Sample = gsub("_combined","",tax_df_long$Sample)

#bubble plot
xx = ggplot(tax_df_long, aes(x = Sample, y = reorder(get(taxonomic_level), desc(get(taxonomic_level))))) + geom_point(aes(colour = Sample, size= Abundance), alpha = 0.7) +theme(legend.key = element_blank(), legend.title = element_text(size = 10), panel.border = element_rect(fill = NA, colour = "grey80"), axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 7, angle = 90, vjust = 0.3, hjust =1), panel.background = element_blank(), panel.grid.major = element_line(colour = "grey94")) + scale_radius(range=c(1,8), breaks = c(1,10,30,50)) + labs(x = "", y = "", colour = taxonomic_level) + scale_colour_manual(values = colours) + guides(colour = "none")
 xx
 #save bubble plot
 ggsave(paste0("bubble_plot_",taxonomic_level,".png"), height = 6, width = 5.5)
 
 #top phyla 
 tax_df$max = apply(tax_df[,2:ncol(tax_df)], 1, FUN = max, na.rm = TRUE)
 #select top 10 most abundant taxa, based on abundance in one sample
 tax_df2 <- tax_df[order(-tax_df$max),][1:10,]
 
 colours = colorRampPalette(c('brown', 'red',"orange", 'gold',  'forestgreen', 'turquoise', 'dodgerblue', 'navy', 'purple', 'pink',  'black'))(10)

 
 tax_df2_long = tax_df2 %>% select(-max) %>% pivot_longer(!taxonomic_level, names_to = "Sample", values_to = "Abundance")

tax_df2_long$Sample = gsub("_combined","",tax_df2_long$Sample)

#bar plot of most abundant phyla
gg = ggplot(tax_df2_long, aes(x = Sample, y = Abundance)) + geom_bar(aes(fill = get(taxonomic_level)), colour = "black", position = "stack", stat = "identity") + scale_fill_manual(values = colours) + labs(x = "", y = "Relative Abundance (%)", fill = taxonomic_level) + theme(panel.background = element_blank(), panel.border = element_rect(fill =NA, colour = "black"), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3), legend.key = element_blank()) + scale_y_continuous(limits = c(0,100), expand = c(0,0))
gg
#save plot
ggsave(paste0("bar_plot_top_",taxonomic_level,".png"), height = 6, width = 5)





