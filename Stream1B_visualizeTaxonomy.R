##This code is meant to be run after completing the Preproccessing.R and Stream1A_assignTaxonomy.R code on your raw fastq files. After the Stream1A_assignTaxonomy.R code you will have a csv file for each sample containing taxonomic information for the Nanopore sequences. 
##This code can be copied and pasted directly into R to summarize visualize the taxonomy assigned in Stream1A_assignTaxonomy.R
##The working directory folder should be set to the folder containing subfolders of fastq files from each barcoded sample. E.g From the Nanopore default output, the "fastq_pass" folder would be the working directory (same as for the Preprocessing and Stream1A steps)
##The end result of this code will be a bubble plot with abundance of all taxa, and a bar plot with the abundance of the top 10 taxa per sample 
##You can change the parameters in the following section before running the code 


######
#parameters to set before running
taxonomic_level = "Phylum" #choose from "Phylum" "Class" "Order" "Family" "Genus" 
sample_number = 12
path_to_working_directory = "." #leave as a "." if you want to set your working directory manually in RStudio "Session"--> "Set Directory" --> "Choose Directory"
######

#in R 
#set working directory 
setwd(path_to_working_directory)

#load tidyverse package
library(tidyverse)


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

#remove "_combined" from sample name
colnames(tax_df) = gsub("_combined","",colnames(tax_df))

#write summary csv of taxonomic level 
write.csv(tax_df, paste0(taxonomic_level,"_summary.csv"), row.names = FALSE)

#convert data to long format
tax_df_long = tax_df %>% pivot_longer(!taxonomic_level, names_to = "Sample", values_to = "Abundance")

#colour scheme for bubble plot
colours = colorRampPalette(c("#2F4858", "#33658A", "#86BBD8", "#830689", "#F5A614", "#F26419", "#BB3551",  "#C1D7AE", "#68AC5D", "#EBDDAD"))(sample_number)


#bubble plot
xx = ggplot(tax_df_long, aes(x = Sample, y = reorder(get(taxonomic_level), desc(get(taxonomic_level))))) + geom_point(aes(colour = Sample, size= Abundance), alpha = 0.7) +theme(legend.key = element_blank(), legend.title = element_text(size = 10), panel.border = element_rect(fill = NA, colour = "grey80"), axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 7, angle = 90, vjust = 0.3, hjust =1), panel.background = element_blank(), panel.grid.major = element_line(colour = "grey94")) + scale_radius(range=c(1,8), breaks = c(1,10,30,50)) + labs(x = "", y = "", colour = taxonomic_level) + scale_colour_manual(values = colours) + guides(colour = "none")
xx

#save bubble plot
ggsave(paste0("bubble_plot_",taxonomic_level,".png"), height = 9, width = 5.5)

#select top 10 most abundant taxa, based on maximum abundance in data set
tax_df$max = apply(tax_df[,2:ncol(tax_df)], 1, FUN = max, na.rm = TRUE)
tax_df2 <- tax_df[order(-tax_df$max),][1:10,]
 
#colour scheme for bar plot
colours = c("#2F4858", "#33658A", "#86BBD8", "#830689", "#F5A614", "#F26419", "#BB3551",  "#C1D7AE", "#68AC5D", "#EBDDAD")

#convert data to long format
tax_df2_long = tax_df2 %>% select(-max) %>% pivot_longer(!taxonomic_level, names_to = "Sample", values_to = "Abundance")

#remove "_combined" from sample name
tax_df2_long$Sample = gsub("_combined","",tax_df2_long$Sample)

#bar plot of most abundant phyla
gg = ggplot(tax_df2_long, aes(x = Sample, y = Abundance)) + geom_bar(aes(fill = get(taxonomic_level)), colour = "black", position = "stack", stat = "identity") + scale_fill_manual(values = colours) + labs(x = "", y = "Relative Abundance (%)", fill = taxonomic_level) + theme(panel.background = element_blank(), panel.border = element_rect(fill =NA, colour = "black"), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3), legend.key = element_blank()) + scale_y_continuous(limits = c(0,100), expand = c(0,0))
gg

#save plot
ggsave(paste0("bar_plot_top_",taxonomic_level,".png"), height = 6, width = 5)
