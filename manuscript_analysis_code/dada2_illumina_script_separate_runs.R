##dada2 analysis of samples 2AT 175NW (E46), 2B-1 TinyBubbles (C18), 2A-2 Purple Haze (D52), 2B-1 Clamshell (B22), and 5-1 Kilo (E21) for Seaquencing 

#set working directory
setwd("~/University of Calgary/PostDoc/Atlantic Condor 2021/UofC_Analysis/Seaquencing/Illumina_reads/run_separation")

#load packages
library(dada2); packageVersion("dada2")
library(seqinr)


#Workflow
path = getwd()
fnFs = sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs = sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))

sample.names = sapply(strsplit(basename(fnFs), "_"), `[`, 1)

###Inspect read quality profiles

#Forward read plot
png("run202222qualityF.png") 
plotQualityProfile(fnFs[1:1])
#Turn the plotting device off.
dev.off()

#Reverse read plot
png("run202222qualityR.png") 
plotQualityProfile(fnRs[1:1])
dev.off() 


#Assign the filenames for the filtered fastq.gz files.
#Place filtered files in filtered/ subdirectory.
filtFs = file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs = file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) = sample.names
names(filtRs) = sample.names

#Double check that all the file pairs are there
if(length(fnFs) != length(fnRs)) stop("Forward and reverse files do not match!") 

#filter and trim sequences
out = filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250,240), trimLeft = c(19,20), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE)


#Check how many reads are lost. 
#Adjust 'filterandTrim' parameters as appropriate
head(out)

#save filterandtrim stats
write.csv(out, "filterandtrim_seaquencing.csv")

###Learn the Error Rates (time-consuming step)
#Needs to be done for each run separately - Run1
##############################################
filtpathF <- "~/University of Calgary/PostDoc/Atlantic Condor 2021/UofC_Analysis/Seaquencing/Illumina_reads/run_separation/run1" # CHANGE ME to the directory containing your filtered forward fastqs
filtpathR <- "~/University of Calgary/PostDoc/Atlantic Condor 2021/UofC_Analysis/Seaquencing/Illumina_reads/run_separation/run1" # CHANGE ME ...
filtFs <- list.files(filtpathF, pattern="F_filt.fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="R_filt.fastq.gz", full.names = TRUE)
sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names
set.seed(100)

errF = learnErrors(filtFs, randomize = TRUE, multithread = TRUE)
errR = learnErrors(filtRs, randomize = TRUE, multithread = TRUE)

#visualize the estimated error rates
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
    derepF <- derepFastq(filtFs[[sam]])
    ddF <- dada(derepF, err=errF, multithread=TRUE)
    derepR <- derepFastq(filtRs[[sam]])
    ddR <- dada(derepR, err=errR, multithread=TRUE)
    merger <- mergePairs(ddF, derepF, ddR, derepR)
    mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)
# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, "~/University of Calgary/PostDoc/Atlantic Condor 2021/UofC_Analysis/Seaquencing/Illumina_reads/run_separation/run1/seqtab_run1.rds") 
################################################################

#Needs to be done for each run separately - Run2
filtpathF <- "~/University of Calgary/PostDoc/Atlantic Condor 2021/UofC_Analysis/Seaquencing/Illumina_reads/run_separation/run2" # CHANGE ME to the directory containing your filtered forward fastqs
filtpathR <- "~/University of Calgary/PostDoc/Atlantic Condor 2021/UofC_Analysis/Seaquencing/Illumina_reads/run_separation/run2" # CHANGE ME ...
filtFs <- list.files(filtpathF, pattern="F_filt.fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="R_filt.fastq.gz", full.names = TRUE)
sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names
set.seed(100)

errF = learnErrors(filtFs, randomize = TRUE, multithread = TRUE)
errR = learnErrors(filtRs, randomize = TRUE, multithread = TRUE)

#visualize the estimated error rates
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
    derepF <- derepFastq(filtFs[[sam]])
    ddF <- dada(derepF, err=errF, multithread=TRUE)
    derepR <- derepFastq(filtRs[[sam]])
    ddR <- dada(derepR, err=errR, multithread=TRUE)
    merger <- mergePairs(ddF, derepF, ddR, derepR)
    mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)
# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, "~/University of Calgary/PostDoc/Atlantic Condor 2021/UofC_Analysis/Seaquencing/Illumina_reads/run_separation/run2/seqtab_run2.rds") 

##########################################################
###old
#sample inference
#dadaFs = dada(filtFs, err=errF, multithread = TRUE)
#dadaRs = dada(filtRs, err=errR, multithread = TRUE)
#saving workspace: 
#save.image(file = "seaquencing_dada2_inference.RData")
#merge forward and reverse reads
#mergers = mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
#head(mergers[[1]])
###Construct sequence table
#seqtab = makeSequenceTable(mergers)
#dim(seqtab) #Number of sequence variants
###Inspect distribution of sequence lengths
#table(nchar(getSequences(seqtab)))

#####


######################################################
###Merge runs here 
st1 <- readRDS("~/University of Calgary/PostDoc/Atlantic Condor 2021/UofC_Analysis/Seaquencing/Illumina_reads/run_separation/run1/seqtab_run1.rds")
st2 <- readRDS("~/University of Calgary/PostDoc/Atlantic Condor 2021/UofC_Analysis/Seaquencing/Illumina_reads/run_separation/run2/seqtab_run2.rds")
st.all <- mergeSequenceTables(st1, st2)
###Remove chimeras
seqtab.nochim <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)

#Sequences remaining after chimera removal.
dim(seqtab.nochim) 
#Frequency of chimeric sequences
sum(seqtab.nochim)/sum(seqtab) 


###Sequences length distribution/filtering
#Sequence variant length distribution after chimera filtering.
table(nchar(getSequences(seqtab.nochim))) 

#Select sequence variants between 251-256 bp (filter spurious #sequences) – change depending on primers used
seqtab2.nochim = seqtab.nochim[,nchar(colnames(seqtab.nochim)) %in% seq(251,256)] 

#Sequences remaining after length filtering.
dim(seqtab2.nochim) 

#Sequence variant length distribution after length filtering.
table(nchar(getSequences(seqtab2.nochim))) 


#save sequence table to CSV file 
write.csv(seqtab2.nochim, "seaquencing_ASVtable_seqs.csv")
saveRDS(seqtab2.nochim, "seaquencing_ASVtable_seqs.rds")



#track reads through steps
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
##write to csv
write.csv(track, "seaquencing_track_reads.csv")


######Sequence table manipulation
#Assign unique ID (ASV1, ASV2,...) to each sequence variant.
seqnum = paste0("ASV", seq(ncol(seqtab2.nochim))) 

#create a list of the sequences
uniqueSeqs = as.list(colnames(seqtab2.nochim))

#Transpose the matrix (ASVs in rows, samples in columns)
seqtab2.nochim.transposed = t(seqtab2.nochim) 

#Change rownames to ASV IDs.
rownames(seqtab2.nochim.transposed) = as.character(seqnum) 

#Write a .fasta file containing sequence variants.
write.fasta(uniqueSeqs, seqnum, "ASVseqs.fasta") 


#assign taxonomy
 taxa = assignTaxonomy(seqtab2.nochim, "../../silva_nr99_v138.1_train_set.fa.gz", multithread = TRUE)

#taxa2 = addSpecies(taxa, "silva_species_assignment_v138.1.fa.gz")
#memory error with addSpecies, might not be necessary anyway

#save taxonomic assignments to csv
write.csv(taxa, "seaquencing_tax_ASV.csv")

#Change rownames to ASV IDs.
rownames(taxa) = as.character(seqnum) 
write.csv(taxa, "seaquencing_tax_ASVID.csv")


###Add taxonomy to ASV table

#Confirm equal number of rows in seqtab and taxa.
nrow(seqtab2.nochim.transposed)
#nrow(taxid)
nrow(taxa)


#Combine columns from the two dataframes.
ASVseqtabtaxa = data.frame(seqtab2.nochim.transposed, taxa) 

#write to csv file
write.csv(ASVseqtabtaxa, "seaquencing_ASVseq_taxa.csv")


#saving workspace: 
save.image(file = "seaquencig_dada2_taxonomy.RData")


#################
###Make Illumina phylum level summary
setwd("~/University of Calgary/PostDoc/Atlantic Condor 2021/UofC_Analysis/Seaquencing/16S_Nanopore/Seaquencing_all_fastq_pass")
#removed extra sample manually 
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

#for calculating average abundance of "Unknown" classifications
y = abund4 %>% filter(Phylum == "Unknown") 
mean(as.matrix(y[,2:ncol(y)]))
sd(as.matrix(y[,2:ncol(y)]))

tax_df = abund4
colnames(tax_df) = gsub("JZO.2021.Condor.", "", colnames(tax_df))
colnames(tax_df) = gsub(".Univ.20220222", "", colnames(tax_df))
colnames(tax_df) = gsub("TSU.2021.Condor.", "", colnames(tax_df))
colnames(tax_df) = gsub(".20211213", "", colnames(tax_df))

colnames(tax_df)[2:ncol(tax_df)] = paste("illumina",  colnames(tax_df)[2:ncol(tax_df)], sep = "-")

write.csv(tax_df, "../../Illumina_reads/run_separation/Illumina_Phylum_summary.csv")



###Illumina phylum level bubble and bar plot
setwd("~/University of Calgary/PostDoc/Atlantic Condor 2021/UofC_Analysis/Seaquencing/Illumina_reads/run_separation")
tax_df = read.csv("Illumina_Phylum_summary.csv", header = TRUE)
colnames(tax_df) = gsub("X", "", colnames(tax_df))

tax_df_long = tax_df %>% pivot_longer(!Phylum, names_to = "Sample", values_to = "Abundance")

tax_df_long$Sample = gsub("4.8", "04.08", tax_df_long$Sample)
tax_df_long$Sample = gsub("8.12", "08.12", tax_df_long$Sample)
tax_df_long$Sample = gsub("1.4", "0.4", tax_df_long$Sample)

#colours
colours = colorRampPalette(c("#2F4858", "#33658A", "#86BBD8", "#830689", "#F5A614", "#F26419", "#BB3551",  "#C1D7AE", "#68AC5D", "#EBDDAD"))(sample_number)

xx = ggplot(tax_df_long, aes(x = Sample, y = reorder(Phylum, desc(Phylum)))) + geom_point(aes(colour = Sample, size= Abundance), alpha = 0.7) +theme(legend.key = element_blank(), legend.title = element_text(size = 10), panel.border = element_rect(fill = NA, colour = "grey80"), axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 7, angle = 90, vjust = 0.3, hjust =1), panel.background = element_blank(), panel.grid.major = element_line(colour = "grey94")) + scale_radius(range=c(1,8), breaks = c(1,10,30,50)) + labs(x = "", y = "") + scale_colour_manual(values = colours) + guides(colour = "none")
ggsave("Illumina_phylum_bubble_plot.png", height = 7, width = 7)

#bar plot 
#tax_df = illum
 tax_df$max = apply(tax_df[,2:ncol(tax_df)], 1, FUN = max, na.rm = TRUE)
 #select top 10 most abundant taxa, based on abundance in one sample
 tax_df2 <- tax_df[order(-tax_df$max),][1:10,]
 
colours = c("#2F4858", "#33658A", "#86BBD8", "#830689", "#F5A614", "#F26419", "#BB3551",  "#C1D7AE", "#68AC5D", "#EBDDAD")

tax_df2_long = tax_df2 %>% select(-max) %>% pivot_longer(!Phylum, names_to = "Sample", values_to = "Abundance")
tax_df2_long$Sample = gsub("4.8", "04.08", tax_df2_long$Sample)
tax_df2_long$Sample = gsub("8.12", "08.12", tax_df2_long$Sample)
tax_df2_long$Sample = gsub("1.4", "0.4", tax_df2_long$Sample)

gg = ggplot(tax_df2_long, aes(x = Sample, y = Abundance)) + geom_bar(aes(fill = Phylum),  position = "stack", stat = "identity", colour = "white", size = 0.1) + scale_fill_manual(values = colours) + labs(x = "", y = "Relative Abundance (%)") + theme(panel.background = element_blank(), panel.border = element_rect(fill =NA, colour = "black"), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3), legend.key = element_blank()) + scale_y_continuous(limits = c(0,100), expand = c(0,0))
ggsave("Illumina_phylum_bar_plot.png", height = 6, width =8)


#facet bar plot
#select taxa in both illumina and nano top 10
tax_df2 <- tax_df %>% select(-max) %>% filter(Phylum == "Acidobacteriota" | Phylum == "Bacteroidota" | Phylum == "Caldatribacteriota" | Phylum == "Campylobacterota" | Phylum == "Chloroflexi" | Phylum == "Cyanobacteria" | Phylum == "Desulfobacterota" | Phylum == "Latescibacterota" | Phylum == "Methylomirabilota" | Phylum == "NB1-j" | Phylum == "Planctomycetota" | Phylum == "Proteobacteria" | Phylum == "Unknown")
#add other row to df
oth = c("Other", 100 -colSums(tax_df2[,2:ncol(tax_df2)]))
test = rbind(tax_df2, oth)
tax_df2_long = test %>% pivot_longer(!Phylum, names_to = "Sample", values_to = "Abundance")
tax_df2_long$Abundance = as.numeric(tax_df2_long$Abundance)
tax_df2_long$Sample = gsub("4.8", "04.08", tax_df2_long$Sample)
tax_df2_long$Sample = gsub("8.12", "08.12", tax_df2_long$Sample)
tax_df2_long$Sample = gsub("1.4", "0.4", tax_df2_long$Sample)
#colours
colours = c("#2F4858", "#33658A", "#86BBD8", "#830689", "#F5A614", "#F26419", "#BB3551",  "#C1D7AE", "#68AC5D", "#31493C","grey80","#B07156","#EBDDAD",  "grey40")

tax_df2_long2 = tax_df2_long %>% separate(Sample, into = c("Site","Core", "Subsite", "Depth1", "Depth2"), sep = "\\.")
tax_df2_long2$Depth3 = paste0(tax_df2_long2$Depth1,"-", tax_df2_long2$Depth2)
tax_df2_long2$Subsite = gsub("Transect", "", tax_df2_long2$Subsite)

#vertical bars
gg = ggplot(tax_df2_long2, aes(x = Depth3, y = Abundance)) + geom_bar(aes(fill = Phylum),  position = "stack", stat = "identity", colour = "white", size = 0.1) + scale_fill_manual(values = colours) + labs(x = "", y = "Relative Abundance (%)") + theme(panel.background = element_blank(), panel.border = element_rect(fill =NA, colour = "black"),strip.background = element_rect(fill = "grey90", colour = "black"), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3), legend.key = element_blank()) + scale_y_continuous(limits = c(0,100), expand = c(0,0)) + facet_grid(.~Subsite, scales = "free", space = "free")

#horizontal bars
gg = ggplot(tax_df2_long2, aes(y = Depth3, x = Abundance)) + geom_bar(aes(fill = Phylum),  position = "stack", stat = "identity", colour = "white", size = 0.1) + scale_fill_manual(values = colours) + labs(y = "", x = "Relative Abundance (%)") + theme(panel.background = element_blank(), panel.border = element_rect(fill =NA, colour = "black"),strip.text.y = element_text(angle = 0), strip.background = element_rect(fill = "grey90", colour = "black"), axis.text.x = element_text(colour = "black"), legend.key = element_blank()) + scale_x_continuous(limits = c(0,100.1), expand = c(0,0)) + facet_grid(Subsite~., scales = "free", space = "free") + scale_y_discrete(limits=rev) 


#ggsave("Illumina_phylum_bar_plot_facet.png", height = 6, width =8.5)
#ggsave("Illumina_phylum_bar_plot_facet_vertical.png", height = 5.5, width =7)


