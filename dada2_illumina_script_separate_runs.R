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
