# Se(a)quencing
The goal of these workflows is to enable the easiest and fastest possible offline analysis of micrbiome data without sacrificing accuracy. Designed to be implemented by researchers with any level of bioinformatics experience using a standard spec laptop (**windows 10?**) 

# Installation 
Before running the workflows offline, the following programs need to be installed. Note that R is required for both stream 1 and stream 2 analyses. 

### Stream 1: Assign taxonomy using standard 16S database (e.g. Silva)
1. R, Rstudio, and packages *dada2*, *tidyverse*, and *ShortRead* 
2. Taxonomy database(s) (e.g. *Silva*: https://zenodo.org/record/4587955#.YfxAfOrMI2w )

### Stream 2: Blast sequences against custom database
1. Windows substem for linux (wsl)
2. cutadapt
3. blast
4. Optional: R, Rstudio, and package *tidyverse* for visualization of blast results 


## Install R, R-Studio, and packages (Stream 1 and 2)

Install R (https://cran.rstudio.com/) and R-Studio (https://www.rstudio.com/products/rstudio/download/#download) and load *dada2* and *tidyverse* packages

```
#install tidyverse
install.packages("tidyverse")

#install dada2 (https://benjjneb.github.io/dada2/dada-installation.html) 
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("dada2", version = "3.14")

```

## Install windows subsystem for linux (Stream 2)
Install windows subsystem for linux (wsl) (https://docs.microsoft.com/en-us/windows/wsl/install) and open from windows power shell by entering ```wsl```
 

## Install cutadapt (optional) (Stream 2)
Open wsl (by typing ```wsl``` in powershell) and install *cutadapt* (https://cutadapt.readthedocs.io/en/stable/)
*do you need to install python as well, I can't remember...*

```
#install cutadapt using pip
python3 -m pip install --user --upgrade cutadapt

```


## Install local Blast (for windows) (Stream 2)

Follow this link and select the *win64.exe* option for download:
https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/

Once downloaded, run the application to finish the installation.



# Analysis Workflow - Stream 1

Analysis workflow for assigning taxonomy to sequences in R. 


## Setup
Open RStudio, set working directory to location containing barcode sample folders, and load packages in this order to avoid masking issues
```
#set working directory
setwd("~/University of Calgary/PostDoc/Atlantic Condor 2021/UofC_Analysis/Seaquencing/16S_Nanopore/fastq_pass_combined/")

#load packages
library(ShortRead)
library(dada2)
library(tidyverse)
```

## Concatenate sequence files in R
Concatenate all individual sequence files from each sample into one combined fastq file per sample. The code expects the input file structure to resemble the default output from MinKnow with the sequences of each sample in separate folders called "barcode01", "barcode02", etc.   

```
folders <- list.files(pattern = "barcode..$" )

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
```

## Filter and trim reads using dada2 *filterandtrim* command
Use the dada2 *filterandtrim* command to remove primer sequences (trimLeft, trimRight) and filter out any sequences that are shorter or longer than expected based on the intended amplicon target (minLen, maxLen: 1200-1800 bp).  

```
#save path to object
path = getwd()

#Forward and fastq filenames have format: 
#barcode01_combined.fastq 
fnFs = sort(list.files(path, pattern="_combined.fastq.gz", full.names = TRUE))

#Extract sample names, assuming filenames have format: #samplename_XXX.fastq.gz
sample.names = sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#filter and trim reads- create new paths for new files 
filtFs <- file.path(path, "filtered", paste0(sample.names, "_filt.fastq.gz"))
names(filtFs) = sample.names

#filter and trim command - originally 1350 and 1850, try 1200 and 1800?
out = filterAndTrim(fnFs, filtFs, trimLeft = 100, trimRight = 100, maxLen = 1800, minLen = 1200,  truncQ = 0, compress = FALSE)
#see how many reads were lost
head(out,12) 

#can see distribution of sequence read lengths: 
table(nchar(getSequences('filtered/barcode01_filt.fastq.gz')))
plot(table(nchar(getSequences('filtered/barcode01_filt.fastq.gz'))))

```

## Assign taxonomy to reads using dada2 *assignTaxonomy* command
Taxonomic assignment is the most time consuming and computationally expenisve part of this workflow. To reduce time and computational costs, only a subset of 1000 sequences are used. This number can be changed in the code below if desired. This code outputs a csv file for each sample with the taxonomic information of each sequence. 

```
#for loop for getting sequences and assigning taxonomy - with subsetting to 1000 reads
for (fastq in filtFs) {
print(fastq)
seqs = getSequences(fastq)
sub = sample(1:length(seqs), 1000, replace=FALSE) 
seq2 = seqs[sub]
tax_rc = assignTaxonomy(seq2, "../../silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE, tryRC = TRUE)
base = basename(fastq)
samples = gsub("_filt.fastq.gz", "", base)
write.csv(tax_rc, paste('tax', samples, 'csv', sep = '.' ))
}
```


## Analyze and visualize results
This code creates a bubble plot of the abundances of all phyla in all samples, and a stacked bar plot of the top 10 most abundant phyla across all samples. Change the "colours" object to update the colour scheme as desired.  

```
#read in newly made csv files 
temp = list.files(pattern="tax.*.csv")
temp_list = list()
for (i in 1:length(temp)) {
    sample = gsub(".csv", "", temp[[i]])
    sample2 = gsub("tax.","",sample)
    new = read.csv(temp[i]) %>% filter(Kingdom == "Bacteria") %>% group_by(Phylum) %>% summarise(n = n()) %>% mutate(abund = n/(colSums(as.matrix(n)))*100) %>% select(-n)
    colnames(new) = c("Phylum", sample2)
temp_list[[length(temp_list) + 1]] <- new }


#merge all data frames in list
tax_df = temp_list %>% reduce(full_join, by='Phylum')

#long format
tax_df_long = tax_df %>% pivot_longer(!Phylum, names_to = "Sample", values_to = "Abundance")

#colour scheme - 12 colours for max 12 barcodes
colours = c('brown', 'red',"orange", 'gold',  'forestgreen', 'turquoise', 'lightblue', 'navy', 'purple', 'pink', 'grey', 'black')


#bubble plot
 xx = ggplot(tax_df_long, aes(x = Sample, y = reorder(Phylum, desc(Phylum)))) + geom_point(aes(colour = Sample, size= Abundance), alpha = 0.7) +theme(legend.key = element_blank(), legend.title = element_text(size = 10), panel.border = element_rect(fill = NA, colour = "grey80"), axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 8, angle = 315, vjust = 0, hjust =0), panel.background = element_blank(), panel.grid.major = element_line(colour = "grey94")) + scale_radius(range=c(1,8), breaks = c(1,10,30,50)) + labs(x = "", y = "") + scale_colour_manual(values = colours )
 xx
#save bubble plot
ggsave("bubble_plot_phyla.png", height = 6, width = 5.5)
 
#top phyla 
 tax_df$max = apply(tax_df[,2:ncol(tax_df)], 1, FUN = max, na.rm = TRUE)
 #select top 10 most abundant taxa, based on abundance in one sample
 tax_df2 <- tax_df[order(-tax_df$max),][1:10,]
 
tax_df2_long = tax_df2 %>% select(-max) %>% pivot_longer(!Phylum, names_to = "Sample", values_to = "Abundance")

#bar plot of most abundant phyla
gg = ggplot(tax_df2_long, aes(x = Sample, y = Abundance)) + geom_bar(aes(fill = Phylum), colour = "black", position = "stack", stat = "identity") + scale_fill_manual(values = colours) + labs(x = "", y = "Relative Abundance (%)") + theme(panel.background = element_blank(), panel.border = element_rect(fill =NA, colour = "black"), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3), legend.key = element_blank()) + scale_y_continuous(limits = c(0,100), expand = c(0,0))
gg
#save plot
ggsave("bar_plot_top_phyla.png", height = 5, width = 5)

```













##################################################################################################
# Code backup:
#longread UMI protocol 

conda activate longread_umi

```
cd /home/jackie/anaconda3/longread_umi

#combine all fastq sequences into one file
cat fastq_pass_UMI_Take2_28July2021/*.fastq > fastq_pass_UMI_Take2_28July2021/fastq_pass_UMI_Take2_28July2021_all.fastq 


#longread_umi pipeline command 
# - f forward adaptor sequence; -F forward primer sequence; -r reverse adaptor #sequence; -R reverse primer sequence; -s check start of read up to s bp for UMIs; -e check end of read up to e bp for UMIs
longread_umi nanopore_pipeline -d fastq_pass_UMI_Take2_28July2021/fastq_pass_UMI_Take2_28July2021_all.fastq -v 30 -o fastq_pass_UMI_Take2_28July2021/UMI_Take2_results/ -s 90 -e 90 -m 1000 -M 2500 -f CAAGCAGAAGACGGCATACGAGAT -F AGRGTTYGATYMTGGCTCAG -r AATGATACGGCGACCACCGAGATC -R GACGGGCGGTGWGTRCA -c 3 -p 1 -q r941_min_high_g330 -t 1


#longread_umi pipeline command for archaea primers 
# - f forward adaptor sequence; -F forward primer sequence; -r reverse adaptor #sequence; -R reverse primer sequence; -s check start of read up to s bp for UMIs; -e check end of read up to e bp for UMIs
longread_umi nanopore_pipeline -d fastq_pass_UMI_Take2_28July2021/fastq_pass_UMI_Take2_28July2021_all.fastq -v 30 -o fastq_pass_UMI_Take2_28July2021/UMI_Take2_results_archaea/ -s 90 -e 90 -m 1000 -M 2500 -f CAAGCAGAAGACGGCATACGAGAT -F TTCCGGTTGATCCYGCCGGA -r AATGATACGGCGACCACCGAGATC -R GGYTACCTTGTTACGACTT -c 3 -p 1 -q r941_min_high_g330 -t 1



#longread_umi pipeline command for JZ July 29 2021
longread_umi nanopore_pipeline -d JZ_29july2021_UMI_all.fastq -v 30 -o UMI_JZ_results -s 100 -e 100 -m 900 -M 2000 -f CAAGCAGAAGACGGCATACGAGAT -F AGRGTTYGATYMTGGCTCAG -r AATGATACGGCGACCACCGAGATC -R GACGGGCGGTGWGTRCA -c 3 -p 1 -q r941_min_high_g330 -t 1


#convert fastq to fasta 
sed -n '1~4s/^@/>/p;2~4p' INFILE.fastq > OUTFILE.fasta

#with nanopore seqs  
sed -n '1~4s/^@/>/p;2~4p' fastq_pass_UMI_Take2_28July2021/fastq_pass_UMI_Take2_28July2021_all.fastq  > fastq_pass_UMI_Take2_28July2021/fastq_pass_UMI_Take2_28July2021_all.fasta

cat infile.fq | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > file.fa


#quick workaround for ship
#step 1 remove first and last 100 bps of sequence to remove primers
#Step 2 filter out reads less than 1000 bp and over 1600 bp
#step 3 move to local computer and use R (dada2) to assign taxonomy using Silva database 

#remove first 100 bps
cutadapt -u 100 -o fastq_pass_UMI_Take2_28July2021/fastq_pass_UMI_Take2_28July2021_all_trimmed.fastq fastq_pass_UMI_Take2_28July2021/fastq_pass_UMI_Take2_28July2021_all.fastq

#remove last 100 bps 
cutadapt -u -100 -o fastq_pass_UMI_Take2_28July2021/fastq_pass_UMI_Take2_28July2021_all_trimmed2.fastq fastq_pass_UMI_Take2_28July2021/fastq_pass_UMI_Take2_28July2021_all_trimmed.fastq


#filter reads that are less than 1000 bp or more than 1600 bp
cutadapt -m 1000 -M 1600 -o  fastq_pass_UMI_Take2_28July2021/fastq_pass_UMI_Take2_28July2021_all_trimmed_1000-1600.fastq fastq_pass_UMI_Take2_28July2021/fastq_pass_UMI_Take2_28July2021_all_trimmed2.fastq 

=== Summary ===

Total reads processed:                  26,000
Reads with adapters:                         0 (0.0%)
Reads that were too short:              10,001 (38.5%)
Reads that were too long:                3,130 (12.0%)
Reads written (passing filters):        12,869 (49.5%)


#cluster reads at x% identity
userarch -cluster_fast query.fasta -id 0.7 --strand both -centroids centroids.fasta -uc clusters.uc


###convert to single line fasta 
awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' JZ_29july2021_UMI_all_16S_clust70.fasta > JZ_29july2021_UMI_all_16S_clust70_singleline.fasta

#change fastq to fasta file with this code: 
sed -n '1~4s/^@/>/p;2~4p' concat_trimmed_1200-1600.fastq > concat_trimmed_1200-1600.fasta

#convert fasta file into a single line fasta file (this is not necessary if your fasta files are already in single line format)
awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' concat_trimmed_1200-1600.fasta > concat_trimmed_1200-1600_singleline.fasta

#make blast database with indicator species 
makeblastdb -in Indicator_species_forDB_Carmen.fasta -out Indicator_species_forDB_Carmen -dbtype nucl

#blastn 70% identity, 2 hits max
blastn -query Take2_1000-1600_centroids70_singleline.fasta -db Indicator_species_forDB_Carmen -outfmt 6 -out Take2_indicator_species_blast.tbl -max_target_seqs 1


####can use windows subsystem for linux
Activate wsl by opening windows powershell and typing: wsl


###assign taxonomy in R using dada2
set.seed(100) # Initialize random number generator for reproducibility
fast = read.csv("Take2_28July2021_clustered70_singleline.fasta", header = FALSE, sep = ">" )
fast1 = as.data.frame(fast[,1]) 
new_fast = fast1[seq(2, nrow(fast1), 2), ]
new_fastv = as.vector(new_fast)

take2_tax_rc_132 = assignTaxonomy(new_fastv, "../silva_nr_v132_train_set.fa.gz", multithread=FALSE, tryRC = TRUE)

```


Install NanoFilt (https://github.com/wdecoster/nanofilt)

```
pip install nanofilt


NanoFilt -h
```

## Step 0: Set working directory, load packages, and concatenate all individual fastq files belonging to one sample into one file  

Set working directory in R, and load *ShortRead*, *dada2*, and *tidyverse* packages (load packages in this order to avoid masking issues). 

```
setwd("~/University of Calgary/PostDoc/Atlantic Condor 2021/UofC_Analysis/Seaquencing/16S_Nanopore/fastq_pass_combined/")
library(ShortRead)
library(dada2)
library(tidyverse)
```

Keeping the same directory and file format as the MinKnow sofware (the fastq files of each sample are in separate folders with the format barcode01, barcode 02, etc - *is this the same for all basecalling programs?*), concatenate all fastq files from the same sample into one large fastq file. 

```
#concatenating nanopore files in R
folders <- list.files(pattern = "barcode..$" )

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
```

## Step 2: Import data into R and assign taxonomy

```
#load required libraries
library(tidyverse)
library(dada2)

#Initialize random number generator for reproducibility
set.seed(100) 

#read in trimmed and filtered sequences - can use fastq or fasta files
seqs = getSequences('concat_trimmed_1200-1600_singleline.fasta')

#assign taxonomy with silva database, change location of database to match location in your machine
#this step will take a while, possibility to subset sequences  
#include tryRC parameter to account for sequences being read in opposite direction by nanopore
tax_rc = assignTaxonomy(seqs, "../silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE, tryRC = TRUE)

#out of 100 sequences, 5 had NA at kingdom level, 18 had NA at phylum level, 15 different phyla 
#out of 500 sequences, 16 (3.2%) had NA at kingdom level, 117 (23%) had NA at phylum level, 21 different phyla
#out of 500 sequences with no "RC", 32 had NA at kingdom level, 249 had NA at phylum level
#using Carmen's mock community, out of 1000 sequences, 0 had NA at phylum and kingdom levels! 
#seems that NA in sediment samples are caused by rarer taxa


#optional step to write taxonomy to csv
write.csv(tax_rc, "tax.csv")

```


## Step 3: Visualize results 

```
#Convert taxonomy matrix into a data frame
tax_df = as.data.frame(tax_rc) 

#Keep only Bacteria 
tax_df2 = tax_df %>% filter(Kingdom == "Bacteria")

#Summarize at Phylum level
phylum = tax_df2 %>% group_by(Phylum) %>% summarise(n = n())

#Calculate relative abundance
phylum$abund = (phylum$n)/(colSums(as.matrix(phylum$n)))*100

#Plot as scatter/bubble plot 
gg = ggplot(phylum, aes(y = Phylum, x = abund)) + 
    geom_point(aes(size = abund, colour = abund)) + 
    labs(x = "Relative Abundance (%)", y = "", size = "", colour = "Relative Abundance (%)") + 
    theme(legend.position = "top", legend.key = element_blank(), 
    panel.border = element_rect(fill = NA, colour = "grey80"), 
    panel.background = element_blank(), panel.grid.major = element_line(colour = "grey94"), 
    legend.title = element_text(size = 10))
gg

#Save file
ggsave("tax_500_test.png", height = 5.5, width = 7.5)

```

## Multiple samples 

## Step 0: Concatenate all individual fastq files belonging to one sample into one file  
Nanopore splits sequences into individual file chunks of a pre-determined number of sequences. It is easiest going forward if these little files are concatenated into one large file.  

```
for i in barcode*; do echo $i; cat $i/*.fastq.gz > ${i}_cat.fastq.gz;done
```

## Step 1:  Filter reads by length 




## Additional Blast search against custom database

Make blast database with sequences of interest:

```
makeblastdb -in sequences_of_interest.fasta -out custom_DB -dbtype nucl
```

Run Blastn with Nanopore 16S sequences as queries and sequences of interest as database

```
blastn -query nanopore_16S.fasta -db custom_DB -outfmt 6 -out blast_results.tbl -max_target_seqs 1

#address: -max_target_seqs 1 simply returns the first good hit found in the database, not the best hit as one would assume. Also depends on db order
#potentially can increase threshold in order to make sure best hit is found 

```
blast-qc: https://environmentalmicrobiome.biomedcentral.com/articles/10.1186/s40793-020-00361-y



### Old stuff
## Step 0: Concatenate all individual fastq files belonging to one sample into one file   
Nanopore splits sequences into individual file chunks of a pre-determined number of sequences. It is easiest going forward if these little files are concatenated into one large file.  

``` 
cat *.fastq > concat.fastq
```

## Step 1: Filter reads by length 

Use cutadapt to remove primer sequences and filter reads by length. Cutadapt can be run from a linux system, or if using a windows operating system, through windows subsystem for linux (wsl). Using a filter cutoff of between 1350 bp and 1650 bp.  

*can we use nanopore basecalling to remove barcodes and primers?*

*barcdoe seems to be ~83-87 bp based on blast hits, don't know if this also includes primer or not*

*need to decide if using the actual primer sequence, or just cutting off the number of base pairs? Also need to optimize filtering length* 

```
#cutadapt 

#remove first 100 bps - change '-u' parameter to primer length
cutadapt -u 100 -o concat_trim1.fastq concat.fastq

#remove last 100 bps - change '-u' paramater to primer length
cutadapt -u -100 -o concat_trim2.fastq concat_trim1.fastq

#filter out reads that are less than 1350 bp or more than 1650 bp - need to determine optimum lengths... 
cutadapt -m 1350 -M 1650 -o  concat_trimmed_1350-1650.fastq concat_trim2.fastq 

#Possible to do with dada2, although using minLen ended up removing all reads? Maybe a bug? 
#it was issue with "truncQ" parameter. It is optional, but automatically defaults to 2, which ended up truncating most of the reads way too soon... need to add truncQ=0 in #order to avoid this
out = filterAndTrim('apt_103_test1.fastq', 'aptseqs_trim_dada2.fastq', minLen = 600, maxLen = 1800, trimLeft = 100, trimRight = 100, truncQ = 0)
#another example, which gives the same result as the previous cutadapt commands:
 out = filterAndTrim(fwd = '../../CLI_Nanopore_test/CLI_nanopore_test.fastq', filt =  'cliseqs_trunc_test_dada2.fastq', trimLeft = 100, trimRight = 100, maxLen = 1850, minLen = 1350,  truncQ = 0, compress = FALSE)

```


## Time info - delete later 
#Time for for loop of 5 samples (all sequences): system.time()
# user   system  elapsed 
#50500.41    69.48  7178.46 

#Time for for loop of 5 samples (1000 sequence subset) ~7 minutes/sample: system.time()
# user   system  elapsed 
#12473.63    30.50  2162.59

#Time for for loop of 10 samples (1000 sequence subset) 
# user   system  elapsed 
#25717.75    34.67  4300.14 
