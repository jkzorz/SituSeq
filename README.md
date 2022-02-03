# Seaquencing
Workflows for offline analysis of 16S rRNA Nanopore data. Using a standard spec **windows 10** laptop 

# Installation 
Before running the workflow offline, the following programs need to be installed: 

1. Windows substem for linux (wsl)
2. cutadapt/or NanoFilt?
3. blast
4. Rstudio and packages *dada2* and *tidyverse* 
5. Optional: taxonomy database(s) (e.g. *Silva*: https://zenodo.org/record/4587955#.YfxAfOrMI2w )


*something to address: MinKnow automatically calculates "pass" or "fail" for reads based on average quality scores (https://bioinformatics.stackexchange.com/questions/8735/how-does-minknow-classify-1d-reads-as-pass-or-fail). Is this quality score high enough for our purposes?* 

## Install windows subsystem for linux 
Install windows subsystem for linux (wsl) (https://docs.microsoft.com/en-us/windows/wsl/install) and open from windows power shell

Open wsl (by typing ```wsl``` in powershell) and install *cutadapt* (https://cutadapt.readthedocs.io/en/stable/) 

```
#install using pip
python3 -m pip install --user --upgrade cutadapt

```

Install NanoFilt (https://github.com/wdecoster/nanofilt)

```
pip install nanofilt


NanoFilt -h
```


## Install local Blast (for windows)

Follow this link and select the *win64.exe* option for download:
https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/

Once downloaded, run the application to finish the installation.

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




## Install R, R-Studio, and packages

Install R (https://cran.rstudio.com/) and R-Studio (https://www.rstudio.com/products/rstudio/download/#download) and load *dada2* and *tidyverse* packages

```
#install tidyverse
install.packages("tidyverse")
library(tidyverse)

#install dada2 (https://benjjneb.github.io/dada2/dada-installation.html) 
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("dada2", version = "3.14")

library(dada2)

```


# Analysis Workflow

*at what point are we splitting reads into samples based on barcodes? First I assume?* 

## Step 0: Concatenate all individual fastq files belonging to one sample into one file (may need to write for loop for barcodes?)  

``` 
cat *.fastq > concat.fastq
```

## Step 1: Filter reads by length 

Use cutadapt to remove primer sequences and filter reads by length. Cutadapt can be run from a linux system, or if using a windows operating system, through windows subsystem for linux (wsl). Using a filter cutoff of between 1350 bp and 1650 bp.  

*need to decide if using the actual primer sequence, or just cutting off the number of base pairs?* 

```
#cutadapt 

#remove first 100 bps - change '-u' parameter to primer length
cutadapt -u 100 -o concat_trim1.fastq concat.fastq

#remove last 100 bps - change '-u' paramater to primer length
cutadapt -u -100 -o concat_trim2.fastq concat_trim1.fastq


#filter out reads that are less than 1350 bp or more than 1650 bp - need to determine optimum lengths... 
cutadapt -m 1350 -M 1650 -o  concat_trimmed_1350-1650.fastq concat_trim2.fastq 

```





Code backup:
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
