# SituSeq
*SituSeq* is a workflow for the analysis of Nanopore-generated 16S rRNA amplicon data. The main workflow (Stream1) only requires R, and can be conducted entirely offline. It was designed to be implemented by researchers with any level of bioinformatics experience using a standard spec laptop. All you have to do is copy and paste the code!  

# Installation 
Before running the workflows offline, the following programs need to be installed. Note that R is required for both stream 1 and stream 2 analyses. 

### Stream 1: Assign taxonomy using standard 16S database (e.g. Silva)
**Downloads**
1. R, Rstudio, and packages *dada2*, *tidyverse*, and *ShortRead* 
2. Taxonomy database(s) (e.g. *Silva*: https://zenodo.org/record/4587955#.YfxAfOrMI2w )

code is here: https://github.com/jkzorz/Seaquencing/blob/main/Stream1_dada2_assignTaxonomy.R

### Stream 2: Blast sequences against custom database
**Downloads**
1. R, Rstudio, and packages *dada2* and *tidyverse*
2. Windows substem for linux (wsl) (*if using a windows machine, if not skip*)
3. Blast

code is here: https://github.com/jkzorz/Seaquencing/blob/main/Stream2_blast_database.sh

## Stream 1 and 2: Install R, R-Studio, and packages

Install R (https://cran.rstudio.com/) and R-Studio (https://www.rstudio.com/products/rstudio/download/#download) and load *dada2*, *ShortRead*, and *tidyverse* packages

```
#install tidyverse
install.packages("tidyverse")

#install ShortRead
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ShortRead")

#install dada2 (https://benjjneb.github.io/dada2/dada-installation.html) 
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("dada2", version = "3.14")

```

**Download Silva taxonomy database** (or other preferred database) for use in Stream1: https://zenodo.org/record/4587955#.YfxAfOrMI2w.
Select "silva_nr99_v138.1_train_set.fa.gz" or "silva_nr99_v138.1_wSpecies_train_set.fa.gz" to download

**Now you can run the code for Stream 1:** 
https://github.com/jkzorz/Seaquencing/blob/main/Stream1_dada2_assignTaxonomy.R

## Stream 2: Install windows subsystem for linux
Install windows subsystem for linux (wsl) (https://docs.microsoft.com/en-us/windows/wsl/install) and open from windows power shell by entering ```wsl```
If not using a Windows system, skip this step.  

## Stream 2: Install local Blast (for windows) 

Follow this link and select the *win64.exe* option for download:
https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/

Once downloaded, run the application to finish the installation.

**Now you can run the code for Stream 2:** 
code is here: https://github.com/jkzorz/Seaquencing/blob/main/Stream2_blast_database.sh
