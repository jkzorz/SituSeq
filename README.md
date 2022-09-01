# SituSeq
*SituSeq* is a workflow for the analysis of Nanopore-generated 16S rRNA amplicon data. The main workflow only requires R, and can be conducted entirely offline. It was designed to be implemented by researchers with any level of bioinformatics experience using a standard spec laptop. All you have to do is copy and paste the code!  

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
3. blast

code is here: https://github.com/jkzorz/Seaquencing/blob/main/Stream2_blast_database.sh

## Stream 1 and 2: Install R, R-Studio, and packages

Install R (https://cran.rstudio.com/) and R-Studio (https://www.rstudio.com/products/rstudio/download/#download) and load *dada2* and *tidyverse* packages

```
#install tidyverse
install.packages("tidyverse")

#install dada2 (https://benjjneb.github.io/dada2/dada-installation.html) 
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("dada2", version = "3.14")

```

## Stream 2: Install windows subsystem for linux
Install windows subsystem for linux (wsl) (https://docs.microsoft.com/en-us/windows/wsl/install) and open from windows power shell by entering ```wsl```
 

## Stream 2: Install local Blast (for windows) 

Follow this link and select the *win64.exe* option for download:
https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/

Once downloaded, run the application to finish the installation.
