# SituSeq
*SituSeq* is a workflow for the remote and offline analysis of Nanopore-generated 16S rRNA amplicon data. 

The first step is [**Preprocessing**](https://github.com/jkzorz/SituSeq/blob/main/Preprocessing.R) of your raw reads, which includes concatenating fastq files, removing primers, and filtering sequences for length. 

Next, there are two options: 
- [**Stream 1**](https://github.com/jkzorz/Seaquencing/blob/main/Stream1_dada2_assignTaxonomy.R): Assign taxonomy to 16S rRNA amplicon data. This method only requires R 
- [**Stream 2**](https://github.com/jkzorz/Seaquencing/blob/main/Stream2_blast_database.sh): Perform a BLAST search of your Nanopore sequences against a custom built database containing sequences of interest

*SituSeq* was designed to be implemented by researchers with any level of bioinformatics experience using a standard spec laptop. All you have to do is copy and paste the code!  

# Installation 

#### The following downloads are required for all analyses

**R:** Install the R programming language (https://cran.rstudio.com/) 

**RStudio:** Install R-Studio, an integrated development environment for easier use of the R language  (https://www.rstudio.com/products/rstudio/download/#download) 

**The R packages *tidyverse*, *ShortRead*, and *dada2*:** Install the *dada2*, *ShortRead*, and *tidyverse* packages by copying and pasting the code below

```
#install tidyverse
install.packages("tidyverse")

#install ShortRead
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ShortRead")

#install dada2 (https://benjjneb.github.io/dada2/dada-installation.html) 
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("dada2")

```

The *dada2* package requires R version 4.2 or higher. If you've already installed R, but it is an older version. Use the following code to update (for windows): 

```
install.packages("installr")
library(installr)
updateR()
```


#### Stream 1: Assign taxonomy using standard 16S database (e.g. Silva)

To perform [**Stream 1**](https://github.com/jkzorz/Seaquencing/blob/main/Stream1_dada2_assignTaxonomy.R), you will additionally need to download a taxonomy database that is compatible with the *assignTaxonomy* function from *dada2*
The **Silva database** (Select "silva_nr99_v138.1_train_set.fa.gz") is a good option: https://zenodo.org/record/4587955#.YfxAfOrMI2w 


#### Stream 2: Blast sequences against custom database

To perform [**Stream 2**](https://github.com/jkzorz/Seaquencing/blob/main/Stream2_blast_database.sh), you will additionally need to install a local copy of **BLAST**. 

If on a Windows system, I also recommend downloading the Windows Subsystem for Linux (wsl) (https://docs.microsoft.com/en-us/windows/wsl/install).
If not using a Windows system, skip this step.  

**BLAST+ installation**
Follow this link and select the option matching your operating system for download (e.g. *win64.exe* for windows):
https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/

Once downloaded, run the application to finish the installation.

