# SituSeq
*SituSeq* is a workflow for the remote and offline analysis of Nanopore-generated 16S rRNA amplicon data. 

- [**Preprocessing**](https://github.com/jkzorz/SituSeq/blob/main/Preprocessing.R): The first step of the workflow is to preprocess your raw reads, which includes concatenating fastq files, removing primers, and filtering sequences for length. Preprocessing only requires R. 

Next, there are two options: 
- [**Stream 1**](https://github.com/jkzorz/SituSeq/blob/main/Stream1A_assignTaxonomy.R): Assign taxonomy to 16S rRNA amplicon data. This method only requires R. Use [Stream1B](https://github.com/jkzorz/SituSeq/blob/main/Stream1B_visualizeTaxonomy.R) for summary and visualization of the taxonomic classification.  
- [**Stream 2**](https://github.com/jkzorz/SituSeq/blob/main/Stream2_BLAST_database_search.R): Perform a BLAST search of your Nanopore sequences against a custom built database containing sequences of interest. This method only requires R.

*SituSeq* was designed to be implemented by researchers with any level of bioinformatics experience using a standard spec laptop. All you have to do is copy and paste the code!  

# Installation 

### The following downloads are required for all analyses (Preprocessing, Stream1, and Stream2)

**R:** Install the [R programming language](https://cran.rstudio.com/) 

**RStudio:** Install [R-Studio](https://www.rstudio.com/products/rstudio/download/#download), an integrated development environment for easier use of the R language 

**The R packages *tidyverse*, *ShortRead*, and *dada2*:** Install these packages by copying and pasting the code below in R

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


### Stream 1: Assign taxonomy using standard 16S database (e.g. Silva)

To perform [**Stream 1**](https://github.com/jkzorz/SituSeq/blob/main/Stream1A_assignTaxonomy.R), you will additionally need to download a taxonomy database that is compatible with the *assignTaxonomy* function from *dada2*.

The [**Silva database**](https://zenodo.org/record/4587955#.YfxAfOrMI2w ) (Select "silva_nr99_v138.1_train_set.fa.gz") is a good option.


### Stream 2: Blast sequences against custom database

To perform [**Stream 2**](https://github.com/jkzorz/SituSeq/blob/main/Stream2_BLAST_database_search.R), you will need to install **rBLAST**

```
if (!require("BiocManager", quietly = TRUE))
  +     install.packages("BiocManager")

BiocManager::install("Biostrings")
install.packages('rBLAST', repos = 'https://mhahsler.r-universe.dev')
```


