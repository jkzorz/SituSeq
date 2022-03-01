# NF-core Ampliseq for analysis of Illumina 16S reads 

Ampliseq pipeline for analysis of Illumina 16S 


```
#create nextflow conda environment
conda create -n nextflow

conda activate nextflow

#install nextflow
conda install -c bioconda nextflow 

#needed to update nextflow for some reason 
nextflow self-update


###singularity is already installed on the server!!! 
nextflow run nf-core/ampliseq -profile test,singularity
```


## primer sequences used: 

```
341F (Bact) : CCTACGGGAGGCAGCAG
785R (Bact): GACTACHVGGGTATCTAATCC
519F (Arch): CAGCMGCCGCGGTAA
915R (Arch): GTGCTCCCCCGCCAATTCCT
515F (Univ):  GTGYCAGCMGCCGCGGTAA
806R (Univ): GGACTACNVGGGTWTCTAAT
```


