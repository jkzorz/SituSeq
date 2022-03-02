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


Ampliseq batch command used:

The --skip-ancom parameter didn't seem to result in ancom skipping... 

```
#!/bin/bash
###### Reserve computing resources ######
#SBATCH --mail-user=jacqueline.zorz@ucalgary.ca
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=180GB
#SBATCH --time=16:00:00
#SBATCH --partition=bigmem,cpu2019,cpu2021,cpu2021-bf24

###### Set environment variables ######
echo "Starting run at : 'date'"
source /home/jacqueline.zorz/software/miniconda3/etc/profile.d/conda.sh
conda activate nextflow

cd /work/ebg_lab/gm/gapp/jzorz/Atlantic_Condor/UofC/

###### Run your script ######

nextflow run nf-core/ampliseq --input "/work/ebg_lab/gm/gapp/jzorz/Atlantic_Condor/UofC/" --metadata "/work/ebg_lab/gm/gapp/jzorz/Atlantic_Condor/UofC/condor_metadata3.txt" --FW_primer "GTGYCAGCMGCCGCGGTAA" --RV_primer "GGACTACNVGGGTWTCTAAT" --skip-ancom -profile singularity

##
echo "Job finished with exit code $? at: 'date'"
##
```
