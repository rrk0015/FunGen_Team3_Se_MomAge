#! /bin/bash

######### Script 1: Download SRA (Sequence Read Archives) files and Perform Preliminary FastQc ############
## Purpose: The purpose of this script is to
##	make a temporary scratch directory to read in raw data SRA sequences 
##	download data from NCBI SRA using the SRAtoolkit and the SRA run IDs: https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/
##	use FASTQC to evaluate the quality of the data: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

## Download from SRA Pipepline: 
##  Input Data: N/A
##  Output: Downloaded read files, R1 and R2 files for each sample if paired-end data (FASTQ)

## FASTQC Pipeline:
##	Input: Downloaded SRA files .fastq
##  Output: Folder for each input .fastq SRA file; The last line of this script will make a tarball of the output directory (PreCleanQuality) to bring back to your computer

## For running the script on the Alabama Super Computer.
##	For more information: https://hpcdocs.asc.edu/content/slurm-queue-system
##	After you have this script in your ASC directory, Download_qc, and you have made it executable using "chmod +x [script name]",
##	then run the script by using "run_script [script name]"
##	suggested paramenters are below to submit this script.
##              queue: class
##              core: 1
##              time limit (HH:MM:SS): 04:00:00
##              Memory: 4gb
##              run on dmc
###############################################


########## Load Modules
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load sra
module load fastqc/0.10.1

##########  Define variables and make directories
## Replace the numbers in the brackets with Your specific information
  ## make variable for your ASC ID so the directories are automatically made in YOUR directory
MyID=team3_finalproj          ##For this project, scripts were run on the ASC in a joint directory with other team projects, hence the specificity for team 3. 

## Make variables that represent YOUR working directory (WD) and your Raw data directory (DD) in scratch, and the precleaned status directory (CS) to be tar-balled at the end.
DD=/scratch/$MyID/PracticeRNAseq/RawData                        ## Example: DD=/scratch/$MyID/PracticeRNAseq/RawData
WD=/scratch/$MyID/PracticeRNAseq                                ## Example: WD=/scratch/$MyID/PracticeRNAseq
CS=PreCleanQuality

##  make the directories in SCRATCH for holding the raw data
## -p tells it to make any upper level directories that are not there. Notice how this will also make the WD.
mkdir -p $DD
## move to the Data Directory, where you will be downloading your SRA sequence files 
cd $DD

##########  Download data files from NCBI: SRA using the Run IDs
  ### from SRA use the SRA tool kit - see NCBI website https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/
        ## this downloads the SRA file and converts to fastq
        ## -F   Defline contains only original sequence name.
        ## -I   Append read id after spot id as 'accession.spot.readid' on defline.
        ## splits the files into R1 and R2 (forward reads, reverse reads)

## These samples are from Bioproject PRJNA438822. An experiment on Daphnia pulex, looking to see if maternal age has an effect on how offspring respond to 
##  selenium toxicity with respect to reproduction, lifespan, and resistance to heat-induced stress
## https://www.ncbi.nlm.nih.gov/bioproject?LinkName=biosample_bioproject&from_uid=8729509

vdb-config --interactive

##The three SRA sequence files shown below are meant to exemplify how to download given SRA sequences of interest. For our project, there were a total of 18 samples, not fully shown here for simplicity's sake.   

fastq-dump -F --split-files SRR6853331 #Sample 13; a pool a three individuals from Experimental Lab Population 2. Their mothers were young (8 days old) when they were born (2nd-3rd brood). The sample individuals were reared in 408 ug of Seleno-L-methionine / L of COMBO media.
fastq-dump -F --split-files SRR6853332 #Sample 14; a pool a three individuals from Experimental Lab Population 2. Their mothers were young (8 days old) when they were born (2nd-3rd brood). The sample individuals were reared in 102 ug of Seleno-L-methionine / L of COMBO media.
fastq-dump -F --split-files SRR6853336 #Sample 18; a pool a three individuals from Experimental Lab Population 1. Their mothers were young (8 days old) when they were born (2nd-3rd brood). The sample individuals were reared in 0 ug of Seleno-L-methionine / L of COMBO media.

##### Extra ####
## If you are downloading data from a sequencing company instead of NCBI, using wget for example, then calculate the md5sum values of all the files in the folder (./*), and read into a text file.
## then, you can compare the values in this file with the ones provided by the company.
## md5sum ./* > md5sum.txt

##### Extra ####
## If you data comes with multiple R1 and R2 files per individual, you can contatenate them together using "cat" before running FASTQC.
## See examples below for one file. You will probably want to use a loop to process through all of the files.
#cat SRR6819014*_R1_*.fastq.gz > SRR6819014_All_R1.fastq.gz
#cat SRR6819014*_R2_*.fastq.gz > SRR6819014_All_R2.fastq.gz


############## FASTQC to assess quality of the sequence data
## FastQC: run on each of the data files that have 'All' to check the quality of the data
## The output from this analysis is a folder of results and a zipped file of results
## NOTE: Make sure that your $CS directory has already been made before proceeding... or else make it now! 

fastqc *.fastq --outdir=$WD/$CS

#######  Tarball the directory containing the FASTQC results so we can easily bring it back to our computer to evaluate (your CS directory)
cd $WD/$CS
tar cvzf $CS.tar.gz $WD/$CS/*

## when finished use scp or rsync to bring the tarballed .gz results file to YOUR computer and open the .html file to evaluate the quality of your raw data.
