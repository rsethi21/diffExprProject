# diffExprProject

## Description
This is a project for Loyola University Chicago's computational biology class (COMP 383).
This project presents a pipeline that enables a user to determine the most differentially expressed genes between different samples of sequencing reads and search a genetic database for other species of a specified family that express those genes. This is a highly dynamic tool that enables a user to customize the pipeline to their needs.
The project requires the following requirements preinstalled:
1. **kallisto** - command line application to quantify transcripts per million (TPMs) from a set of reads mapped to a reference genome for assembly
2. R packages - **sleuth**, **dplyr** for statistical analysis of TPMs to find the most differentially expressed transcripts between the samples
3. Python packages - the following set-up section

## Set-up
1. First clone this repository locally with the following command line arguments:
  git clone [HTTPS]
2. Once cloned, move into the main project folder such that your current working directory is the diffExprProject folder
3. Install python requirements
  a. create a virtual environment, conda environment, or install locally; below is an example of using virtual environments
    python3 -m venv venv -> this create a virtual environment
    source venv/bin/activate -> this activates it, should see a (venv)
  b. install requirements from requirements.txt
    pip3 install -r requirements.txt

## Test Run
1. python3 run.py -e [EMAIL] -t ./testFolder --> make sure email is one that can be used to access NCBI sequences; ensure this is run from the main project directory

The test folder contains sample fastq files, a Betaherpesvirinae genome database folder to blast the most significantly expressed genes across, and a metatable that stores information about sample data. 
The links folder is are the SRA links from where the data comes from.

## More Complicated Run
