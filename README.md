# diffExprProject

## Description
- This is a project for Loyola University Chicago's computational biology class (COMP 383).
- This project presents a pipeline that enables a user to determine the most differentially expressed genes between different samples of sequencing reads and search a genetic database for other species of a specified family that express those genes. This is a highly dynamic tool that enables a user to customize the pipeline to their needs.
- The project requires the following requirements preinstalled:
1. **kallisto** - command line application to quantify transcripts per million (TPMs) from a set of reads mapped to a reference genome for assembly
2. **ncbi-blast+** - command line application to run a blast search of a query against a genome database locally downloaded
3. R packages - **sleuth**, **dplyr** for statistical analysis of TPMs to find the most differentially expressed transcripts between the samples
4. Python packages - the following set-up section

## Set-up
1. First clone this repository locally with the following command line arguments:
    - git clone [HTTPS]
2. Once cloned, move into the main project folder such that your current working directory is the diffExprProject folder
3. Install python requirements
    - create a virtual environment, conda environment, or install locally; below is an example of using virtual environments
        - python3 -m venv venv -> this create a virtual environment
        - source venv/bin/activate -> this activates it, should see a (venv)
    - install requirements from requirements.txt
        - pip3 install -r requirements.txt

## Test Run
### RUN
1. python3 run.py -e [EMAIL] -t ./testFolder
- make sure email is one that can be used to access NCBI sequences; ensure this is run from the main project directory
- The test folder contains sample fastq files, a Betaherpesvirinae genome database folder to blast the most significantly expressed genes across, and a metatable that stores information about sample data. 
- The links folder is are the SRA links from where the data comes from.
- For more information on the flags and arguments:
  - python3 run.py --help
### Output
- The outputs will be saved in PipelineProject_Rohan_Sethi directory in which results are saved in results and data for input into various functions is saved in data
- The log file contains information from the test run

## More Complicated Run
- all flags in brackets are optional and default to the follwoing:
    - -s = testData/links/fileLinks.txt
    - -i = NC_006273.2
    - -e = no default, need to specify one
    - -m = testData/metatable.tsv
    - -l = ./PipelineProject_Rohan_Sethi/PipelineProject.log
    - -n = Betaherpesvirinae
    - -u = 10
    - -t = None

- usage: run.py [-h] [-s INPUT] [-i INDEX] -e EMAIL [-m METATABLE] [-l LOGFILE] [-n NAME] [-u NUMSELECT] [-t TESTDATA]

- options:
  - -h, --help            show this help message and exit
  - -s INPUT, --input INPUT
                        input file with NCBI links for the SRA sample data to download from
  - -i INDEX, --index INDEX
                        input accession id for index to assemble the reads
  - -e EMAIL, --email EMAIL
                        input email for NCBI access via biopython
  - -m METATABLE, --metatable METATABLE
                        metatable tab deliminated containing information about the samples; an example in testData
  - -l LOGFILE, --logfile LOGFILE
                        name/path of log file to store important output information tab delimited
  - -n NAME, --name NAME  
                        name of species to blast against to see what other species the most differentially expressed genes are expressed in
  - -u NUMSELECT, --numSelect NUMSELECT
                        number of blast results to store from the blast search
  - -t TESTDATA, --testData TESTDATA
                        input test data folder name; only if wanting to run test run, else ignore
