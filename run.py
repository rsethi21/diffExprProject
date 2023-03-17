# imports
import os
import argparse
import pandas as pd


'''
This script essentially orchestrates all the other scripts required to perform the differential gene expression analysis. It is required to be run from the home directory of this project.
'''


# command line arguments defined
parser = argparse.ArgumentParser(description='differential expression pipeline') # defining the argument parser
parser.add_argument('-s', '--input', help='input file with NCBI links', required=False, default='testData/links/fileLinks.txt') # this is the input file with NCBI links if you want to download fastq files from scratch
parser.add_argument('-i', '--index', help='input accession id for index', required=False, default='NC_006273.2') # this is the id for the reference database which defaults to the one in the problem but can be substituted
parser.add_argument('-e', '--email', help='input email for NCBI access', required=True) # this is required and is the email for the ability to use Entrez
parser.add_argument('-m', '--metatable', help='metatable tab deliminated', required=False, default='./testData/metatable.tsv') # this is the metadata table that is not required but can be substituted for the default one
parser.add_argument('-l', '--logfile', help='name/path of log file', required=False, default='./PipelineProject_Rohan_Sethi/PipelineProject.log') # this is the name/path of where the log file will be created defaults into the folder
parser.add_argument('-n', '--name', help='name of reference to blast against', required=False, default='Betaherpesvirinae') # this is the name of the species that the sequences will be blasted against; defaults to the one on the problem but can be substituted
parser.add_argument('-u', '--numSelect', help='number of blast results to store', required=False, default=10) # this is the argument to determine the number of significant blast results to store in the final csv file; defaults to 10 like in the problem
parser.add_argument('-t', '--testData', help='input test data folder name', required=False, default=None) # this is the folder of the testdata to run the test; need to input this


if __name__ == '__main__': # this will run everything below only if run from the command line and not imported as a package

    args = parser.parse_args() # creating an instance of arguments
    
    if args.testData == None: # if no user input of the testData folder, then runs from beginning; not the test run that is in the else conditional

        os.system('chmod +x ./scripts/data.sh') # this is making the bash script executable from the command lined

        os.system(f'./scripts/data.sh {args.logfile} {args.input}') # this is a bash file that will create all the necessary output folders and download the sequences from the links provided

        os.system('python3 ./scripts/decompressFastQ.py -i ./PipelineProject_Rohan_Sethi/data/raw -o ./PipelineProject_Rohan_Sethi/data/fastq') # this is a python script that will run fastqdump on all the SRA files downladed

        os.system(f'python3 ./scripts/retreiveIndex.py -e {args.email} -a {args.index} -o ./PipelineProject_Rohan_Sethi/data/index/raw.fasta -p ./PipelineProject_Rohan_Sethi/data/index/protein.fasta -l {args.logfile}') # this will retreive the fasta files for the reference genome (protein and nucleotides because will use the most significantly expressed protein sequence for blast in the last step)

        os.system(f'time kallisto index -i ./PipelineProject_Rohan_Sethi/data/index/index.idx ./PipelineProject_Rohan_Sethi/data/index/raw.fasta') # this will created the index out of the nucleotide sequences using kallisto

        meta = pd.read_table(args.metatable) # this opens the metatable
        names = list(meta['sample']) # this captures all the SRA file names

        for fn in names: # for every SRA file
            f1 = f'./PipelineProject_Rohan_Sethi/data/fastq/{fn}_1.fastq' # save the path of the forward read
            f2 = f'./PipelineProject_Rohan_Sethi/data/fastq/{fn}_2.fastq' # save the path of the reverse read
            os.system(f'time kallisto quant -i ./PipelineProject_Rohan_Sethi/data/index/index.idx -o ./PipelineProject_Rohan_Sethi/results/{fn} -b 30 -t 4 {f1} {f2}') # quantify using kallisto the transcripts per million to determine the abundance of each gene expressed; the gene is determined by aligning to the index and the number of alignments = the transcripts which is normalized to get TPMs

        os.system(f'python3 ./scripts/quantify.py -m {args.metatable} -r ./PipelineProject_Rohan_Sethi/results -l {args.logfile}') # here I save the summary of the TPMS per sample and store it in the logfile

        os.system(f'Rscript ./scripts/diffExpAnalysis.R {args.metatable} {args.logfile} ./PipelineProject_Rohan_Sethi/results/sigDiffExp.tsv') # here I run an R script I wrote to find the most significantly differentially expressed TPMS between samples; the ids are stored in a file that will be accessed by other scripts

        os.system(f'python3 ./scripts/sigExpSpecies.py -i ./PipelineProject_Rohan_Sethi/data/index/protein.fasta -s ./PipelineProject_Rohan_Sethi/results/sigDiffExp.tsv -o ./PipelineProject_Rohan_Sethi/data/blast/mostDifferentiallyExpressed.fasta') # here I find the single most signficant differentially expressed gene from the step above and I store its protein sequence by extracting it from the protein fasta I made of the reference using the protein id

        os.system('rm ./PipelineProject_Rohan_Sethi/results/sigDiffExp.tsv') # here I remove an intermediate file that is not necessary anymore

        os.system(f'python3 ./scripts/blastReference.py -n {args.name} -e {args.email} -o ./PipelineProject_Rohan_Sethi/data/blast/{args.name}.fasta') # here I extract the genome files of the Betaherpesvirinae from NCBI; this can be changed to a different one but by default uses Betaherpesvirinae

        os.system('chmod +x ./scripts/blast.sh') # here I make this shell script executable

        os.system(f'./scripts/blast.sh {args.name} ./PipelineProject_Rohan_Sethi/data/blast/{args.name}.fasta ./PipelineProject_Rohan_Sethi/data/blast/mostDifferentiallyExpressed.fasta ./PipelineProject_Rohan_Sethi/data/blast/{args.name} ./PipelineProject_Rohan_Sethi/results/blastResults.csv') # here I run the blast tblastn, blasting the most significant protein against the Betavirinae genome sequences

        os.system(f'python3 ./scripts/selectN.py -i ./PipelineProject_Rohan_Sethi/results/blastResults.csv -n {args.numSelect} -l {args.logfile}') # here I select the top 10 matches from the blast search of Betavirinae strains found for the protein in the blast serach; only the top 10 most significant mathces (HPSs) are selected from as the tblastn command in the blast.sh and testBlast.sh file deal with this using the flag max_hps=1


    else: # this will run if using the testData for the test run; it uses the same scripts with some exceptions so the comments are the same as above and mentioned as such
        os.system('chmod +x ./scripts/testData.sh') # I wrote a separate shell file that doesn't require downloading fastq files since samples will be provided; this script will create the necessary folders for the testrun and extract the data needed from the testData folder provided

        os.system(f'./scripts/testData.sh {args.logfile} {args.testData}') # this will run the script explained above which is made executable above

        os.system(f'python3 ./scripts/retreiveIndex.py -e {args.email} -a {args.index} -o ./PipelineProject_Rohan_Sethi/data/index/raw.fasta -p ./PipelineProject_Rohan_Sethi/data/index/protein.fasta -l {args.logfile}')# here we retreive the index like above

        os.system(f'time kallisto index -i ./PipelineProject_Rohan_Sethi/data/index/index.idx ./PipelineProject_Rohan_Sethi/data/index/raw.fasta') # here we create the index like above
        
        # the steps below run kallisto quanitification for each of the samples like mentioned above

        meta = pd.read_table(args.metatable)
        names = list(meta['sample'])

        for fn in names:
            f1 = f'./PipelineProject_Rohan_Sethi/data/fastq/{fn}_1.fastq'
            f2 = f'./PipelineProject_Rohan_Sethi/data/fastq/{fn}_2.fastq'
            os.system(f'time kallisto quant -i ./PipelineProject_Rohan_Sethi/data/index/index.idx -o ./PipelineProject_Rohan_Sethi/results/{fn} -b 30 -t 4 {f1} {f2}')

    # the three scripts below serve the same function they did above

        os.system(f'python3 ./scripts/quantify.py -m ./PipelineProject_Rohan_Sethi/data/metatable.tsv -r ./PipelineProject_Rohan_Sethi/results -l {args.logfile}')

        os.system(f'Rscript ./scripts/diffExpAnalysis.R ./PipelineProject_Rohan_Sethi/data/metatable.tsv {args.logfile} ./PipelineProject_Rohan_Sethi/results/sigDiffExp.tsv')
       
        os.system(f'python3 ./scripts/sigExpSpecies.py -i ./PipelineProject_Rohan_Sethi/data/index/protein.fasta -s ./PipelineProject_Rohan_Sethi/results/sigDiffExp.tsv -o ./PipelineProject_Rohan_Sethi/data/blast/mostDifferentiallyExpressed.fasta') # here I find the most signficant differentially expressed gene from the step above and I store its protein sequence by extracting it  from the protein fasta I made of the reference using the protein id

    # the scripts below run the blast search using the predownloaded blast db of Betaherpesvirinae in the testData folder, which is why I use a slightly different script for the blasting of test data

        os.system('chmod +x ./scripts/testBlast.sh')

        os.system(f'./scripts/testBlast.sh ./PipelineProject_Rohan_Sethi/data/blast/mostDifferentiallyExpressed.fasta ./PipelineProject_Rohan_Sethi/data/blast/{args.name} ./PipelineProject_Rohan_Sethi/results/blastResults.csv')

   # this serves the same purpose as it did above

        os.system(f'python3 ./scripts/selectN.py -i ./PipelineProject_Rohan_Sethi/results/blastResults.csv -n {args.numSelect} -l {args.logfile}')
