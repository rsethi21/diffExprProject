# imports
from Bio import SeqIO
import argparse
import pandas as pd
from tqdm import tqdm

# arguments and flags whcih are explained further in the help and description
parser = argparse.ArgumentParser(description='select a subset of fastq records')
parser.add_argument('-i', '--input', help='input folder for reads', required=True)
parser.add_argument('-a', '--amount', help='amount of fastq records to keep', required=True)
parser.add_argument('-o', '--output', help='output folder', required=False, default='.')
parser.add_argument('-m', '--metafile', help='input metafile path', required=True)

def extractSequences(fq1, fq2):
    '''
    given a set a paired end reads, extracting fastq records in a list
    '''
    with open(fq1, 'r') as f1: # openning the forward read
        fastq1 = list(SeqIO.parse(f1, 'fastq'))

    with open(fq2, 'r') as f2: # opening the reverse read
        fastq2 = list(SeqIO.parse(f2, 'fastq'))
    
    return fastq1, fastq2

def takeSample(seqs1, seqs2, amount, out1, out2):
    '''
    given a set of paired end reads, extracting subset of fastq records and storing in separate test data folder
    '''
    with open(out1, 'w') as o1:
        SeqIO.write(seqs1[:int(amount)], o1, 'fastq') # saves the sequences from forward reads by splicing the list to grab only the amount specified

    with open(out2, 'w') as o2:
        SeqIO.write(seqs2[:int(amount)], o2, 'fastq') # saves the sequences from reverse reads by splicing the list to grab only the amount specified


if __name__ == '__main__':
    # runs all in the command line
    args = parser.parse_args() # creating command line arguments instance
    m = pd.read_table(args.metafile) # reading in the metadata table
    names = m['sample'] # grabbing the sample names from this table
    for n in tqdm(names): # for each sample, open the assosicated fastq files already fastq-dumped (tqdm creates a progress bar to see how far in the loop from the command line)
        forward = f'{args.input}/{n}_1.fastq' # forward read path _1
        reverse = f'{args.input}/{n}_2.fastq' # reverse read path _2
        fq1, fq2 = extractSequences(forward, reverse) # calling extract function
        forwardOut = f'{args.output}/{n}_1.fastq' # output path forward reads subset
        reverseOut = f'{args.output}/{n}_2.fastq' # output path reverse reads subset
        takeSample(fq1, fq2, args.amount, forwardOut, reverseOut) # calling on extractng sample function above
