# imports
import os
import argparse

# inputs as flags and arguments from command line
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', required=True) # input file for fastq
parser.add_argument('-o', '--output', required=True) # output folder/file for fastq after decompressed from SRA

def decompress(i, o):
    '''
    decompress function allows for each SRA file into its independent fastq files reverse and forward reads
    '''
    inpath = i
    outpath = o

    rawFiles = list(os.listdir(inpath)) # list of the SRA file names

    for rawFile in rawFiles:
        os.system(f'fastq-dump -I --split-files {inpath}/{rawFile}') # call the fastq-dump function to split files into forward and reverse reads
    os.system(f'mv *.fastq {outpath}') # move all files into specified output folder

if __name__ == '__main__': # run from command line
    args = parser.parse_args() # create args objects
    decompress(args.input, args.output) # call on decompress function
