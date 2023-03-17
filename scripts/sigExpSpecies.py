# imports
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
import pandas as pd

# command line arguments and flags; explained in the help string
parser = argparse.ArgumentParser(description='blast species with significant diff. expression levels')
parser.add_argument('-i', '--input', help='input fasta file of the CDS proteins', required=True)
parser.add_argument('-s', '--significant', help='input tsv with the significantly differentially expressed ids', required=True)
parser.add_argument('-o', '--output', help='input output path for differentially expressed proteins', required=True)

def openFasta(inputFasta):
    '''
    openning the fasta sequences of the reference protein fasta file downloaded earlier
    '''
    with open(inputFasta, 'r') as file:
        sequences = list(SeqIO.parse(file, 'fasta')) # storing fasta records into a list
    return sequences

def extractSigProteins(sequence, significantFile, output):
    '''
    take in all the sequences and extract the most significant record
    '''
    df = pd.read_table(significantFile, sep='\t') # openning the file of signficant differentially expressed genes coming from the sleuth statistical analysis
    significantIDs = list(df['target_id']) # extracting the target ids
    significantSeqRecord = [s for s in sequence if s.id == significantIDs[0]][0] # the records are organized by significance so the most significant one is the top records; thus look for the first significant id in the protein fasta reference to extract the record and save for later blast step
    with open(output, 'w') as file:
        SeqIO.write(significantSeqRecord, file, 'fasta') # storing for blast query later

# runs the functions above using commands from the command line when called from the command line
if __name__ == '__main__':
    args = parser.parse_args()
    seqs = openFasta(args.input)
    sigSeqs = extractSigProteins(seqs, args.significant, args.output)
