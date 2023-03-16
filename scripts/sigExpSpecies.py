from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='blast species with significant diff. expression levels')
parser.add_argument('-i', '--input', help='input fasta file of the CDS proteins', required=True)
parser.add_argument('-s', '--significant', help='input tsv with the significantly differentially expressed ids', required=True)
parser.add_argument('-o', '--output', help='input output path for differentially expressed proteins', required=True)

def openFasta(inputFasta):
    with open(inputFasta, 'r') as file:
        sequences = list(SeqIO.parse(file, 'fasta'))
    return sequences

def extractSigProteins(sequence, significantFile, output):
    df = pd.read_table(significantFile, sep='\t')
    significantIDs = list(df['target_id'])
    significantSeqRecord = [s for s in sequence if s.id == significantIDs[0]][0]
    with open(output, 'w') as file:
        SeqIO.write(significantSeqRecord, file, 'fasta')

if __name__ == '__main__':
    args = parser.parse_args()
    seqs = openFasta(args.input)
    sigSeqs = extractSigProteins(seqs, args.significant, args.output)
