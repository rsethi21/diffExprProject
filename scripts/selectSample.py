from Bio import SeqIO
import argparse
import pandas as pd
from tqdm import tqdm

parser = argparse.ArgumentParser(description='select a subset of fastq records')
parser.add_argument('-i', '--input', help='input folder for reads', required=True)
parser.add_argument('-a', '--amount', help='amount of fastq records to keep', required=True)
parser.add_argument('-o', '--output', help='output folder', required=False, default='.')
parser.add_argument('-m', '--metafile', help='input metafile path', required=True)

def extractSequences(fq1, fq2):
    with open(fq1, 'r') as f1:
        fastq1 = list(SeqIO.parse(f1, 'fastq'))

    with open(fq2, 'r') as f2:
        fastq2 = list(SeqIO.parse(f2, 'fastq'))
    
    return fastq1, fastq2

def takeSample(seqs1, seqs2, amount, out1, out2):
    with open(out1, 'w') as o1:
        SeqIO.write(seqs1[:int(amount)], o1, 'fastq')

    with open(out2, 'w') as o2:
        SeqIO.write(seqs2[:int(amount)], o2, 'fastq')


if __name__ == '__main__':
    args = parser.parse_args()
    m = pd.read_table(args.metafile)
    names = m['sample']
    for n in tqdm(names):
        forward = f'{args.input}/{n}_1.fastq'
        reverse = f'{args.input}/{n}_2.fastq'
        fq1, fq2 = extractSequences(forward, reverse)
        forwardOut = f'{args.output}/{n}_1.fastq'
        reverseOut = f'{args.output}/{n}_2.fastq'
        takeSample(fq1, fq2, args.amount, forwardOut, reverseOut)
