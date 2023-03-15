import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser(description='TPM calculations for each sample')
parser.add_argument('-m', '--metatable', help='metatable of samples', required=True)
parser.add_argument('-r', '--resultsFolder', help='path to results folder', required = True)
parser.add_argument('-l', '--logfile', help='input logfile path', required=True)

def output(metaFile, resultsFolder):
    meta = pd.read_table(metaFile)
    outputDict = {}
    for sampleID in os.listdir(resultsFolder):
        sample = sampleID
        condition = int(meta[meta['sample'] == sampleID]['condition'])
        sampleAbundance = pd.read_table(f'{resultsFolder}/{sampleID}/abundance.tsv')
        tpms = sampleAbundance['tpm']
        outputDict[sampleID] = {'sample': sample, 'condition': condition, 'min_tpm': round(tpms.min(), 1), 'med_tpm': round(tpms.median(), 1), 'mean_tpm': round(tpms.mean(), 1), 'max_tpm': round(tpms.max(), 1)}
    return pd.DataFrame.from_dict(outputDict, orient='index')
if __name__ == '__main__':
    args = parser.parse_args()
    df = output(args.metatable, args.resultsFolder)
    with open(args.logfile, 'a') as file:
        file.write('\n')
    df.to_csv(args.logfile, sep='\t', header=True, index=False, mode='a', encoding='ascii')
    with open(args.logfile, 'a') as file:
        file.write('\n')
