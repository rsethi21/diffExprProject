# imports
import pandas as pd
import argparse
import os

# arguments from the command line; information about each explained in the help
parser = argparse.ArgumentParser(description='TPM calculations for each sample')
parser.add_argument('-m', '--metatable', help='metatable of samples', required=True)
parser.add_argument('-r', '--resultsFolder', help='path to results folder', required = True)
parser.add_argument('-l', '--logfile', help='input logfile path', required=True)

def output(metaFile, resultsFolder):
    '''
    take the TPM output from the quant step of kallisto and calculate the minimum, maximum, mean, median TPM counts from each sample
    '''
    meta = pd.read_table(metaFile)
    outputDict = {}
    for sampleID in os.listdir(resultsFolder): # for each sample take the output abundances and TPMs
        sample = sampleID # extract id of sample
        condition = int(meta[meta['sample'] == sampleID]['condition']) # extract the condition from meta table
        sampleAbundance = pd.read_table(f'{resultsFolder}/{sampleID}/abundance.tsv') # read the abundance TPM values; using the sampleID to extract the appropriate abundance values
        tpms = sampleAbundance['tpm'] # extract the tpms from this folder
        outputDict[sampleID] = {'sample': sample, 'condition': condition, 'min_tpm': round(tpms.min(), 1), 'med_tpm': round(tpms.median(), 1), 'mean_tpm': round(tpms.mean(), 1), 'max_tpm': round(tpms.max(), 1)}
    return pd.DataFrame.from_dict(outputDict, orient='index') # returning the table for storage into log file

# run the following code if run from the command line
if __name__ == '__main__':
    args = parser.parse_args()
    df = output(args.metatable, args.resultsFolder) # call on the function above
    # open log file and append the table output into the log file and add spaces in between for readability
    with open(args.logfile, 'a') as file:
        file.write('\n')
    df.to_csv(args.logfile, sep='\t', header=True, index=False, mode='a', encoding='ascii') 
    with open(args.logfile, 'a') as file:
        file.write('\n')
