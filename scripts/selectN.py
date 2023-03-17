# imports
import pandas as pd
import argparse

# command line arguments and flags; description and help explains use of each flag
parser = argparse.ArgumentParser(description='select n of the top queries from blast results')
parser.add_argument('-i', '--input', help='results from blast search tsv file', required=True)
parser.add_argument('-n', '--number', help='number of queries to store', required=True)
parser.add_argument('-l', '--logfile', help='input path to logfile', required=True)

def selectIt(inputpath, number, logfile):
    '''
    selecting the top n matches/alignments from the blast output
    '''
    columnHeaders = ['sacc', 'pident', 'length', 'qstart', 'qend', 'sstart', 'send', 'bitscore', 'evalue', 'stitle'] # adding these column headers to the blast output tsv
    df = pd.read_table(inputpath, sep='\t') # opening the blast output
    df.columns = columnHeaders # creating header
    sortedDF = df.sort_values(by=['evalue']) # sorting the dataframe from lowest e-score (best alignment) to highest (worst alignment)
    finalDF = sortedDF[:10] # selecting the top 10 from storted (the part where it must be only the best HSP per query-subject pair is handled by max-hsps=1 flag in the blast command run in the command line by the blast.sh file)
    # storing the top 10 in the log file
    with open(args.logfile, 'a') as file:
        file.write('\n')
    finalDF.to_csv(args.logfile, sep='\t', header=True, index=False, mode='a', encoding='ascii')

# run the above in the command line
if __name__ == '__main__':
    args = parser.parse_args()
    selectIt(args.input, args.number, args.logfile)
