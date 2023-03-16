import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='select n of the top queries from blast results')
parser.add_argument('-i', '--input', help='results from blast search tsv file', required=True)
parser.add_argument('-n', '--number', help='number of queries to store', required=True)
parser.add_argument('-l', '--logfile', help='input path to logfile', required=True)

def selectIt(inputpath, number, logfile):
    columnHeaders = ['sacc', 'pident', 'length', 'qstart', 'qend', 'sstart', 'send', 'bitscore', 'evalue', 'stitle']
    df = pd.read_table(inputpath, sep='\t')
    df.columns = columnHeaders
    sortedDF = df.sort_values(by=['evalue'])
    finalDF = sortedDF[:10]
    with open(args.logfile, 'a') as file:
        file.write('\n')
    finalDF.to_csv(args.logfile, sep='\t', header=True, index=False, mode='a', encoding='ascii')


if __name__ == '__main__':
    args = parser.parse_args()
    selectIt(args.input, args.number, args.logfile)
