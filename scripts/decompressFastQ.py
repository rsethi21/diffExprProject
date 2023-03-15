import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', required=True)
parser.add_argument('-o', '--output', required=True)

args = parser.parse_args()

def decompress(i, o):
    inpath = args.input
    outpath = args.output

    rawFiles = list(os.listdir(inpath))

    for rawFile in rawFiles:
        print(rawFile)
        os.system(f'fastq-dump -I --split-files {inpath}/{rawFile}')
    os.system(f'mv *.fastq {outpath}')

if __name__ == '__main__':
    decompress(args.input, args.output)
