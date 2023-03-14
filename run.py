import os
import argparse
import pdb

parser = argparse.ArgumentParser(description='differential expression pipeline')
parser.add_argument('-s', '--input', help='input file with NCBI links', required=False, default='./data/links/testInput.txt')
parser.add_argument('-i', '--index', help='input accession id for index', required=False, default='NC_006273.2')
parser.add_argument('-e', '--email', help='input email for NCBI access', required=True)

args = parser.parse_args()

# os.system('chmod +x ./data/data.sh')
# os.system('./data/data.sh')
# os.system('python3 ./data/decompressFastQ.py -i ./data/raw -o ./data/fastq')
os.system(f'python3 ./data/retreiveIndex.py -e {args.email} -a {args.index} -o ./data/index/raw.fasta')
