from Bio import Entrez
from Bio import SeqIO
import argparse
import pdb

parser = argparse.ArgumentParser(description='extract local database for species')
parser.add_argument('-n', '--name', help='name of species to blast against', required=True)
parser.add_argument('-e', '--email', help='email to extract using', required=True)
parser.add_argument('-o', '--output', help='path to save sequences for local database', required=True)

def extract(name, email, output):
    Entrez.email = email
    handle = Entrez.efetch(db='nucleotide', id=name, retmode='fasta')
    records = Entrez.read(handle)
    print(len(records))
    with open(output, 'w') as file:
        SeqIO.write(records, file, 'fasta')

if __name__ == '__main__':
    args = parser.parse_args()
    extract(args.name, args.email, args.output)
