from Bio import Entrez
from Bio import SeqIO
import argparse
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

parser = argparse.ArgumentParser(description='Retreive fasta of CDSs from index')
parser.add_argument('-e', '--email', help='input email address to access sequences with', required=True)
parser.add_argument('-i', '--input', help='input accession', required=True)
parser.add_argument('-o', '--output', help='specify output file path', required=True)

args = parser.parse_args()

Entrez.email = args.email
handle = Entrez.efetch(db='nucleotide', id=args.input, rettype='gb')
record = SeqIO.read(handle, "gb")

cds = []
for feature in record.features:
    if feature.type == 'CDS':
        s = SeqRecord(Seq(feature.extract(record.seq)), id=feature.qualifiers['protein_id'][0], description=feature.qualifiers['protein_id'][0])
        cds.append(s)

with open(args.output, 'w') as file:
    SeqIO.write(cds, file, 'fasta')
