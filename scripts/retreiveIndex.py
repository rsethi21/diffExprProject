from Bio import Entrez
from Bio import SeqIO
import argparse
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def retreiveRaw(email, accession, output):
    Entrez.email = email
    handle = Entrez.efetch(db='nucleotide', id=accession, rettype='gb')
    record = SeqIO.read(handle, "gb")

    cds = []
    for feature in record.features:
        if feature.type == 'CDS':
            s = SeqRecord(Seq(feature.extract(record.seq)), id=feature.qualifiers['protein_id'][0], description=feature.qualifiers['protein_id'][0])
            cds.append(s)

    with open(output, 'w') as file:
        SeqIO.write(cds, file, 'fasta')
    
    return f'The HCMV genome ({accession}) has {len(cds)} CDS.\n'

parser = argparse.ArgumentParser(description='grab raw CDS for index')
parser.add_argument('-e', '--email', required=True)
parser.add_argument('-a', '--accession', required=True)
parser.add_argument('-o', '--output', required=True)
parser.add_argument('-l', '--log', required=True)

if __name__ == '__main__':
    args = parser.parse_args()
    message = retreiveRaw(args.email, args.accession, args.output)
    with open(args.log, 'a') as logFile:
        logFile.write(message)
