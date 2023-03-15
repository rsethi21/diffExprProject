from Bio import Entrez
from Bio import SeqIO
import argparse
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def retreiveRaw(email, accession, output, poutput):
    Entrez.email = email
    handle = Entrez.efetch(db='nucleotide', id=accession, rettype='gb')
    record = SeqIO.read(handle, "gb")

    cds = []
    protein = []
    for feature in record.features:
        if feature.type == 'CDS':
            s = SeqRecord(Seq(feature.extract(record.seq)), id=feature.qualifiers['protein_id'][0], description=feature.qualifiers['protein_id'][0])
            cds.append(s)
            protein.append(SeqRecord(Seq(feature.qualifiers['translation'][0]), id=feature.qualifiers['protein_id'][0], description=feature.qualifiers['protein_id'][0]))

    with open(output, 'w') as file:
        SeqIO.write(cds, file, 'fasta')
    
    with open(poutput, 'w') as pfile:
        SeqIO.write(protein, pfile, 'fasta')
    
    return f'The HCMV genome ({accession}) has {len(cds)} CDS.\n'

parser = argparse.ArgumentParser(description='grab raw CDS for index')
parser.add_argument('-e', '--email', required=True)
parser.add_argument('-a', '--accession', required=True)
parser.add_argument('-o', '--output', required=True)
parser.add_argument('-l', '--log', required=True)
parser.add_argument('-p', '--protein', required=True)

if __name__ == '__main__':
    args = parser.parse_args()
    message = retreiveRaw(args.email, args.accession, args.output, args.protein)
    with open(args.log, 'a') as logFile:
        logFile.write(message)
