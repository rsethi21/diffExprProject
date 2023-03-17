# imports
from Bio import Entrez
from Bio import SeqIO
import argparse
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def retreiveRaw(email, accession, output, poutput):
    '''
    retreive the reference genome using biopython
    '''
    Entrez.email = email
    handle = Entrez.efetch(db='nucleotide', id=accession, rettype='gb') # fetch the genbank record of the reference
    record = SeqIO.read(handle, "gb") # read the byte file

    cds = []
    protein = []
    # loop through the features in the genbank file
    for feature in record.features:
        if feature.type == 'CDS': # look for the CDS feature
            s = SeqRecord(Seq(feature.extract(record.seq)), id=feature.qualifiers['protein_id'][0], description=feature.qualifiers['protein_id'][0]) # extract the CDS and store into a new SeqRecord object with proteinID and sequence
            cds.append(s) # add to cds list
            protein.append(SeqRecord(Seq(feature.qualifiers['translation'][0]), id=feature.qualifiers['protein_id'][0], description=feature.qualifiers['protein_id'][0])) # do the same thing but instead to grab the protein sequence; this is because we will need the protein sequence later and since this is a virus CDS it may have a different amino acid encoding for translation unknown

    # store the output protein and genome fastas in the output files specified
    with open(output, 'w') as file:
        SeqIO.write(cds, file, 'fasta')
    
    with open(poutput, 'w') as pfile:
        SeqIO.write(protein, pfile, 'fasta')
    
    # returning string for the log file
    return f'The HCMV genome ({accession}) has {len(cds)} CDS.\n'

# creating arguments and flags for the command line function
parser = argparse.ArgumentParser(description='grab raw CDS for index')
parser.add_argument('-e', '--email', required=True) # email to access the genome using biopython
parser.add_argument('-a', '--accession', required=True) # accession id for index genome
parser.add_argument('-o', '--output', required=True) # output genome fasta path
parser.add_argument('-l', '--log', required=True) # output log file path
parser.add_argument('-p', '--protein', required=True) # output protein fasta path

# below runs all the code above and stores the message in the log file as specified in the assignment
# runs from the command line
if __name__ == '__main__':
    args = parser.parse_args()
    message = retreiveRaw(args.email, args.accession, args.output, args.protein)
    with open(args.log, 'a') as logFile:
        logFile.write(message)
