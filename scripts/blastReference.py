# imports
from Bio import Entrez
from Bio import SeqIO
import argparse

# below is creating an object for command line arguments
parser = argparse.ArgumentParser(description='extract local database for species')
parser.add_argument('-n', '--name', help='name of species to blast against', required=True)
parser.add_argument('-e', '--email', help='email to extract using', required=True)
parser.add_argument('-o', '--output', help='path to save sequences for local database', required=True)

def extract(name, email, output):
    '''
    This function extracts the genome entries of the species one would like to blast the query against, allowing for a maximum of 20,000 entries
    '''
    # the four lines below is utilizing biopython's esearch functionality
    Entrez.email = email # email to search
    handle = Entrez.esearch(db='nucleotide', term=f'{name}[ORGN]', idtype="acc", retmax=10000) # retreiving the max number of entries from the NCBI nucleotides using the term organism; name from argument
    record = Entrez.read(handle) # reading the byte file
    records = list(record['IdList']) # creating a list of the accession ids for each of the entries

    # below is extra code to grab the next 10,000 accession ids starting from the 10000th entries; in try except in case not enough entries to start from 10000th
    try:
        handleExtra = Entrez.esearch(db='nucleotide', term=f'{name}[ORGN]', idtype="acc", retstart=10000, retmax=10000)
        recordExtra = Entrez.read(handleExtra)
        records.extend(list(recordExtra['IdList']))
        print(len(records))
    except:
        pass
    
    # going through each of the accession ids and fetching the fasta entries and putting into one file
    for ident in records:
        handle2 = Entrez.efetch(db='nucleotide', id=ident, rettype='fasta')
        records2 = SeqIO.read(handle2, 'fasta')
        with open(output, 'a') as file: # appending all the entries into the fasta file specified in the argument
            SeqIO.write(records2, file, 'fasta')

if __name__ == '__main__': # this if statement ensures that whatever in if statement run only when called from the command line
    args = parser.parse_args() # creating argument objects
    extract(args.name, args.email, args.output) # calling the extract function above
