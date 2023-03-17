# below are arguments for this bash file
NAME=$1 # name of output database
REF=$2 # path to fasta file of genomic data for species
QUERY=$3 # fasta of query sequence (most differentially expressed)
OUTdb=$4 # output datat folder for database
OUTresults=$5 # output path of the blast search

makeblastdb -in $REF -out $OUTdb -title $NAME -dbtype nucl # this creates a local database of specified organism/species

# below runs the blast query; in this case I choose to do a tblastn because we are querying a translated protein against a genome of nucleotides for a species
tblastn -query $QUERY -db $OUTdb -max_hsps 1 -out $OUTresults -outfmt "6 sacc pident length qstart qend sstart send bitscore evalue stitle" # the output is formated with specified fields and is tsv (option 6); max-hsps ensures that only the best query, strain alignment is chosen
