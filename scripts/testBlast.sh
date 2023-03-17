# input arguments from command line; query protein sequence, database to query against, and results output path
# this is for the test run so almost exact copy to the entire run one with some small change
QUERY=$1
DB=$2
OUTresults=$3

tblastn -query $QUERY -db $DB -max_hsps 1 -out $OUTresults -outfmt "6 sacc pident length qstart qend sstart send bitscore evalue stitle" # running the blast algorithm from the command line to find strains that may express this significant differentially expressed gene from samples
