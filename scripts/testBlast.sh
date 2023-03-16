QUERY=$1
DB=$2
OUTresults=$3

tblastn -query $QUERY -db $DB -out $OUTresults -outfmt "6 sacc pident length qstart qend sstart send bitscore evalue stitle"
