QUERY=$1
DB=$2
OUTresults=$3

echo $QUERY
echo $DB
echo $OUTresults

tblastn -query $QUERY -db $DB -max_hsps 1 -out $OUTresults -outfmt "6 sacc pident length qstart qend sstart send bitscore evalue stitle"
