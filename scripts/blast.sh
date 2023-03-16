NAME=$1
REF=$2
QUERY=$3
OUTdb=$4
OUTresults=$5

makeblastdb -in $REF -out $OUTdb -title $NAME -dbtype nucl

tblastn -query $QUERY -db $OUTdb -out $OUTresults -outfmt "6 sacc pident length qstart qend sstart send bitscore evalue stitle"
