mkdir ./data
mkdir ./results
mkdir ./data/raw
mkdir ./data/fastq
mkdir ./data/index

LOG=$1
touch $LOG

cd ./data/raw

FILE=$2
FILE+="../../"
while read line; do
	echo 'Downloading' $line'...'
	wget $line
done < $FILE

cd ../..
