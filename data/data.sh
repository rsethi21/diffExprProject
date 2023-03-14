# mkdir ./data/raw
# mkdir ./data/fastq
# mkdir ./data/index
# touch PipelineProject.log
# cd ./data/raw

FILE=$1

while read line; do
	echo 'Downloading' $line'...'
	wget $line
done < $FILE

# cd ../..
