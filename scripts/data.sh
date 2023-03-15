mkdir ./PipelineProject_Rohan_Sethi
mkdir ./PipelineProject_Rohan_Sethi/data
mkdir ./PipelineProject_Rohan_Sethi/results
mkdir ./PipelineProject_Rohan_Sethi/data/raw
mkdir ./PipelineProject_Rohan_Sethi/data/fastq
mkdir ./PipelineProject_Rohan_Sethi/data/index
mkdir ./PipelineProject_Rohan_Sethi/data/blast

LOG=$1
touch $LOG

cd ./PipelineProject_Rohan_Sethi/data/raw

FILE=$2
append='../../../'
FILE=$append$FILE
while read line; do
	echo 'Downloading' $line'...'
	wget $line
done < $FILE

cd ../../..
