mkdir ./PipelineProject_Rohan_Sethi
mkdir ./PipelineProject_Rohan_Sethi/data
mkdir ./PipelineProject_Rohan_Sethi/results
mkdir ./PipelineProject_Rohan_Sethi/data/fastq
mkdir ./PipelineProject_Rohan_Sethi/data/index
mkdir ./PipelineProject_Rohan_Sethi/data/blast

LOG=$1
touch $LOG

cd ./PipelineProject_Rohan_Sethi/data/

FOLDER=$2
append='../../'
fqFOLDER=$append$FOLDER'/fastq'
blastdb=$append$FOLDER'/blastdb'

cp -r fqFOLDER ./fastq
cp blastdb/* ./blast

cd ../..
