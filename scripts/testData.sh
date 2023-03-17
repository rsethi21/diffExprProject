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
fqFOLDER=$append$FOLDER'/fastq/*'
echo $fqFOLDER
blastdb=$append$FOLDER'/blastdb/*'
echo $blastdb
meta=$append$FOLDER'/metatable.tsv'
echo $meta

cp -r $fqFOLDER ./fastq
cp $blastdb ./blast
cp $meta .

cd ../..
