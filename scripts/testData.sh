# creating and organizing all the data output under the appropriate output folder
# created a separate instance for test data because allows for the minor changes to be addressed clearly; its close cousin the data.sh script is for the entire run from start to finish while this is for test run
mkdir ./PipelineProject_Rohan_Sethi
mkdir ./PipelineProject_Rohan_Sethi/data
mkdir ./PipelineProject_Rohan_Sethi/results
mkdir ./PipelineProject_Rohan_Sethi/data/fastq
mkdir ./PipelineProject_Rohan_Sethi/data/index
mkdir ./PipelineProject_Rohan_Sethi/data/blast

# input of log file to create it
LOG=$1
touch $LOG

# going into the directory to capture all test data in the folders created above to ensure coherence among paths for all scripts
cd ./PipelineProject_Rohan_Sethi/data/

FOLDER=$2 # input of folder name where test data is stored
append='../..'

# creating paths for each of the essential test data stored for test run 
fqFOLDER=$append$FOLDER'/fastq/*'

blastdb=$append$FOLDER'/blastdb/*'

meta=$append$FOLDER'/metatable.tsv'

# copying test data into the organized folder for coherence among all the scripts in test run
cp -r $fqFOLDER ./fastq
cp $blastdb ./blast
cp $meta .

cd ../.. # returning to project directory
