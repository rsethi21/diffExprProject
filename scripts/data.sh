# this script is run toward the beginning to ensure that the output is organized and predictablly stored for access by other modules that work together in this script
# this one will run if running the pipeline from the beginning, otherwise its close cousin the testData.sh will do a similar function for the test run
mkdir ./PipelineProject_Rohan_Sethi
mkdir ./PipelineProject_Rohan_Sethi/data
mkdir ./PipelineProject_Rohan_Sethi/results
mkdir ./PipelineProject_Rohan_Sethi/data/raw
mkdir ./PipelineProject_Rohan_Sethi/data/fastq
mkdir ./PipelineProject_Rohan_Sethi/data/index
mkdir ./PipelineProject_Rohan_Sethi/data/blast

# input log file as argument in command line
LOG=$1
touch $LOG # creates the logfile specified

cd ./PipelineProject_Rohan_Sethi/data/raw # move into raw file

FILE=$2 # file input argument to access links to download SRAs from
append='../../../'
FILE=$append$FILE # creating path to links
# look below goes through each link in the file and downloads using wget
while read line; do
	echo 'Downloading' $line'...'
	wget $line
done < $FILE

cd ../../.. # returning back to home directory
