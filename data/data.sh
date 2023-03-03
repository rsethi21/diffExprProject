mkdir ./data/raw
mkdir ./data/fastq
cd ./data/raw
echo 'Downloading SRR5660030'
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030
echo 'Downloading SRR5660033'
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033
echo 'Downloading SRR5660044'
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660044/SRR5660044
echo 'Downloading SRR5660045'
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660045/SRR5660045
cd ../..
