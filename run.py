import os

os.system('chmod +x ./data/data.sh')
os.system('./data/data.sh')
os.system('python3 ./data/decompressFastQ.py -i ./data/raw -o ./data/fastq')
