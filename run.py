import os

os.system('chmod +x ./data/data.sh')
os.system('./data/data.sh')
os.system('python3 decompressFastQ.py -i ./data/raw -o ./data/fastq')
