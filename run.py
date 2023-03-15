import os
import argparse

parser = argparse.ArgumentParser(description='differential expression pipeline')
parser.add_argument('-s', '--input', help='input file with NCBI links', required=False, default='./data/links/testInput.txt')
parser.add_argument('-i', '--index', help='input accession id for index', required=False, default='NC_006273.2')
parser.add_argument('-e', '--email', help='input email for NCBI access', required=True)

args = parser.parse_args()

# os.system('chmod +x ./data/data.sh')
# os.system(f'./data/data.sh {args.input}')
# os.system('python3 ./data/decompressFastQ.py -i ./data/raw -o ./data/fastq')
# os.system(f'python3 ./data/retreiveIndex.py -e {args.email} -a {args.index} -o ./data/index/raw.fasta -l ./PipelineProject.log')
# os.system(f'time kallisto index -i ./data/index/index.idx ./data/index/raw.fasta')
for fn in os.listdir('./data/raw'):
    f1 = f'./data/fastq/{fn}_1.fastq'
    f2 = f'./data/fastq/{fn}_2.fastq'
    os.system(f'time kallisto quant -i ./data/index/index.idx -o ./results/{fn} -b 30 -t 4 {f1} {f2}')
