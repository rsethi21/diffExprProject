import os
import argparse

parser = argparse.ArgumentParser(description='differential expression pipeline')
parser.add_argument('-s', '--input', help='input file with NCBI links', required=False, default='testData/links/fileLinks.txt')
parser.add_argument('-i', '--index', help='input accession id for index', required=False, default='NC_006273.2')
parser.add_argument('-e', '--email', help='input email for NCBI access', required=True)
parser.add_argument('-m', '--metatable', help='metatable tab deliminated', required=False, default='./testData/metatable.tsv')
parser.add_argument('-l', '--logfile', help='name/path of log file', required=False, default='./PipelineProject_Rohan_Sethi/PipelineProject.log')
parser.add_argument('-n', '--name', help='name of reference to blast against', required=False, default='Betaherpesvirinae')


if __name__ == '__main__':

    args = parser.parse_args()
    
    os.system('chmod +x ./scripts/data.sh')

    os.system(f'./scripts/data.sh {args.logfile} {args.input}')

    os.system('python3 ./scripts/decompressFastQ.py -i ./PipelineProject_Rohan_Sethi/data/raw -o ./PipelineProject_Rohan_Sethi/data/fastq')

    os.system(f'python3 ./scripts/retreiveIndex.py -e {args.email} -a {args.index} -o ./PipelineProject_Rohan_Sethi/data/index/raw.fasta -l {args.logfile}')

    os.system(f'time kallisto index -i ./PipelineProject_Rohan_Sethi/data/index/index.idx ./PipelineProject_Rohan_Sethi/data/index/raw.fasta')

    for fn in os.listdir('./data/raw'):
        f1 = f'./PipelineProject_Rohan_Sethi/data/fastq/{fn}_1.fastq'
        f2 = f'./PipelineProject_Rohan_Sethi/data/fastq/{fn}_2.fastq'
        os.system(f'time kallisto quant -i ./PipelineProject_Rohan_Sethi/data/index/index.idx -o ./PipelineProject_Rohan_Sethi/results/{fn} -b 30 -t 4 {f1} {f2}')

    os.system(f'python3 ./scripts/quantify.py -m {args.metatable} -r ./PipelineProject_Rohan_Sethi/results -l {args.logfile}')

    os.system(f'Rscript ./scripts/diffExpAnalysis.R {args.metatable} {args.logfile} ./PipelineProject_Rohan_Sethi/results/sigDiffExp.tsv')

    # step 5
    os.system(f'python3 ./scripts/sigExpSpecies.py -i ./PipelineProject_Rohan_Sethi/data/index/raw.fasta -s ./PipelineProject_Rohan_Sethi/results/sigDiffExp.tsv -o ./PipelineProject_Rohan_Sethi/results/mostDifferentiallyExpressed.fasta')
    os.system('rm ./PipelineProject_Rohan_Sethi/results/sigDiffExp.tsv')
    '''
    os.system(f'python3 ./scripts/blastReference.py -n {args.name} -e {args.email} -o ./PipelineProject_Rohan_Sethi/data/blast/betaherpesvirinae.fasta')
    '''
