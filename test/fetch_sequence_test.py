import sys
import os

# 프로젝트 루트 디렉토리를 sys.path에 추가
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(project_root)

from sequence_fetcher.sequence_fetcher import SequenceFetcher
import pandas as pd

def main():
    # 각 데이터베이스별 유전자 리스트 불러오기
    db_names = ['uniprot', 'ncbi', 'pdb', 'ensembl']
    genes_dict = {}

    # 데이터베이스별 파일을 로드하여 딕셔너리에 저장
    for db_name in db_names:
        filename = f'data/genes_{db_name}.txt'
        with open(filename, 'r') as f:
            genes_dict[db_name] = [line.strip() for line in f]

    email = 'your_email@example.com'
    fetcher = SequenceFetcher()

    # 각 데이터베이스별 isoform 처리
    fetcher.process_isoforms('uniprot', genes_dict['uniprot'][:10])
    fetcher.process_isoforms('ncbi', genes_dict['ncbi'][:10], email=email)
    fetcher.process_isoforms('pdb', genes_dict['pdb'][:10])
    fetcher.process_isoforms('ensembl', genes_dict['ensembl'][:10], species='homo_sapiens', max_workers=5)

if __name__ == "__main__":
    main()