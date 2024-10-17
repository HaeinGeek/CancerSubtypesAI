import sys
import os
import pandas as pd

# 프로젝트 루트 디렉토리를 sys.path에 추가
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(project_root)

from sequence_processor import ProteinDatabase, add_mutated_sequences

def main():
    # 1. 데이터 불러오기 
    filepath = 'data\processed\mutant_seq_not_found.csv' 
    unique_mutations_df = pd.read_csv(filepath)

    # 2. ProteinDatabase 인스턴스 생성
    db_names = ['uniprot', 'ncbi', 'pdb', 'ensembl']
    protein_db = ProteinDatabase(db_names)

    # 3. 각 데이터베이스별로 변이 서열 생성
    for db_name in db_names:
        print(f"\nProcessing database: {db_name}")
        protein_dict = protein_db.get_protein_dict(db_name)
        unique_mutations_df = add_mutated_sequences(unique_mutations_df, protein_dict, db_name)

    # 결과 저장
    output_filename = f'data/processed/mutant_seq_test.csv'
    unique_mutations_df.to_csv(output_filename, index=False)
    print(f"Mutated sequences saved to {output_filename}")

if __name__ == "__main__":
    main()