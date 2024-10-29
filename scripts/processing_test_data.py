import sys
import os
import pandas as pd
import numpy as np

# 프로젝트 루트 디렉토리를 sys.path에 추가
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(project_root)

from sequence_processor import ProteinDatabase, add_mutated_sequences
from encoder.mutation_encoder import (
    classify_and_aggregate_mutations, 
    parse_multiple_mutations,
    process_mutation_features
)

##### generate mutants #####
# 1. 데이터 불러오기 
test_filepath = 'data/test.csv' 
amino_acid_filepath = 'data/amino_acid_features.csv'

test_df = pd.read_csv(test_filepath)
test_df = test_df.astype('str')
amino_acid_features = pd.read_csv(amino_acid_filepath).set_index('amino acid')

# 2. 유일한 유전자와 돌연변이 조합 추출
melted_df = test_df.melt(id_vars=['ID'], var_name='gene', value_name='mutation_str')
unique_mutations_df = melted_df[['gene', 'mutation_str']].drop_duplicates(subset = ['gene', 'mutation_str'])

# 3. ProteinDatabase 인스턴스 생성
db_names = ['uniprot', 'ncbi', 'pdb', 'ensembl']
protein_db = ProteinDatabase(db_names)

# 4. 각 데이터베이스별로 변이 서열 생성
for db_name in db_names:
    print(f"\nProcessing database: {db_name}")
    protein_dict = protein_db.get_protein_dict(db_name)
    unique_mutations_df = add_mutated_sequences(unique_mutations_df, protein_dict, db_name)

# 결과 저장
output_filename = f'data/processed/test/mutant_seq_unique.csv'
unique_mutations_df.to_csv(output_filename, index=False)
print(f"Mutated sequences saved to {output_filename}")

notfound_filename = f'data/processed/test/mutant_seq_not_found.csv'
mutant_seq_not_found = unique_mutations_df.loc[unique_mutations_df.wt_seq.isna(),['gene', 'mutation_str']]
mutant_seq_not_found.to_csv(notfound_filename, index=False)
print(f"Sequences not found saved to {notfound_filename}")

##### encode mutation #####
# 변이 인코딩 데이터프레임 생성
mutation_encoding_df = unique_mutations_df.copy()

# 변이 종류 분류
mutation_encoding_df['type'] = unique_mutations_df['mutation_str'].apply(classify_and_aggregate_mutations)

# 변이 아미노산 및 변이 위치 처리
mutation_encoding_df['origin'], mutation_encoding_df['position'], mutation_encoding_df['mutant'] = zip(
    *unique_mutations_df['mutation_str'].apply(parse_multiple_mutations)
    )

# 변이 위치의 개수 계산
mutation_encoding_df['mut_num'] = mutation_encoding_df['position'].apply(lambda x: len(x) if isinstance(x, list) else 0)


# 아미노산 특성 차이 계산
mutation_encoding_df['feature_changes'] = mutation_encoding_df.apply(
    lambda row: process_mutation_features(row, amino_acid_features), 
    axis=1
)

# 아미노산 특성 차이를 개별 열로 분리하고 통계 계산
feature_columns = ['hydrophobicity', 'polarity', 'mw', 'pI', 'charge']
for col in feature_columns:
    mutation_encoding_df[f'{col}_min'] = mutation_encoding_df['feature_changes'].apply(lambda x: np.min(x[col]) if x[col] else np.nan)
    mutation_encoding_df[f'{col}_max'] = mutation_encoding_df['feature_changes'].apply(lambda x: np.max(x[col]) if x[col] else np.nan)
    mutation_encoding_df[f'{col}_mean'] = mutation_encoding_df['feature_changes'].apply(lambda x: np.mean(x[col]) if x[col] else np.nan)
    mutation_encoding_df[f'{col}_std'] = mutation_encoding_df['feature_changes'].apply(lambda x: np.std(x[col]) if len(x[col]) > 1 else 0)

mutation_encoding_df['status'] = mutation_encoding_df['feature_changes'].apply(lambda x: x['status'])

# 'feature_changes' 열 제거
mutation_encoding_df = mutation_encoding_df.drop('feature_changes', axis=1)

# 결과 CSV로 저장
output_filepath = 'data/processed/test/test_mutation_encoding.csv'
mutation_encoding_df.to_csv(output_filepath, index=False)
print("Encoding completed successfully with amino acid feature differences.")