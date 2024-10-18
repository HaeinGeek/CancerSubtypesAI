import sys
import os
import pandas as pd
import numpy as np

# 프로젝트 루트 디렉토리를 sys.path에 추가
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(project_root)

from encoder.mutation_encoder import (
    classify_and_aggregate_mutations, 
    parse_multiple_mutations,
    process_mutation_features
)

def main():
    import os

    # 입력 및 출력 파일 경로 설정
    input_filepath = 'data/train.csv'
    amino_acid_filepath = 'data/amino_acid_features.csv'
    output_filepath = 'data/processed/train/train_mutation_encoding.csv'

    # 파일 존재 여부 확인
    if not os.path.exists(input_filepath):
        print(f"Input file '{input_filepath}' not found.")
        return

    if not os.path.exists(amino_acid_filepath):
        print(f"Amino acid feature file '{amino_acid_filepath}' not found.")
        return

    # 데이터 로드
    train_df = pd.read_csv(input_filepath)
    amino_acid_features = pd.read_csv(amino_acid_filepath).set_index('amino acid')

    # 데이터 변환 및 중복 제거
    melted_df = train_df.melt(id_vars=['ID', 'SUBCLASS'], var_name='gene', value_name='mutation_str')
    unique_mutations_df = melted_df[['gene', 'mutation_str']].drop_duplicates(subset=['gene', 'mutation_str'])

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
    mutation_encoding_df.to_csv(output_filepath, index=False)
    print("Encoding completed successfully with amino acid feature differences.")

if __name__ == "__main__":
    main()