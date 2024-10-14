import sys
import os

# Add the project root directory to sys.path
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(project_root)

from encoder.mutation_encoder import classify_and_aggregate_mutations
from encoder.mutation_encoder import parse_multiple_mutations
import pandas as pd

def main():
    import os

    # CSV 파일 읽기
    input_filepath = 'data/train.csv'
    if not os.path.exists(input_filepath):
        print(f"Input file '{input_filepath}' not found.")
        return

    train_df = pd.read_csv(input_filepath)

    # 데이터 변환 및 중복 제거
    melted_df = train_df.melt(id_vars=['ID', 'SUBCLASS'], var_name='gene', value_name='mutation_str')
    unique_mutations_df = melted_df[['gene', 'mutation_str']].drop_duplicates(subset=['gene', 'mutation_str'])

    # 변이 인코딩 데이터프레임 생성
    mutation_encoding_df = unique_mutations_df.copy()

    # 변이 종류
    mutation_encoding_df['type'] = unique_mutations_df['mutation_str'].apply(classify_and_aggregate_mutations)

    # 변이 아미노산 및 변이 위치 처리
    mutation_encoding_df['origin'], mutation_encoding_df['position'], mutation_encoding_df['mutant'] = zip(
        *unique_mutations_df['mutation_str'].apply(parse_multiple_mutations)
        )
    
    # 변이 위치의 개수
    mutation_encoding_df['mut_num'] = mutation_encoding_df['position'].apply(lambda x: len(x) if isinstance(x, list) else 0)


    # 결과 저장
    output_filepath = 'data/processed/train_mutation_encoding.csv'
    mutation_encoding_df.to_csv(output_filepath, index=False)

    print("Encoding completed successfully.")

if __name__ == "__main__":
    main()