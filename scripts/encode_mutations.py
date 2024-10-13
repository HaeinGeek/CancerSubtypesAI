import sys
import os

# Add the project root directory to sys.path
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(project_root)

from encoder.mutation_encoder import process_mutations
import pandas as pd

def main():
    import os

    # CSV 파일 읽기
    input_filepath = 'data/train.csv'  # 실제 파일 경로로 수정하세요
    if not os.path.exists(input_filepath):
        print(f"Input file '{input_filepath}' not found.")
        return

    train_df = pd.read_csv(input_filepath)

    train_type = train_df.copy()
    train_position = train_df.copy()

    for column in train_df.columns:
        if column not in ['ID', 'SUBCLASS']:  # ID와 SUBCLASS 열 제외
            # process_mutations 함수가 변이 종류와 변이 위치 두 개의 값을 반환하므로, 이를 각각 저장
            train_type[column], train_position[column] = zip(*train_df[column].apply(process_mutations))

    # 결과 저장
    train_type.to_csv('data/processed/train_mut_type_encoding.csv', index=False)
    train_position.to_csv('data/processed/train_mut_position_encoding.csv', index=False)

    # 변이 위치의 개수 인코딩
    df_pos_len = train_position.copy()
    for column in df_pos_len.columns:
        if column not in ['ID', 'SUBCLASS']:  # ID와 SUBCLASS 열 제외
            df_pos_len[column] = df_pos_len[column].apply(lambda x: len(x) if isinstance(x, list) else 0)

    df_pos_len.to_csv('data/processed/train_mut_position_len.csv', index=False)

    print("Encoding completed successfully.")

if __name__ == "__main__":
    main()