import os
import sys
import logging
import pandas as pd
import ast

def setup_environment():
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    sys.path.append(project_root)

def load_logger():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    logger = logging.getLogger(__name__)
    return logger

def load_full_feature_data(train=True):
    if train:
        mut_seq_filepath = 'data/processed/train/mutant_seq_unique.csv'
        mut_encoding_filepath = 'data/processed/train/train_mutation_encoding.csv'
    else:
        mut_seq_filepath = 'data/processed/test/mutant_seq_unique.csv'
        mut_encoding_filepath = 'data/processed/test/test_mutation_encoding.csv'

    mut_seq_df = pd.read_csv(mut_seq_filepath)
    mut_encoding_df = pd.read_csv(mut_encoding_filepath)

    columns = ['origin', 'position','mutant']
    for column in columns:
        mut_encoding_df[column] = mut_encoding_df[column].apply(ast.literal_eval)

    mutation_df = mut_encoding_df.merge(mut_seq_df, on=['gene','mutation_str'], how='left')

    # 계산 불가 -> staus_prot = 1
    mutation_df['status_prot'] = 0
    mutation_df.loc[mutation_df['type'] == 'Frameshift', 'status_prot'] = 1
    mutation_df.loc[(mutation_df['type'] == 'Complex_mutation') & (mutation_df['mutant'].apply(lambda x: 'fs' in x)), 'status_prot'] = 1
    mutation_df.loc[(mutation_df['isoform_id'].isna()), 'status_prot'] = 1

    return mutation_df

def load_full_seq_data():
    seq_train = load_full_feature_data()
    seq_test = load_full_feature_data(train=False)

    mutation_df = pd.concat([seq_train, seq_test]).drop_duplicates(subset=['gene', 'mutation_str'])
    mutation_df.reset_index(drop=True, inplace=True)

    wt_seq_df = mutation_df.loc[(mutation_df.mutation_str == 'WT'),['gene', 'mutation_str', 'type','isoform_id', 'wt_seq']]
    wt_seq_df.columns = ['gene', 'mutation_str', 'type','isoform_id', 'seq']
    mut_seq_df = mutation_df.loc[(mutation_df.mutation_str != 'WT'),['gene', 'mutation_str', 'type', 'isoform_id', 'mut_seq']]
    mut_seq_df.columns = ['gene', 'mutation_str', 'type','isoform_id', 'seq']

    seq_df = pd.concat([wt_seq_df, mut_seq_df])
    seq_df = seq_df.drop_duplicates(subset=['gene', 'mutation_str'])

    return seq_df

def download_file(file_path):
    from google.colab import files
    import os
    """
    Google Colab에서 특정 경로의 파일을 다운로드합니다.

    :param file_path: 다운로드할 파일의 경로
    """
    if os.path.exists(file_path):
        files.download(file_path)
        # print(f"{file_path} 파일이 다운로드되었습니다.")
    else:
        print(f'Error: "{file_path}" not found.')

def txt_to_list(file_path):
    # 파일에서 읽어서 리스트에 저장
    li = []
    with open(file_path, 'r') as f:
        # 각 줄을 읽어서 리스트에 추가
        li = [line.strip() for line in f]
        # 또는
        # for line in f:
        #     high_freq_genes.append(line.strip())

    print(f"Length of list: {len(li)}")
    print(f"First 5 components: {li[:5]}")

    return li

def list_to_txt(li, file_path):
    """
    리스트의 각 항목을 텍스트 파일의 각 줄에 저장하는 함수
    
    Parameters:
    li (list): 저장할 리스트
    file_path (str): 저장할 파일 경로
    """
    # 파일을 쓰기 모드로 열고 리스트의 각 항목을 줄별로 저장
    with open(file_path, 'w') as f:
        for item in li:
            f.write(str(item) + '\n')
    
    print(f"File saved: {file_path}")
    print(f"Length of list: {len(li)}")