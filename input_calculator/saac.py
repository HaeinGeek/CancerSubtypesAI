import os
import glob
import numpy as np
import pandas as pd
from collections import defaultdict


aa_codes = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 
            'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'fs']

def calculate_saac(seq_df, cutoff_ratio = 0.2):
    # 디렉토리가 존재하는지 확인하고 없으면 생성
    # output_dir = 'data/processed/saac'
    # os.makedirs(output_dir, exist_ok=True)

    # 새로운 열을 미리 생성
    for region in ['n_term', 'mid', 'c_term']:
        for code in aa_codes:
            seq_df[f'{region}_{code}'] = 0

    for index, row in seq_df.iterrows():
        seq = row['seq']
        # 문자열인지 확인
        if isinstance(seq, str):
            seq_len = len(seq)
            mut_str = row['mutation_str']
            # 길이가 5보다 작으면 모두 N-term에 포함
            if seq_len < 5:
                n_term_seq = list(seq)
                mid_seq = []
                c_term_seq = []
            
            else:
                # fs 변이인 경우
                if isinstance(mut_str, str) and 'fs' in mut_str:
                    term_len = round(cutoff_ratio * seq_len)
                    seq_list = list(seq)
                    n_term_seq = seq_list[:term_len]
                    mid_seq = seq_list[term_len:-term_len]
                    c_term_seq = seq_list[-term_len:] + ['fs']

                # fs 아닌 경우
                else:
                    term_len = round(cutoff_ratio * seq_len)
                    seq_list = list(seq)
                    n_term_seq = seq_list[:term_len]
                    mid_seq = seq_list[term_len:-term_len]
                    c_term_seq = seq_list[-term_len:] 

            # n_term, mid, c_term 아미노산 빈도 계산
            for code in aa_codes:
                # 분모가 0인 경우 0을 반환
                seq_df.at[index, f'n_term_{code}'] = float((n_term_seq.count(code) / len(n_term_seq))*100) if len(n_term_seq) > 0 else 0.0
                seq_df.at[index, f'mid_{code}'] = float((mid_seq.count(code) / len(mid_seq))*100) if len(mid_seq) > 0 else 0.0
                seq_df.at[index, f'c_term_{code}'] = float((c_term_seq.count(code) / len(c_term_seq))*100) if len(c_term_seq) > 0 else 0.0

        # 문자열이 아닌 경우 모든 빈도 0
        else:
            for region in ['n_term', 'mid', 'c_term']:
                for code in aa_codes:
                    seq_df.at[index, f'{region}_{code}'] = 0.0

    # CSV로 저장
    # file_name = f'{output_dir}/saac_seq_cutoff_{int(cutoff_ratio*100)}.csv'
    # seq_df.to_csv(file_name, index=False)
    return seq_df

# def save_saac(train_df, saac_df, gene_list, output_dir):
#     # file_dir = 'data/processed/saac/train/'
#     # file_path = f'data/processed/saac/train/saac_{gene}.csv'
    
#     # 디렉토리가 존재하는지 확인하고 없으면 생성
#     os.makedirs(output_dir, exist_ok=True)
#     file_path = output_dir + f'saac_{gene}.csv'
#     saac_array = saac_df.drop(['type', 'isoform_id', 'seq'], axis=1).values
#     for gene in gene_list:
#         # 데이터프레임 초기화
#         df = train_df[['ID']].copy()
        
#         # 열 생성
#         df[f'{gene}_seq_len'] = 0
#         for region in ['n_term', 'mid', 'c_term']:
#             for code in aa_codes:
#                 df[f'{gene}_{region}_{code}'] = 0

#         selected_cols = df.columns[1:]

#         # 유전자와 돌연변이 문자열에 대한 SAAC 값을 미리 계산
#         saac_dict = defaultdict(lambda: np.zeros(len(selected_cols)))
#         gene_mask = (saac_array[:, 0] == gene)
#         for mut_str in train_df[gene].unique():
#             mut_mask = (saac_array[:, 1] == mut_str)
#             arr = saac_array[gene_mask & mut_mask][:, 2:]
#             if len(arr) > 0:
#                 saac_dict[mut_str] = arr[0]

#         # 데이터프레임 채우기
#         for idx, row in train_df[[gene]].iterrows():
#             mut_str = row[gene]
#             df.loc[idx, selected_cols] = saac_dict[mut_str]

#         df.to_csv(file_path, index=False)
#         print(f'Save {gene} file as {file_path}.')

#     # 모든 SAAC CSV 파일 목록 가져오기
#     all_csv_files = glob.glob(os.path.join(output_dir, 'saac_*.csv'))

#     # 첫 번째 CSV 파일을 기준으로 데이터프레임 초기화
#     base_df = pd.read_csv(all_csv_files[0])

#     # 나머지 CSV 파일들을 순회하며 병합
#     for file in all_csv_files[1:]:
#         df = pd.read_csv(file)
#         # ID 열을 기준으로 병합
#         base_df = pd.merge(base_df, df, on='ID', how='outer')

#     base_df = pd.merge(train_df[['ID', 'SUBCLASS']], base_df, on='ID', how='outer')

#     # 결과 저장
#     # base_df.to_csv('data/processed/saac/train_saac_features.csv', index=False)
#     # print(f"Creat merged CSV file.\nTotal lines: {len(base_df)}\nTotal columns: {len(base_df.columns)}")
#     return base_df

def get_saac_features(input_df, saac_df, gene_list, output_dir, is_train=True):
    """
    Generate and save SAAC features for given dataset
    
    Parameters:
    -----------
    input_df : pandas.DataFrame
        Input dataframe (can be either train or test dataset)
    saac_df : pandas.DataFrame
        SAAC dataframe containing sequence information
    gene_list : list
        List of genes to process
    output_dir : str
        Directory to save output files
    is_train : bool, default=True
        Flag to indicate if input_df is training dataset
    
    Returns:
    --------
    pandas.DataFrame
        Merged dataframe with SAAC features
    """
    # 디렉토리가 존재하는지 확인하고 없으면 생성
    os.makedirs(output_dir, exist_ok=True)
    
    for gene in gene_list:
        file_path = output_dir + f'saac_{gene}.csv'
        saac_array = saac_df.drop(['type', 'isoform_id', 'seq'], axis=1).values
        
        # 데이터프레임 초기화
        df = input_df[['ID']].copy()
        
        # 열 생성
        df[f'{gene}_seq_len'] = 0
        for region in ['n_term', 'mid', 'c_term']:
            for code in aa_codes:
                df[f'{gene}_{region}_{code}'] = 0
        selected_cols = df.columns[1:]
        
        # 유전자와 돌연변이 문자열에 대한 SAAC 값을 미리 계산
        saac_dict = defaultdict(lambda: np.zeros(len(selected_cols)))
        gene_mask = (saac_array[:, 0] == gene)
        for mut_str in input_df[gene].unique():
            mut_mask = (saac_array[:, 1] == mut_str)
            arr = saac_array[gene_mask & mut_mask][:, 2:]
            if len(arr) > 0:
                saac_dict[mut_str] = arr[0]
        
        # 데이터프레임 채우기
        for idx, row in input_df[[gene]].iterrows():
            mut_str = row[gene]
            df.loc[idx, selected_cols] = saac_dict[mut_str]
        
        df.to_csv(file_path, index=False)
        print(f'Save {gene} file as {file_path}.')
    
    # 모든 SAAC CSV 파일 목록 가져오기
    all_csv_files = glob.glob(os.path.join(output_dir, 'saac_*.csv'))
    
    # 첫 번째 CSV 파일을 기준으로 데이터프레임 초기화
    base_df = pd.read_csv(all_csv_files[0])
    
    # 나머지 CSV 파일들을 순회하며 병합
    for file in all_csv_files[1:]:
        df = pd.read_csv(file)
        # ID 열을 기준으로 병합
        base_df = pd.merge(base_df, df, on='ID', how='outer')
    
    # train 데이터셋인 경우에만 SUBCLASS 열 추가
    if is_train:
        base_df = pd.merge(input_df[['ID', 'SUBCLASS']], base_df, on='ID', how='outer')
    
    return base_df
