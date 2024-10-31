import os
import glob
from tqdm.auto import tqdm
import numpy as np
import pandas as pd
from collections import defaultdict


aa_codes = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 
            'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'fs']

def calculate_saac(seq_df, cutoff_ratio=0.2):
    # 복사본 생성
    result_df = seq_df.copy()
    
    # 새로운 열 생성
    new_columns = [f'{region}_{code}' for region in ['n_term', 'mid', 'c_term'] for code in aa_codes]
    result_df[new_columns] = 0.0
    
    # tqdm으로 진행상황 표시
    for index, row in tqdm(result_df.iterrows(), 
                          total=len(result_df), 
                          desc="Calculating SAAC features",
                          position=0):
        seq = row['seq']
        # 문자열인지 확인
        if isinstance(seq, str):
            seq_len = len(seq)
            mut_str = row['mutation_str']
            
            # 길이가 5보다 작은 경우 처리
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
                result_df.loc[index, f'n_term_{code}'] = float((n_term_seq.count(code) / len(n_term_seq))*100) if len(n_term_seq) > 0 else 0.0
                result_df.loc[index, f'mid_{code}'] = float((mid_seq.count(code) / len(mid_seq))*100) if len(mid_seq) > 0 else 0.0
                result_df.loc[index, f'c_term_{code}'] = float((c_term_seq.count(code) / len(c_term_seq))*100) if len(c_term_seq) > 0 else 0.0
    
    return result_df

def get_saac_features(input_df, saac_df, gene_list, output_dir, is_train=True):
    # 디렉토리가 존재하는지 확인하고 없으면 생성
    os.makedirs(output_dir, exist_ok=True)
    
    # 전체 진행상황을 보여주는 tqdm
    for gene in tqdm(gene_list, desc="Overall Progress", position=0):
        file_path = os.path.join(output_dir, f'saac_{gene}.csv')
        saac_array = saac_df.drop(['type', 'isoform_id', 'seq'], axis=1).values
        
        # 데이터프레임 초기화
        df = input_df[['ID']].copy()
        
        # 열 생성
        df[f'{gene}_seq_len'] = 0
        new_columns = [f'{region}_{code}' for region in ['n_term', 'mid', 'c_term'] for code in aa_codes]
        df[new_columns] = 0.0
        selected_cols = df.columns[1:]

        # 유전자와 돌연변이 문자열에 대한 SAAC 값을 미리 계산
        saac_dict = defaultdict(lambda: np.zeros(len(selected_cols)))
        gene_mask = (saac_array[:, 0] == gene)
        
        # 유전자별 돌연변이 처리 진행상황을 보여주는 tqdm
        unique_mutations = input_df[gene].unique()
        for mut_str in tqdm(unique_mutations, 
                           desc=f"Processing {gene} mutations", 
                           position=1, 
                           leave=False):
            mut_mask = (saac_array[:, 1] == mut_str)
            arr = saac_array[gene_mask & mut_mask][:, 2:]
            if len(arr) > 0:
                saac_dict[mut_str] = arr[0]
        
        # 데이터프레임 채우기 진행상황을 보여주는 tqdm
        for idx, row in tqdm(input_df[[gene]].iterrows(), 
                           desc=f"Filling {gene} dataframe",
                           position=1,
                           leave=False,
                           total=len(input_df)):
            mut_str = row[gene]
            df.loc[idx, selected_cols] = saac_dict[mut_str]
        
        df.to_csv(file_path, index=False)
        # print(f'Save {gene} file as {file_path}')
    
    # 모든 SAAC CSV 파일 목록 가져오기
    all_csv_files = glob.glob(os.path.join(output_dir, 'saac_*.csv'))
    
    # 첫 번째 CSV 파일을 기준으로 데이터프레임 초기화
    base_df = pd.read_csv(all_csv_files[0])
    
    # 나머지 CSV 파일들을 순회하며 병합
    for file in tqdm(all_csv_files[1:], 
                    desc="Merging CSV files",
                    position=0):
        df = pd.read_csv(file)
        # ID 열을 기준으로 병합
        base_df = pd.merge(base_df, df, on='ID', how='outer')
    
    # train 데이터셋인 경우에만 SUBCLASS 열 추가
    if is_train:
        base_df = pd.merge(input_df[['ID', 'SUBCLASS']], base_df, on='ID', how='outer')
    
    return base_df