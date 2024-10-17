import os
import sys
import pandas as pd
import numpy as np
import ast
from transformers import BertTokenizer, BertModel
import torch

# 프로젝트 루트 디렉토리를 sys.path에 추가
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(project_root)

# 데이터 로드
mut_seq_filepath = 'data/processed/mutant_seq_unique.csv'
mut_encoding_filepath = 'data/processed/train_mutation_encoding.csv'
mut_seq_df = pd.read_csv(mut_seq_filepath)
mut_encoding_df = pd.read_csv(mut_encoding_filepath)

# 특정 열의 문자열을 리스트로 변환
columns = ['origin', 'position','mutant']
for column in columns:
    mut_encoding_df[column] = mut_encoding_df[column].apply(ast.literal_eval)

# 데이터프레임 병합
merged_df = mut_encoding_df.merge(mut_seq_df, on=['gene','mutation_str'], how='left')

# 'status_prot' 열 초기화 및 업데이트
merged_df['status_prot'] = 0
status_prot_1_types = ['Frameshift']

for type in status_prot_1_types:
    merged_df.loc[merged_df['type'] == type, 'status_prot'] = 1

# 'Complex_mutation' 타입에 대한 추가 처리
rows = merged_df.loc[merged_df['type'] == 'Complex_mutation'].iterrows()
for idx, row in rows:
    if 'fs' in row['mutant']:
        merged_df.at[idx, 'status_prot'] = 1

# 필요한 데이터만 추출
merged_df_prot = merged_df.loc[(merged_df.type != 'WT') & (merged_df.type != 'Silent_Missense') & (merged_df.type != 'Silent_Nonsense')]
merged_df_prot = merged_df_prot.loc[(merged_df.status_prot == 0)]

# WT 단백질 데이터 추출
wt_prot = merged_df_prot.drop_duplicates(subset=['isoform_id'])
wt_prot = wt_prot[['gene', 'isoform_id', 'wt_seq']]

# 변이 단백질 데이터 추출
mut_prot = merged_df_prot.drop_duplicates(subset=['isoform_id', 'mutation_str'])
mut_prot = mut_prot[['gene', 'mutation_str', 'isoform_id', 'mut_seq']]

# ProtBERT 토크나이저와 모델 로드
tokenizer = BertTokenizer.from_pretrained('Rostlab/prot_bert')
model = BertModel.from_pretrained('Rostlab/prot_bert')

# 임베딩 추출 함수
def extract_embedding(sequence):
    inputs = tokenizer(sequence, 
                       return_tensors='pt', 
                       max_length=1024,
                       padding=True, 
                       truncation=True)
    with torch.no_grad():
        embedding = model(**inputs).last_hidden_state
    return embedding

# WT 임베딩 추출
wt_embedding_dict = {}
for idx, row in wt_prot.iterrows():
    isoform_id = row['isoform_id']
    wt_embedding = extract_embedding(row['wt_seq'])
    wt_embedding_dict[isoform_id] = wt_embedding

# WT 임베딩 저장
torch.save(wt_embedding_dict, 'data/processed/wt_embedding_tensor.pt')

# 변이 임베딩 추출
mut_embedding_dict = {}
for idx, row in mut_prot.iterrows():
    mut_str = row['mutation_str']
    isoform_id = row['isoform_id']
    mut_embedding = extract_embedding(row['mut_seq'])

    if isoform_id not in mut_embedding_dict:
        mut_embedding_dict[isoform_id] = {}
    
    mut_embedding_dict[isoform_id][mut_str] = mut_embedding

# 변이 임베딩 저장
torch.save(mut_embedding_dict, 'data/processed/mut_embedding_tensor.pt')

print("Embeddings extraction completed and saved.")