import pandas as pd
import numpy as np
import torch
from tqdm import tqdm
import logging
from transformers import BertTokenizer, BertModel
import h5py
import os

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def extract_data_for_embedding(mutation_df):
    # 임베딩 계산할 항목만 추출
    df_prot = mutation_df.loc[(mutation_df.type != 'WT') &
                    (mutation_df.type != 'Silent_Missense') &
                    (mutation_df.type != 'Silent_Nonsense') &
                    (mutation_df.status_prot == 0)
                    ]
    # 변이 전(wt), 후(mut) 서열 분리
    wt_prot = df_prot.drop_duplicates(subset=['isoform_id'])[['gene', 'isoform_id', 'wt_seq']]
    mut_prot = df_prot.drop_duplicates(subset=['isoform_id', 'mutation_str'])[['gene', 'mutation_str', 'isoform_id', 'mut_seq']]

    return wt_prot, mut_prot

def load_protbert():
    # ProtBERT tokenizer & model 로딩
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    logger.info(f"Using device: {device}")

    tokenizer = BertTokenizer.from_pretrained('Rostlab/prot_bert', do_lower_case=False)
    model = BertModel.from_pretrained("Rostlab/prot_bert").to(device)
    
    return tokenizer, model, device

# def extract_embeddings(sequences, tokenizer, model, device, batch_size=32):
#     embeddings = []
#     for i in range(0, len(sequences), batch_size):
#         batch = sequences[i:i+batch_size]
#         inputs = tokenizer(batch, 
#                            return_tensors='pt', 
#                            max_length=1024,
#                            padding='max_length', 
#                            truncation=True)
        
#         inputs = {k: v.to(device) for k, v in inputs.items()}
        
#         with torch.no_grad():
#             outputs = model(**inputs)
#             last_hidden_states = outputs.last_hidden_state
#             attention_mask = inputs['attention_mask']
        
#         embeddings.extend([
#             {
#                 'last_hidden_states': lhs.cpu(),
#                 'attention_mask': am.cpu()
#             }
#             for lhs, am in zip(last_hidden_states, attention_mask)
#         ])
    
#     return embeddings

# def process_embeddings(df, tokenizer, model, device, is_mutant=False, batch_size=32):
#     embedding_dict = {}
#     skipped_count = 0
    
#     sequences = df['mut_seq' if is_mutant else 'wt_seq'].tolist()
#     isoform_ids = df['isoform_id'].tolist()
    
#     if is_mutant:
#         mutation_strs = df['mutation_str'].tolist()
    
#     for i in tqdm(range(0, len(sequences), batch_size), desc="Processing batches"):
#         batch_sequences = sequences[i:i+batch_size]
#         batch_isoform_ids = isoform_ids[i:i+batch_size]
        
#         try:
#             batch_embeddings = extract_embeddings(batch_sequences, tokenizer, model, device, batch_size)
            
#             for j, embedding in enumerate(batch_embeddings):
#                 isoform_id = batch_isoform_ids[j]
                
#                 if is_mutant:
#                     mut_str = mutation_strs[i+j]
#                     if isoform_id not in embedding_dict:
#                         embedding_dict[isoform_id] = {}
#                     embedding_dict[isoform_id][mut_str] = embedding
#                 else:
#                     embedding_dict[isoform_id] = embedding
                    
#         except Exception as e:
#             logger.error(f"Error processing batch starting at index {i}: {str(e)}")
#             skipped_count += len(batch_sequences)
    
#     logger.info(f"Skipped {skipped_count} sequences due to errors")

#     return embedding_dict

# def save_embeddings_to_zip(embedding_dict, output_file):
#     with zipfile.ZipFile(output_file + '.zip', 'w', zipfile.ZIP_DEFLATED) as zf:
#         buffer = io.BytesIO()
#         torch.save(embedding_dict, buffer, _use_new_zipfile_serialization=False)
#         zf.writestr('embeddings.pt', buffer.getvalue())
    
#     logger.info(f"Saved {len(embedding_dict)} embeddings to {output_file}.zip")

# def load_compressed_embeddings(file_path):
#     with zipfile.ZipFile(file_path, 'r') as zf:
#         with zf.open('embeddings.pt') as f:
#             buffer = io.BytesIO(f.read())
#             return torch.load(buffer)

def apply_pooling(last_hidden_states, attention_mask, pooling_strategy='mean_max_cls'):
    if pooling_strategy == 'mean_max_cls':
        mask = attention_mask.unsqueeze(-1).expand(last_hidden_states.size()).float()
        masked_embeddings = last_hidden_states * mask

        # Mean pooling
        summed = torch.sum(masked_embeddings, dim=1)
        counts = torch.clamp(mask.sum(dim=1), min=1e-9)
        mean_pooled = summed / counts

        # Max pooling
        masked_embeddings[mask == 0] = -1e9  # 패딩 위치에 매우 작은 값 할당
        max_pooled, _ = torch.max(masked_embeddings, dim=1)

        # CLS 토큰
        cls_token = last_hidden_states[:, 0, :]

        # Mean, Max, CLS 토큰 결합
        final_embedding = torch.cat([mean_pooled, max_pooled, cls_token], dim=1)

    elif pooling_strategy == 'preserve_length':
        # 시퀀스 길이를 유지하며 각 위치에 대해 Mean 및 Max 풀링 적용
        mean_pooled = last_hidden_states.mean(dim=-1)  # Shape: (batch_size, seq_len)
        max_pooled, _ = last_hidden_states.max(dim=-1)  # Shape: (batch_size, seq_len)
        final_embedding = torch.cat([mean_pooled, max_pooled], dim=1)  # Shape: (batch_size, seq_len * 2)

    elif pooling_strategy == 'full':
        # 전체 last_hidden_states 반환
        final_embedding = last_hidden_states.view(last_hidden_states.size(0), -1)  # Flatten
    else:
        raise ValueError("Invalid pooling strategy")

    return final_embedding.cpu().numpy()

def extract_embeddings(sequences, tokenizer, model, device, pooling_strategy='mean_max_cls'):
    inputs = tokenizer(
        sequences,
        return_tensors='pt',
        max_length=1024,
        padding='max_length',
        truncation=True
    )
    inputs = {k: v.to(device) for k, v in inputs.items()}

    with torch.no_grad():
        outputs = model(**inputs)
        last_hidden_states = outputs.last_hidden_state  # (batch_size, seq_len, hidden_size)
        attention_mask = inputs['attention_mask']       # (batch_size, seq_len)

    # 풀링 적용
    pooled_embeddings = apply_pooling(last_hidden_states, attention_mask, pooling_strategy)
    return pooled_embeddings

def process_embeddings(df, tokenizer, model, device, output_file_base, is_mutant=False, pooling_strategy='mean_max_cls', batch_size=32):
    sequences = df['mut_seq' if is_mutant else 'wt_seq'].tolist()
    isoform_ids = df['isoform_id'].tolist()

    if is_mutant:
        mutation_strs = df['mutation_str'].tolist()

    total_size = 0  # 저장된 데이터의 총 크기 추적
    file_index = 1  # 파일 인덱스
    output_file = f"{output_file_base}_{file_index}.h5"
    h5f = h5py.File(output_file, 'w')  # 새로운 파일 생성 시 'w' 모드 사용

    for i in tqdm(range(0, len(sequences), batch_size), desc="Processing batches"):
        batch_sequences = sequences[i:i+batch_size]
        batch_isoform_ids = isoform_ids[i:i+batch_size]

        if is_mutant:
            batch_mutation_strs = mutation_strs[i:i+batch_size]

        try:
            batch_embeddings = extract_embeddings(
                batch_sequences, tokenizer, model, device, pooling_strategy
            )
            batch_embeddings = np.array(batch_embeddings)

            # 임시 파일에 저장하여 크기 계산
            temp_h5f = h5py.File('temp.h5', 'w')
            for j, embedding in enumerate(batch_embeddings):
                if is_mutant:
                    key = f"{batch_isoform_ids[j]}_{batch_mutation_strs[j]}"
                else:
                    key = batch_isoform_ids[j]
                temp_h5f.create_dataset(key, data=embedding, compression='gzip', compression_opts=9)
            temp_h5f.close()
            temp_file_size = os.path.getsize('temp.h5')

            # 파일 크기 제한 체크
            # if total_size + temp_file_size > 10 * 1024 ** 3:  # 10GB 초과 시 새로운 파일 생성
            if total_size + temp_file_size > 100 * 1024 ** 2:  # 100MB 초과 시 새로운 파일 생성
                h5f.close()
                os.remove('temp.h5')  # 임시 파일 삭제
                file_index += 1
                output_file = f"{output_file_base}_{file_index}.h5"
                h5f = h5py.File(output_file, 'w')
                total_size = 0

                # 다시 임시 파일에 저장
                temp_h5f = h5py.File('temp.h5', 'w')
                for j, embedding in enumerate(batch_embeddings):
                    if is_mutant:
                        key = f"{batch_isoform_ids[j]}_{batch_mutation_strs[j]}"
                    else:
                        key = batch_isoform_ids[j]
                    temp_h5f.create_dataset(key, data=embedding, compression='gzip', compression_opts=9)
                temp_h5f.close()
                temp_file_size = os.path.getsize('temp.h5')

            # 임시 파일의 내용을 실제 HDF5 파일로 복사
            temp_h5f = h5py.File('temp.h5', 'r')
            for key in temp_h5f.keys():
                temp_h5f.copy(key, h5f)
            temp_h5f.close()
            os.remove('temp.h5')  # 임시 파일 삭제

            total_size += temp_file_size

        except Exception as e:
            logger.error(f"Error processing batch starting at index {i}: {str(e)}")

    h5f.close()
    logger.info(f"Embeddings saved to {output_file}")