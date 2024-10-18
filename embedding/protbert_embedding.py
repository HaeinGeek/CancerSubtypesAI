import pandas as pd
from tqdm import tqdm
import logging
from transformers import BertTokenizer, BertModel
import torch
import zipfile
import io

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

# def extract_full_embedding(sequence, tokenizer, model, device):
#     try:
#         inputs = tokenizer(sequence, 
#                            return_tensors='pt', 
#                            max_length=1024,
#                            padding='max_length', 
#                            truncation=True)
        
#         inputs = {k: v.to(device) for k, v in inputs.items()}
        
#         with torch.no_grad():
#             outputs = model(**inputs)
#             last_hidden_states = outputs.last_hidden_state
#             attention_mask = inputs['attention_mask']
        
#         return {
#             'last_hidden_states': last_hidden_states.squeeze(0).cpu(),
#             'attention_mask': attention_mask.squeeze(0).cpu()
#         }
#     except Exception as e:
#         logger.error(f"Error in extract_full_embedding: {str(e)}")
#         logger.error(f"Sequence causing error: {sequence}")
#         raise

# def process_embeddings(df, tokenizer, model, device, is_mutant=False):
#     embedding_dict = {}
#     skipped_count = 0
    
#     for idx, row in tqdm(df.iterrows(), total=len(df), desc="Processing sequences"):
#         try:
#             isoform_id = row['isoform_id']
#             sequence = row['mut_seq'] if is_mutant else row['wt_seq']
#             embedding = extract_full_embedding(sequence, tokenizer, model, device)
            
#             if is_mutant:
#                 mut_str = row['mutation_str']
#                 if isoform_id not in embedding_dict:
#                     embedding_dict[isoform_id] = {}
#                 embedding_dict[isoform_id][mut_str] = embedding
#             else:
#                 embedding_dict[isoform_id] = embedding
                
#         except Exception as e:
#             logger.error(f"Error processing row {idx}: {str(e)}")
#             skipped_count += 1
    
#     logger.info(f"Skipped {skipped_count} rows due to errors or NaN values")

#     return embedding_dict

import torch
from tqdm import tqdm

def extract_embeddings(sequences, tokenizer, model, device, batch_size=32):
    embeddings = []
    for i in range(0, len(sequences), batch_size):
        batch = sequences[i:i+batch_size]
        inputs = tokenizer(batch, 
                           return_tensors='pt', 
                           max_length=1024,
                           padding='max_length', 
                           truncation=True)
        
        inputs = {k: v.to(device) for k, v in inputs.items()}
        
        with torch.no_grad():
            outputs = model(**inputs)
            last_hidden_states = outputs.last_hidden_state
            attention_mask = inputs['attention_mask']
        
        embeddings.extend([
            {
                'last_hidden_states': lhs.cpu(),
                'attention_mask': am.cpu()
            }
            for lhs, am in zip(last_hidden_states, attention_mask)
        ])
    
    return embeddings

def process_embeddings(df, tokenizer, model, device, is_mutant=False, batch_size=32):
    embedding_dict = {}
    skipped_count = 0
    
    sequences = df['mut_seq' if is_mutant else 'wt_seq'].tolist()
    isoform_ids = df['isoform_id'].tolist()
    
    if is_mutant:
        mutation_strs = df['mutation_str'].tolist()
    
    for i in tqdm(range(0, len(sequences), batch_size), desc="Processing batches"):
        batch_sequences = sequences[i:i+batch_size]
        batch_isoform_ids = isoform_ids[i:i+batch_size]
        
        try:
            batch_embeddings = extract_embeddings(batch_sequences, tokenizer, model, device, batch_size)
            
            for j, embedding in enumerate(batch_embeddings):
                isoform_id = batch_isoform_ids[j]
                
                if is_mutant:
                    mut_str = mutation_strs[i+j]
                    if isoform_id not in embedding_dict:
                        embedding_dict[isoform_id] = {}
                    embedding_dict[isoform_id][mut_str] = embedding
                else:
                    embedding_dict[isoform_id] = embedding
                    
        except Exception as e:
            logger.error(f"Error processing batch starting at index {i}: {str(e)}")
            skipped_count += len(batch_sequences)
    
    logger.info(f"Skipped {skipped_count} sequences due to errors")

    return embedding_dict

def save_embeddings_to_zip(embedding_dict, output_file):
    with zipfile.ZipFile(output_file + '.zip', 'w', zipfile.ZIP_DEFLATED) as zf:
        buffer = io.BytesIO()
        torch.save(embedding_dict, buffer, _use_new_zipfile_serialization=False)
        zf.writestr('embeddings.pt', buffer.getvalue())
    
    logger.info(f"Saved {len(embedding_dict)} embeddings to {output_file}.zip")

def load_compressed_embeddings(file_path):
    with zipfile.ZipFile(file_path, 'r') as zf:
        with zf.open('embeddings.pt') as f:
            buffer = io.BytesIO(f.read())
            return torch.load(buffer)