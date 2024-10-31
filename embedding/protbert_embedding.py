import pandas as pd
import numpy as np
import torch
from tqdm.auto import tqdm
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

def apply_pooling(last_hidden_states, attention_mask, pooling_strategy='mean_max_cls'):
    if pooling_strategy == 'mean_max_cls':
        mask = attention_mask.unsqueeze(-1).expand(last_hidden_states.size()).float()
        masked_embeddings = last_hidden_states * mask

        # Mean pooling
        summed = torch.sum(masked_embeddings, dim=1)
        counts = torch.clamp(mask.sum(dim=1), min=1e-9)
        mean_pooled = summed / counts

        # Max pooling
        masked_embeddings[mask == 0] = -1e9  # Assign a very small value to padding positions
        max_pooled, _ = torch.max(masked_embeddings, dim=1)

        # CLS token
        cls_token = last_hidden_states[:, 0, :]

        # Combine Mean, Max, and CLS token
        final_embedding = torch.cat([mean_pooled, max_pooled, cls_token], dim=1)

    elif pooling_strategy == 'preserve_length':
        # Apply Mean and Max pooling while preserving sequence length
        mean_pooled = last_hidden_states.mean(dim=-1)  # Shape: (batch_size, seq_len)
        max_pooled, _ = last_hidden_states.max(dim=-1)  # Shape: (batch_size, seq_len)
        final_embedding = torch.cat([mean_pooled, max_pooled], dim=1)  # Shape: (batch_size, seq_len * 2)

    elif pooling_strategy == 'full':
        # Return the entire last_hidden_states
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

    # Apply pooling
    pooled_embeddings = apply_pooling(last_hidden_states, attention_mask, pooling_strategy)
    return pooled_embeddings

def get_existing_keys(embedding_files):
    existing_keys = set()
    for file in embedding_files:
        with h5py.File(file, 'r') as h5f:
            keys = list(h5f.keys())
            existing_keys.update(keys)
    return existing_keys

def process_embeddings(df, tokenizer, model, device, output_file_base, existing_keys=None, is_mutant=False, pooling_strategy='mean_max_cls', batch_size=32):
    # Create 'key' column
    if is_mutant:
        df['key'] = df.apply(lambda row: f"{row['isoform_id']}_{row['mutation_str']}", axis=1)
    else:
        df['key'] = df['isoform_id']

    # If existing_keys is not None, filter out sequences whose keys are in existing_keys
    if existing_keys is not None:
        initial_len = len(df)
        df = df[~df['key'].isin(existing_keys)].reset_index(drop=True)
        logger.info(f"Filtered out {initial_len - len(df)} sequences already in existing embeddings.")

    sequences = df['mut_seq' if is_mutant else 'wt_seq'].tolist()
    keys = df['key'].tolist()

    total_size = 0  # Total size of saved data
    file_index = 1
    output_file = f"{output_file_base}_{file_index}.h5"
    h5f = h5py.File(output_file, 'w')  # 'w' mode to create a new file

    for i in tqdm(range(0, len(sequences), batch_size), desc="Processing batches"):
        batch_sequences = sequences[i:i+batch_size]
        batch_keys = keys[i:i+batch_size]

        try:
            batch_embeddings = extract_embeddings(
                batch_sequences, tokenizer, model, device, pooling_strategy
            )
            batch_embeddings = np.array(batch_embeddings)

            # Save embeddings to temporary file to estimate size
            temp_h5f = h5py.File('temp.h5', 'w')
            for j, embedding in enumerate(batch_embeddings):
                key = batch_keys[j]
                temp_h5f.create_dataset(key, data=embedding, compression='gzip', compression_opts=9)
            temp_h5f.close()
            temp_file_size = os.path.getsize('temp.h5')

            # Check file size limit
            # if total_size + temp_file_size > 10 * 1024 ** 3:  # 10GB limit
            if total_size + temp_file_size > 100 * 1024 ** 2:  # 100MB limit
                h5f.close()
                os.remove('temp.h5')
                file_index += 1
                output_file = f"{output_file_base}_{file_index}.h5"
                h5f = h5py.File(output_file, 'w')
                total_size = 0

                # Save embeddings to temp file again
                temp_h5f = h5py.File('temp.h5', 'w')
                for j, embedding in enumerate(batch_embeddings):
                    key = batch_keys[j]
                    temp_h5f.create_dataset(key, data=embedding, compression='gzip', compression_opts=9)
                temp_h5f.close()
                temp_file_size = os.path.getsize('temp.h5')

            # Copy from temp file to actual HDF5 file
            temp_h5f = h5py.File('temp.h5', 'r')
            for key in temp_h5f.keys():
                temp_h5f.copy(key, h5f)
            temp_h5f.close()
            os.remove('temp.h5')

            total_size += temp_file_size

        except Exception as e:
            logger.error(f"Error processing batch starting at index {i}: {str(e)}")

    h5f.close()
    logger.info(f"Embeddings saved to {output_file}")