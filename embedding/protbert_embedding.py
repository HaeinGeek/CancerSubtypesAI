import pandas as pd
from tqdm import tqdm
import logging
from transformers import BertTokenizer, BertModel
import torch

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

def extract_full_embedding(sequence, tokenizer, model, device):
    try:
        inputs = tokenizer(sequence, 
                           return_tensors='pt', 
                           max_length=1024,
                           padding='max_length', 
                           truncation=True)
        
        inputs = {k: v.to(device) for k, v in inputs.items()}
        
        with torch.no_grad():
            outputs = model(**inputs)
            last_hidden_states = outputs.last_hidden_state
            attention_mask = inputs['attention_mask']
        
        return {
            'last_hidden_states': last_hidden_states.squeeze(0).cpu(),
            'attention_mask': attention_mask.squeeze(0).cpu()
        }
    except Exception as e:
        logger.error(f"Error in extract_full_embedding: {str(e)}")
        logger.error(f"Sequence causing error: {sequence}")
        raise

# Function to process and save embeddings
def process_and_save_embeddings(df, tokenizer, model, device, output_file, is_mutant=False):
    embedding_dict = {}
    skipped_count = 0
    
    for idx, row in tqdm(df.iterrows(), total=len(df), desc="Processing sequences"):
        try:
            isoform_id = row['isoform_id']
            sequence = row['mut_seq'] if is_mutant else row['wt_seq']
            embedding = extract_full_embedding(sequence, tokenizer, model, device)
            
            if is_mutant:
                mut_str = row['mutation_str']
                if isoform_id not in embedding_dict:
                    embedding_dict[isoform_id] = {}
                embedding_dict[isoform_id][mut_str] = embedding
            else:
                embedding_dict[isoform_id] = embedding
                
        except Exception as e:
            logger.error(f"Error processing row {idx}: {str(e)}")
            skipped_count += 1
    
    torch.save(embedding_dict, output_file)
    logger.info(f"Saved {len(embedding_dict)} embeddings to {output_file}")
    logger.info(f"Skipped {skipped_count} rows due to errors or NaN values")

    return embedding_dict