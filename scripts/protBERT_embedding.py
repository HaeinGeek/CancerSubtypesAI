import os
import sys
import pandas as pd
import ast
from transformers import BertTokenizer, BertModel
import torch

def setup_environment():
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    sys.path.append(project_root)

def load_full_feature_data():
    mut_seq_df = pd.read_csv('data/processed/mutant_seq_unique.csv')
    mut_encoding_df = pd.read_csv('data/processed/train_mutation_encoding.csv')
    
    for column in ['origin', 'position', 'mutant']:
        mut_encoding_df[column] = mut_encoding_df[column].apply(ast.literal_eval)
    
    merged_df =  mut_seq_df.merge(mut_encoding_df, on=['gene', 'mutation_str'], how='left')

    merged_df['status_prot'] = 0
    merged_df.loc[merged_df['type'] == 'Frameshift', 'status_prot'] = 1
    merged_df.loc[(merged_df['type'] == 'Complex_mutation') & (merged_df['mutant'].apply(lambda x: 'fs' in x)), 'status_prot'] = 1

    return merged_df

def extract_embedding(sequence, tokenizer, model):
    inputs = tokenizer(sequence, return_tensors='pt', max_length=1024, padding=True, truncation=True)
    with torch.no_grad():
        embedding = model(**inputs).last_hidden_state
    return embedding

def extract_and_save_embeddings(prot_data, tokenizer, model, output_file):
    embedding_dict = {}
    for _, row in prot_data.iterrows():
        isoform_id = row['isoform_id']
        if 'mutation_str' in row:
            mut_str = row['mutation_str']
            sequence = row['mut_seq']
            if isoform_id not in embedding_dict:
                embedding_dict[isoform_id] = {}
            embedding_dict[isoform_id][mut_str] = extract_embedding(sequence, tokenizer, model)
        else:
            sequence = row['wt_seq']
            embedding_dict[isoform_id] = extract_embedding(sequence, tokenizer, model)
    
    torch.save(embedding_dict, output_file)
    print(f"Embeddings were saved to {output_file}")

def main():
    setup_environment()
    df = load_full_feature_data()
    
    df_prot = df.loc[(df.type != 'WT') &
                     (df.type != 'Silent_Missense') &
                     (df.type != 'Silent_Nonsense') &
                     (df.status_prot == 0)
                     ]

    wt_prot = df_prot.drop_duplicates(subset=['isoform_id'])[['gene', 'isoform_id', 'wt_seq']]
    mut_prot = df_prot.drop_duplicates(subset=['isoform_id', 'mutation_str'])[['gene', 'mutation_str', 'isoform_id', 'mut_seq']]

    tokenizer = BertTokenizer.from_pretrained('Rostlab/prot_bert')
    model = BertModel.from_pretrained('Rostlab/prot_bert')

    extract_and_save_embeddings(wt_prot, tokenizer, model, 'data/processed/wt_embedding_tensor.pt')
    extract_and_save_embeddings(mut_prot, tokenizer, model, 'data/processed/mut_embedding_tensor.pt')

    print("All embeddings have been extracted and saved.")

if __name__ == "__main__":
    main()