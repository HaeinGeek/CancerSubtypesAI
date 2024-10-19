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