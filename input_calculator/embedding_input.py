import pandas as pd
import numpy as np
import json
import h5py
import os
import torch
import collections

class EmbeddingManager:
    def __init__(self, embedding_dir, cache_size=1000):
        self.embedding_dir = embedding_dir
        self.embedding_files = []
        self.cache_size = cache_size  # 최대 캐시 크기
        self.embeddings = collections.OrderedDict()  # 임베딩 캐시 (LRU 방식)
        self.file_index = {}  # 키와 파일 매핑 정보
        self._index_embeddings()  # 임베딩 파일 인덱싱

    def _index_embeddings(self):
        """임베딩 파일들을 인덱싱하여 키와 파일의 매핑을 생성"""
        self.embedding_files = [f for f in os.listdir(self.embedding_dir) if f.endswith('.h5')]
        print(f"Loaded {len(self.embedding_files)} embedding files for indexing")

        for file_name in self.embedding_files:
            file_path = os.path.join(self.embedding_dir, file_name)
            with h5py.File(file_path, 'r') as f:
                for key in f.keys():
                    self.file_index[key] = file_path

    def get_embedding(self, key):
        # 캐시에 임베딩이 있는지 확인
        if key in self.embeddings:
            # LRU 캐시 업데이트
            self.embeddings.move_to_end(key)
            return self.embeddings[key]

        # 캐시에 없으면 파일에서 로드
        if key not in self.file_index:
            print(f"Warning: Key {key} not found in any embedding file")
            return None

        file_path = self.file_index[key]
        with h5py.File(file_path, 'r') as f:
            embedding_data = f[key][:]
            embedding = torch.tensor(embedding_data)

        # 캐시에 임베딩 추가
        self.embeddings[key] = embedding
        # 캐시 크기가 초과되면 가장 오래된 항목 제거
        if len(self.embeddings) > self.cache_size:
            self.embeddings.popitem(last=False)

        return embedding

def create_model_input(df, mutation_df, embedding_dir):
    # EmbeddingManager 인스턴스 생성
    embedding_manager = EmbeddingManager(embedding_dir, cache_size=5000)  # 캐시 크기 조절 가능

    # mutation_df를 사전으로 변환하여 조회 속도 향상
    mutation_type_dict = mutation_df.set_index(['gene', 'mutation_str'])['type'].to_dict()
    mutation_info_dict = mutation_df.set_index(['gene', 'mutation_str']).to_dict('index')

    def process_sample(row):
        result = {}

        # SUBCLASS 열이 있는 경우에만 추가
        if 'SUBCLASS' in df.columns:
            result['subclass'] = row['SUBCLASS']

        gene_columns = df.columns.drop(['ID', 'SUBCLASS'], errors='ignore')
        mutated_genes = [gene for gene in gene_columns if row[gene] != 'WT']
        mutation_strs = [row[gene] for gene in mutated_genes]

        # mutation_types 가져오기
        mutation_types = [
            mutation_type_dict.get((gene, mut), 'Unknown') for gene, mut in zip(mutated_genes, mutation_strs)
        ]

        mutation_stats = {
            'num_mutated_genes': len(mutated_genes),
            'mutations': {mut_str: mutation_strs.count(mut_str) for mut_str in set(mutation_strs)},
            'mutation_type_freq': {type_: mutation_types.count(type_) for type_ in set(mutation_types)}
        }

        # sample_mutations 가져오기
        sample_mutations = [
            mutation_info_dict.get((gene, mut))
            for gene, mut in zip(mutated_genes, mutation_strs)
            if (gene, mut) in mutation_info_dict
        ]
        sample_mutations = mutation_df[(mutation_df['gene'].isin(mutated_genes))&(mutation_df['mutation_str'].isin(mutation_strs))]

        if not sample_mutations.empty:
            aa_properties = ['hydrophobicity', 'polarity', 'mw', 'pI', 'charge']
            aa_change_stats = {}
            for prop in aa_properties:
                cols = [f'{prop}_min', f'{prop}_max', f'{prop}_mean', f'{prop}_std']
                stats = sample_mutations[cols].mean().to_dict()
                aa_change_stats[prop] = stats

            embedding_diffs = []
            for _, mutation_info in sample_mutations.iterrows():
                isoform_id = mutation_info['isoform_id']
                mutation_str = mutation_info['mutation_str']
                status_prot = mutation_info['status_prot']

                if mutation_str != 'WT' and status_prot != 1:
                    wt_key = isoform_id
                    mut_key = f"{isoform_id}_{mutation_str}"

                    wt_emb = embedding_manager.get_embedding(wt_key)
                    mut_emb = embedding_manager.get_embedding(mut_key)

                    if wt_emb is not None and mut_emb is not None:
                        embedding_diffs.append(mut_emb - wt_emb)
                    else:
                        print(f"Warning: Missing embedding for isoform: {isoform_id}, mutation: {mutation_str}, but status_prot != 1")
                        embedding_diffs.append(torch.zeros(3072, dtype=torch.float32))

            if embedding_diffs:
                embedding_diffs = torch.stack(embedding_diffs)
                embedding_stats = {
                    'mean': torch.mean(embedding_diffs, dim=0).numpy(),
                    'max': torch.max(embedding_diffs, dim=0)[0].numpy(),
                    'min': torch.min(embedding_diffs, dim=0)[0].numpy(),
                }

                if embedding_diffs.size(0) > 1:
                    embedding_stats['std'] = torch.std(embedding_diffs, dim=0).numpy()
                else:
                    embedding_stats['std'] = np.zeros_like(embedding_stats['mean'])
            else:
                embedding_stats = None

            additional_stats = {
                'avg_mut_num': sample_mutations['mut_num'].mean(),
                'max_mut_num': sample_mutations['mut_num'].max(),
                'status_ratio': sample_mutations['status'].mean(),
                'status_prot_ratio': sample_mutations['status_prot'].mean()
            }
        else:
            # sample_mutations가 비어있는 경우 기본값 설정
            aa_change_stats = None
            embedding_stats = None
            additional_stats = None

        result.update({
            'id': row.ID,
            'mutation_stats': mutation_stats,
            'aa_change_stats': aa_change_stats,
            'embedding_stats': embedding_stats,
            'additional_stats': additional_stats
        })
        return result

    # DataFrame을 순회하며 결과 생성
    results = df.apply(process_sample, axis=1).tolist()
    return results

def convert_to_serializable(obj):
    """
    객체를 JSON 직렬화 가능한 형태로 변환합니다.
    """
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    elif isinstance(obj, torch.Tensor):
        return obj.tolist()
    elif isinstance(obj, dict):
        return {k: convert_to_serializable(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_to_serializable(v) for v in obj]
    else:
        return obj

def save_model_input(model_input, filename):
    """
    모델 입력 데이터를 JSON 파일로 저장합니다.
    """
    serializable_input = convert_to_serializable(model_input)

    with open(filename, 'w') as f:
        json.dump(serializable_input, f, indent=2)

    print(f"모델 입력 데이터가 {filename}에 저장되었습니다.")

def load_model_input(filename):
    """
    JSON 파일에서 모델 입력 데이터를 불러옵니다.
    """
    with open(filename, 'r') as f:
        loaded_input = json.load(f)

    print(f"{filename}에서 모델 입력 데이터를 불러왔습니다.")
    return loaded_input