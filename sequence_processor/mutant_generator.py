import re
import numpy as np
import pandas as pd
from tqdm import tqdm

def select_longest_isoform(isoform_sequences, db_name):
    """
    주어진 데이터베이스 기준으로 가장 긴 단백질 isoform 서열을 선택합니다.

    Parameters:
    - isoform_sequences (dict): 유전자의 isoform ID와 단백질 서열을 매핑한 딕셔너리.
    - db_name (str): 단백질 서열이 속한 데이터베이스 이름 ('uniprot', 'ncbi', 'pdb', 'ensembl').

    Returns:
    - str: 선택된 가장 긴 단백질 isoform 서열.

    Notes:
    - 각 데이터베이스의 형식에 따라 우선 순위가 다릅니다:
      * UniProt: 'P'로 시작하는 Swiss-Prot 서열.
      * NCBI: 'NP_'로 시작하는 RefSeq 서열, 없을 경우 'XP_' 서열.
      * PDB: PDB 식별자 형식 ('4자리 + "_" + 숫자').
      * Ensembl: 최신 버전의 'ENSP'로 시작하는 서열.
    """
    if db_name.lower() == 'uniprot':
        selected_sequences = {k: v for k, v in isoform_sequences.items() if k.startswith('P')}
    elif db_name.lower() == 'ncbi':
        selected_sequences = {k: v for k, v in isoform_sequences.items() if k.startswith('NP_')}
        if not selected_sequences:
            selected_sequences = {k: v for k, v in isoform_sequences.items() if k.startswith('XP_')}
    elif db_name.lower() == 'pdb':
        selected_sequences = {k: v for k, v in isoform_sequences.items() if re.match(r'^[0-9A-Za-z]{4}_\d+$', k)}
    elif db_name.lower() == 'ensembl':
        ensembl_sequences = {k: v for k, v in isoform_sequences.items() if k.startswith('ENSP')}
        if ensembl_sequences:
            latest_version = max(
                ensembl_sequences.items(),
                key=lambda x: int(x[0].split('.')[-1]) if '.' in x[0] else 0
            )
            return latest_version[1]

    if not selected_sequences:
        selected_sequences = isoform_sequences

    return max(selected_sequences.items(), key=lambda x: len(x[1]))[1]

def parse_mutation(mutation):
    """
    돌연변이 문자열을 원래 아미노산, 위치, 돌연변이 아미노산으로 파싱합니다.

    Parameters:
    - mutation (str): 돌연변이 정보가 포함된 문자열 (예: 'A123T', 'Q58*').

    Returns:
    - tuple: (원래 아미노산, 위치(int), 돌연변이 아미노산).

    Raises:
    - ValueError: 유효하지 않은 돌연변이 형식일 경우 예외를 발생시킵니다.

    Example:
    >>> parse_mutation("A123T")
    ('A', 123, 'T')
    """
    pattern = r'^([A-Z])(\d+)([A-Z*]|fs\*\d*|fs)$'
    match = re.match(pattern, mutation)
    if not match:
        raise ValueError(f"Invalid mutation format: {mutation}")

    orig, pos, mut = match.groups()
    return orig, int(pos), mut

def mutate_sequence(wt_sequence, mutation_str):
    """
    단백질의 WT 서열에 돌연변이를 적용합니다.

    Parameters:
    - wt_sequence (str): 변이 적용 전의 WT 단백질 서열.
    - mutation_str (str): 하나 이상의 돌연변이 문자열 (예: 'A123T, Q58*').

    Returns:
    - str: 돌연변이가 적용된 새로운 단백질 서열.

    Raises:
    - ValueError: 돌연변이 위치가 서열 범위를 초과하거나, 원래 아미노산이 일치하지 않는 경우 예외를 발생시킵니다.
    """
    mutations = re.split(r'\s+|,', mutation_str.strip())
    sequence = list(wt_sequence)

    for mutation in mutations:
        if not mutation:
            continue
        orig, pos, mut = parse_mutation(mutation)

        if pos - 1 >= len(sequence):
            raise ValueError(f"Position {pos} is out of range in the sequence.")

        if sequence[pos - 1] != orig:
            raise ValueError(f"Position {pos} in WT sequence does not match the original amino acid {orig}")

        if mut == '*':
            sequence = sequence[:pos - 1]  # 종결 코돈 이후 서열 제거
            sequence.append('*')
            break
        elif mut.startswith('fs'):
            sequence = sequence[:pos - 1]  # 프레임시프트 변이 이후 서열 제거
            break
        else:
            sequence[pos - 1] = mut

    return ''.join(sequence)

def add_mutated_sequences(unique_mutations_df, protein_dict, db_name):
    """
    유전자와 변이 정보를 사용해 변이 또는 WT 서열을 데이터프레임에 추가합니다.

    Parameters:
    - unique_mutations_df (pd.DataFrame): 유전자와 변이 정보가 포함된 데이터프레임.
      * 필수 열: 'gene', 'mutation_str'
    - protein_dict (dict): 유전자명과 해당 유전자 isoform 서열의 딕셔너리.
    - db_name (str): 단백질 서열이 속한 데이터베이스 이름.

    Returns:
    - pd.DataFrame: 'sequence' 열이 추가된 데이터프레임.
    
    Notes:
    - 'sequence' 열이 없으면 NaN 값으로 초기화됩니다.
    - 'sequence'가 NaN인 경우에만 서열이 추가됩니다.
    - 변이 정보를 먼저 처리해 변이 서열을 저장하고,
      변이와 일치하는 isoform의 원본 서열 중 가장 긴 것을 WT 서열로 저장합니다.
    """
    # 'sequence' 열이 없으면 NaN 값으로 초기화
    if 'sequence' not in unique_mutations_df.columns:
        unique_mutations_df['sequence'] = np.nan

    # NaN인 행만 필터링 (이미 처리된 서열은 건너뜀)
    nan_rows = unique_mutations_df[unique_mutations_df['sequence'].isna()]

    print("돌연변이 서열을 생성 중입니다...")
    for row in tqdm(nan_rows.itertuples(index=True), total=len(nan_rows), desc="변이 처리", unit="건"):
        gene = row.gene
        mutation_str = row.mutation_str
        isoform_sequences = protein_dict.get(gene, {})

        if not isoform_sequences:
            continue  # isoform이 없는 경우 넘어감

        # 1. 변이가 있는 경우 처리 (mutation_str != 'WT')
        if mutation_str != 'WT':
            longest_mutated_seq = None  # 가장 긴 변이 서열 추적
            max_mutated_length = 0

            for isoform_id, wt_seq in isoform_sequences.items():
                try:
                    mutated_seq = mutate_sequence(wt_seq, mutation_str)

                    # 가장 긴 변이 서열을 추적
                    if len(mutated_seq) > max_mutated_length:
                        longest_mutated_seq = mutated_seq
                        max_mutated_length = len(mutated_seq)

                except ValueError as e:
                    print(f"Error in gene {gene} with mutation {mutation_str}: {e}")
                    continue

            # 변이 서열이 존재하면 저장하고 다음 행으로 넘어감
            if longest_mutated_seq:
                unique_mutations_df.at[row.Index, 'sequence'] = longest_mutated_seq
                continue  # 다음 행 처리

        # 2. WT 처리 (mutation_str == 'WT')
        # 가장 긴 원본(WT) 서열을 찾기
        longest_wt_seq = None  # 가장 긴 원본 서열 추적
        max_wt_length = 0

        for isoform_id, wt_seq in isoform_sequences.items():
            if len(wt_seq) > max_wt_length:
                longest_wt_seq = wt_seq
                max_wt_length = len(wt_seq)

        # 가장 긴 원본(WT) 서열을 저장
        if longest_wt_seq:
            unique_mutations_df.at[row.Index, 'sequence'] = longest_wt_seq

    return unique_mutations_df