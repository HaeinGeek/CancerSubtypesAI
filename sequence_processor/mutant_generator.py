import re
from collections import Counter
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

def parse_single_mutation(mutation):
    """
    단일 돌연변이 문자열을 원래 아미노산, 위치, 돌연변이 아미노산으로 파싱합니다.

    Parameters:
    - mutation (str): 돌연변이 정보가 포함된 문자열 (예: 'A123T', 'Q58*').

    Returns:
    - tuple: (원래 아미노산, 위치(int), 돌연변이 아미노산).

    Raises:
    - ValueError: 유효하지 않은 돌연변이 형식일 경우 예외를 발생시킵니다.

    Example:
    >>> parse_single_mutation("A123T")
    ('A', 123, 'T')
    """
    patterns = [
        r'^([A-Z*])(\d+)([A-Z*])$',  # 일반적인 아미노산 또는 종결코돈 치환
        r'^([A-Z])(\d+)(del)$',    # 아미노산 삭제
        r'^([A-Z])(\d+)(fs)$',    # 아미노산 시퀀스로 시작하는 프레임시프트
        r'^(-)(\d+)(fs)$',        # 음수 기호로 시작하는 프레임시프트
    ]
    
    for pattern in patterns:
        match = re.match(pattern, mutation)
        if match:
            orig, pos, mut = match.groups()
            return orig, int(pos), mut
    
    raise ValueError(f"Invalid mutation format: {mutation}")

# def parse_multiple_mutation(mutation):
#     """
#     연속 위치 돌연변이 문자열을 파싱합니다.
#     """
#     range_pattern = r'^(\d+)_(\d+)([A-Z]{2,})>([A-Z*]{1,})$'  # 12_13AL>A*
#     range_match = re.match(range_pattern, mutation)
    
#     if range_match:
#         start_pos, end_pos, orig, mut = range_match.groups()
#         return orig, int(start_pos), int(end_pos), mut
    
#     fs_pattern = r'^([A-Z]{2,})(\d+)(fs)$'   # AKL23fs
#     fs_match = re.match(fs_pattern, mutation)

#     if fs_match:
#         orig, end_pos, mut = fs_match.groups()
#         start_pos = len(orig) - 1
#         return orig, int(start_pos), int(end_pos), mut


#     raise ValueError(f"Invalid mutation format: {mutation}")

# def mutate_sequence(wt_sequence, mutation_str):
#     """
#     단백질의 WT 서열에 돌연변이를 적용합니다.

#     Parameters:
#     - wt_sequence (str): 변이 적용 전의 WT 단백질 서열.
#     - mutation_str (str): 하나 이상의 돌연변이 문자열 (예: 'A123T, Q58*').

#     Returns:
#     - str: 돌연변이가 적용된 새로운 단백질 서열.

#     Raises:
#     - ValueError: 돌연변이 위치가 서열 범위를 초과하거나, 원래 아미노산이 일치하지 않는 경우 예외를 발생시킵니다.
#     """
#     # mutations = re.split(r'\s+|,', mutation_str.strip())
#     # mutations = [m.strip() for m in mutations if m.strip()]  # 공백 제거 및 빈 문자열 필터링
#     mutations = mutation_str.split()
    
#     # 중복 확인 및 경고 메시지 출력
#     mutation_counts = Counter(mutations)
#     for mutation, count in mutation_counts.items():
#         if count > 1:
#             # 중복 제거
#             mutations = list(set(mutations))
#             print(f"Warning: Mutation '{mutation}' appears {count} times. Duplicates will be removed.")

#     sequence = list(wt_sequence)
#     for mutation in mutations:
#         if not mutation:
#             continue
#         try:
#             # 단일 변이 처리 시도 
#             orig, pos, mut = parse_single_mutation(mutation)

#             if pos - 1 >= len(sequence):
#                 raise ValueError(f"Position {pos} is out of range in the sequence.")

#             if orig != '-' and sequence[pos - 1] != orig:
#                 raise ValueError(f"Position {pos} in WT sequence does not match the original amino acid {orig}")

#             if mut == '*' or mut == 'fs':
#                 # 전사 시작 위치 이전의 프레임시프트 처리
#                 if orig == '-':
#                     sequence = list(mut)
#                 else:
#                     sequence = sequence[:pos - 1] + [mut[0]]
#                 break  # 종결 코돈 또는 프레임시프트 이후의 변이는 무시

#             else:
#                 sequence[pos - 1] = mut

#         except ValueError:
#             # 단일 위치 변이 처리 실패 시 연속 위치 변이 처리 시도
#             try:
#                 orig, start_pos, end_pos, mut = parse_multiple_mutation(mutation)
                
#                 if end_pos > len(sequence):
#                     raise ValueError(f"Position {end_pos} is out of range in the sequence.")
                
#                 if ''.join(sequence[start_pos-1:end_pos]) != orig:
#                     raise ValueError(f"Positions {start_pos}-{end_pos} in WT sequence do not match the original amino acids {orig}")
                
#                 if mut[0] == '*':
#                     sequence = sequence[:start_pos-1] + ['*']
#                     break  # 종결 코돈 이후의 변이는 무시
#                 elif mut[-1] == '*':
#                     sequence = sequence[:end_pos-1] + ['*']
#                     break  # 종결 코돈 이후의 변이는 무시
#                 else:
#                     sequence[start_pos-1:end_pos] = list(mut)
            
#             except ValueError:
#                 # 두 가지 파싱 방법 모두 실패한 경우
#                 raise ValueError(f"Invalid mutation format: {mutation}")

#     return ''.join(sequence)

# def add_mutated_sequences(unique_mutations_df, protein_dict, db_name):
#     """
#     유전자와 변이 정보를 사용해 변이 또는 WT 서열을 데이터프레임에 추가합니다.

#     Parameters:
#     - unique_mutations_df (pd.DataFrame): 유전자와 변이 정보가 포함된 데이터프레임.
#       * 필수 열: 'gene', 'mutation_str'
#     - protein_dict (dict): 유전자명과 해당 유전자 isoform 서열의 딕셔너리.
#     - db_name (str): 단백질 서열이 속한 데이터베이스 이름.

#     Returns:
#     - pd.DataFrame: 'sequence' 열이 추가된 데이터프레임.
    
#     Notes:
#     - 'isoform_id', 'wt_seq', 'mut_seq' 열이 없으면 NaN 값으로 초기화됩니다.
#     - 'mut_seq' 열 값이 NaN인 경우에만 서열이 추가됩니다.
#     - 변이가 있는 경우 가장 긴 변이 서열과 해당 isoform id, 원본 서열을 저장합니다.
#     - 변이가 없는 경우(WT) 가장 긴 원본(WT) 서열과 해당 isoform id를 저장합니다.
#     """
#     # 'isoform_id', 'wt_seq', 'mut_seq' 열이 없으면 NaN으로 초기화
#     for col in ['isoform_id', 'wt_seq', 'mut_seq']:
#         if col not in unique_mutations_df.columns:
#             unique_mutations_df[col] = np.nan

#     # NaN인 행만 필터링 (이미 처리된 서열은 건너뜀)
#     nan_rows = unique_mutations_df[
#         unique_mutations_df['mut_seq'].isna() | unique_mutations_df['wt_seq'].isna() | unique_mutations_df['isoform_id'].isna()
#     ]

#     print("돌연변이 서열을 생성 중입니다...")
#     for row in tqdm(nan_rows.itertuples(index=True), total=len(nan_rows), desc="변이 처리", unit="건"):
#         gene = row.gene
#         mutation_str = row.mutation_str
#         isoform_sequences = protein_dict.get(gene, {})

#         longest_isoform_id = None
#         longest_wt_seq = None
#         longest_mut_seq = None
#         max_length = 0

#         if not isoform_sequences:
#             continue  # isoform이 없는 경우 넘어감

#         # 1. 변이가 있는 경우 처리 (mutation_str != 'WT')
#         if mutation_str != 'WT':
#             for isoform_id, wt_seq in isoform_sequences.items():
#                 try:
#                     mutated_seq = mutate_sequence(wt_seq, mutation_str)

#                     # 가장 긴 변이 서열을 추적
#                     if len(mutated_seq) > max_length:
#                         longest_isoform_id = isoform_id
#                         longest_wt_seq = wt_seq
#                         longest_mut_seq = mutated_seq
#                         max_length = len(mutated_seq)

#                 except ValueError as e:
#                     print(f"Error in gene {gene} with mutation {mutation_str}: {e}")
#                     continue

#         # 2. WT 처리 (mutation_str == 'WT')
#         # 가장 긴 원본(WT) 서열을 찾기
#         if mutation_str == 'WT':
#             for isoform_id, wt_seq in isoform_sequences.items():
#                 if len(wt_seq) > max_length:
#                     longest_isoform_id = isoform_id
#                     longest_wt_seq = wt_seq
#                     longest_mut_seq = None
#                     max_length = len(wt_seq)

#         # isoform_id, wt_seq, mut_seq를 저장
#         unique_mutations_df.at[row.Index, 'isoform_id'] = longest_isoform_id
#         unique_mutations_df.at[row.Index, 'wt_seq'] = longest_wt_seq
#         unique_mutations_df.at[row.Index, 'mut_seq'] = longest_mut_seq

#     return unique_mutations_df

def parse_single_mutation(mutation):
    """
    단일 돌연변이 문자열을 원래 아미노산, 위치, 변이된 아미노산으로 파싱합니다.

    반환:
    - tuple: (원래 아미노산, 위치 (int), 변이된 아미노산)

    예외:
    - ValueError: 돌연변이 형식이 유효하지 않은 경우 발생합니다.
    """
    patterns = [
        r'^([A-Z\*])(\d+)([A-Z\*])$',        # 일반적인 아미노산 또는 종결 코돈 치환
        r'^([A-Z])(\d+)(del)$',              # 아미노산 결실
        r'^([A-Z\-]?)(\d+)(fs)$',            # 아미노산 또는 '-'로 시작하는 프레임시프트
        r'^([A-Z])(\d+)delins([A-Z]+)$',     # 단일 위치에서의 결실-삽입
    ]
    
    for pattern in patterns:
        match = re.match(pattern, mutation)
        if match:
            groups = match.groups()
            if len(groups) == 3:
                orig, pos, mut = groups
                return orig, int(pos), mut
            else:
                raise ValueError(f"Invalid mutation format: {mutation}")
    
    raise ValueError(f"Invalid mutation format: {mutation}")

def parse_multiple_mutation(mutation):
    """
    범위를 포함하는 돌연변이를 파싱합니다.

    반환:
    - tuple: (원래 아미노산들, 시작 위치 (int), 끝 위치 (int), 변이된 아미노산들)

    예외:
    - ValueError: 돌연변이 형식이 유효하지 않거나 필요한 아미노산 정보가 부족한 경우 발생합니다.
    """
    # 아미노산이 있는 연속 위치 돌연변이 (예: 12_14AA>AG)
    range_pattern = r'^(\d+)_(\d+)([A-Z\*]{2,})>([A-Z\*]+)$'
    match = re.match(range_pattern, mutation)
    if match:
        start_pos, end_pos, orig, mut = match.groups()
        return orig, int(start_pos), int(end_pos), mut

    # 아미노산이 있는 범위 결실 (예: N162_Q172del)
    del_range_aa_pattern = r'^([A-Z])(\d+)_([A-Z])(\d+)del$'
    match = re.match(del_range_aa_pattern, mutation)
    if match:
        orig_start_aa, start_pos, orig_end_aa, end_pos = match.groups()
        orig = orig_start_aa + '...' + orig_end_aa  # 원래 아미노산들의 플레이스홀더
        return orig, int(start_pos), int(end_pos), 'del'

    # 다음 패턴들은 아미노산 정보가 부족하여 처리되지 않습니다
    # 아미노산 없는 범위 결실 (예: 292_293del)
    del_range_pattern = r'^(\d+)_(\d+)del$'
    match = re.match(del_range_pattern, mutation)
    if match:
        # 원래 아미노산 정보가 없으므로 이 돌연변이를 처리할 수 없습니다
        raise ValueError(f"Cannot process mutation without amino acid information: {mutation}")

    # 삽입 돌연변이 (예: C1479_T1480insFND)
    ins_range_pattern = r'^([A-Z])(\d+)_([A-Z])(\d+)ins([A-Z]+)$'
    match = re.match(ins_range_pattern, mutation)
    if match:
        orig_start_aa, start_pos, orig_end_aa, end_pos, inserted_aa = match.groups()
        orig = orig_start_aa + '...' + orig_end_aa
        return orig, int(start_pos), int(end_pos), 'ins' + inserted_aa

    # 두 번째 아미노산 없는 삽입 돌연변이 (예: A123_124insQ)
    ins_simple_pattern = r'^([A-Z])(\d+)_(\d+)ins([A-Z]+)$'
    match = re.match(ins_simple_pattern, mutation)
    if match:
        orig_start_aa, start_pos, end_pos, inserted_aa = match.groups()
        orig = orig_start_aa
        return orig, int(start_pos), int(end_pos), 'ins' + inserted_aa

    # 범위에 걸친 결실-삽입 돌연변이 (예: D401_L713delinsG)
    delins_range_pattern = r'^([A-Z])(\d+)_([A-Z])(\d+)delins([A-Z]+)$'
    match = re.match(delins_range_pattern, mutation)
    if match:
        orig_start_aa, start_pos, orig_end_aa, end_pos, new_aa = match.groups()
        orig = orig_start_aa + '...' + orig_end_aa
        return orig, int(start_pos), int(end_pos), 'delins' + new_aa

    # 아미노산 변화가 있는 복잡한 프레임시프트 돌연변이 (예: D989Tfs)
    complex_fs_pattern = r'^([A-Z])(\d+)([A-Z])fs$'
    match = re.match(complex_fs_pattern, mutation)
    if match:
        orig, pos, mut = match.groups()
        return orig, int(pos), int(pos), mut + 'fs'

    raise ValueError(f"Invalid mutation format: {mutation}")

def mutate_sequence(wt_sequence, mutation_str):
    """
    야생형 단백질 서열에 돌연변이를 적용합니다.

    매개변수:
    - wt_sequence (str): 야생형 단백질 서열
    - mutation_str (str): 돌연변이 문자열, 예: 'A123T Q58*'

    반환:
    - str: 변이된 단백질 서열

    예외:
    - ValueError: 돌연변이 위치가 범위를 벗어나거나 원래 아미노산이 일치하지 않는 경우 발생합니다.
    """
    mutations = mutation_str.split()
    
    # 중복 제거 및 필요시 경고
    mutation_counts = Counter(mutations)
    for mutation, count in mutation_counts.items():
        if count > 1:
            print(f"Warning: Mutation '{mutation}' appears {count} times. Duplicates will be removed.")
    mutations = list(mutation_counts.keys())

    sequence = list(wt_sequence)
    for mutation in mutations:
        if not mutation:
            continue
        try:
            # 단일 돌연변이로 파싱 시도
            orig, pos, mut = parse_single_mutation(mutation)

            if pos - 1 >= len(sequence):
                raise ValueError(f"Position {pos} is out of range in the sequence.")

            if orig != '-' and sequence[pos - 1] != orig and orig != '*':
                raise ValueError(f"Position {pos} in WT sequence does not match the original amino acid {orig}")

            if mut == '*' or mut == 'fs':
                # 종결 코돈 또는 프레임시프트 돌연변이
                sequence = sequence[:pos - 1] + ['*']
                break  # 종결 코돈 또는 프레임시프트 이후의 돌연변이 무시
            elif mut == 'del':
                # 단일 위치 결실
                sequence = sequence[:pos - 1] + sequence[pos:]
            elif mut.startswith('delins'):
                # 단일 위치 결실-삽입
                new_aa = mut[len('delins'):]
                sequence = sequence[:pos - 1] + list(new_aa) + sequence[pos:]
            else:
                # 아미노산 치환
                sequence[pos - 1] = mut

        except ValueError as e_single:
            # 다중 돌연변이로 파싱 시도
            try:
                orig, start_pos, end_pos, mut = parse_multiple_mutation(mutation)
                
                if end_pos > len(sequence):
                    raise ValueError(f"Position {end_pos} is out of range in the sequence.")

                # 서열에서 원래 아미노산 추출
                seq_orig = ''.join(sequence[start_pos - 1:end_pos])
                if orig and '...' not in orig and orig != seq_orig:
                    raise ValueError(f"Positions {start_pos}-{end_pos} in WT sequence do not match the original amino acids {orig}")

                if mut == 'del':
                    # 범위 결실
                    sequence = sequence[:start_pos - 1] + sequence[end_pos:]
                elif mut.startswith('ins'):
                    # 삽입
                    inserted_aa = mut[len('ins'):]
                    sequence = sequence[:end_pos] + list(inserted_aa) + sequence[end_pos:]
                elif mut.startswith('delins'):
                    # 범위 결실-삽입
                    new_aa = mut[len('delins'):]
                    sequence = sequence[:start_pos - 1] + list(new_aa) + sequence[end_pos:]
                elif mut.endswith('fs'):
                    # 프레임시프트 돌연변이
                    sequence = sequence[:start_pos - 1] + [mut[:-2] + 'fs']
                    break
                else:
                    # 범위 치환
                    sequence[start_pos - 1:end_pos] = list(mut)
            except ValueError as e_multiple:
                # 두 파싱 방법 모두 실패하거나 아미노산 정보 부족으로 돌연변이를 처리할 수 없음
                print(f"Skipping mutation '{mutation}': {e_multiple}")
                continue  # 이 돌연변이를 건너뛰고 다음으로 진행

    # 종결 코돈 '*'에서 서열 잘라내기
    if '*' in sequence:
        sequence = sequence[:sequence.index('*')]

    # 프레임시프트 마커 'fs' 제거 (프레임시프트 효과를 시뮬레이션하려면 더 복잡한 로직이 필요할 수 있음)
    sequence = [aa for aa in sequence if not aa.endswith('fs')]

    return ''.join(sequence)

def add_mutated_sequences(unique_mutations_df, protein_dict, db_name):
    """
    유전자와 돌연변이 정보를 사용하여 데이터프레임에 변이 또는 WT 서열을 추가합니다.

    매개변수:
    - unique_mutations_df (pd.DataFrame): 'gene'과 'mutation_str' 열을 포함하는 DataFrame
    - protein_dict (dict): 유전자 이름을 아이소폼 서열에 매핑하는 딕셔너리
    - db_name (str): 단백질 서열이 속한 데이터베이스의 이름

    반환:
    - pd.DataFrame: 'isoform_id', 'wt_seq', 'mut_seq' 열이 업데이트된 DataFrame
    """
    # 열이 존재하지 않으면 초기화
    for col in ['isoform_id', 'wt_seq', 'mut_seq']:
        if col not in unique_mutations_df.columns:
            unique_mutations_df[col] = np.nan

    # 처리가 필요한 행 필터링
    nan_rows = unique_mutations_df[
        unique_mutations_df['mut_seq'].isna() | unique_mutations_df['wt_seq'].isna() | unique_mutations_df['isoform_id'].isna()
    ]

    print("Generating mutated sequences...")
    for row in tqdm(nan_rows.itertuples(index=True), total=len(nan_rows), desc="Processing mutations", unit="entries"):
        gene = row.gene
        mutation_str = row.mutation_str
        isoform_sequences = protein_dict.get(gene, {})

        longest_isoform_id = None
        longest_wt_seq = None
        longest_mut_seq = None
        max_length = 0

        if not isoform_sequences:
            continue  # 이용 가능한 아이소폼이 없으면 건너뛰기

        # 돌연변이 처리
        if mutation_str != 'WT':
            mutation_processed = False
            for isoform_id, wt_seq in isoform_sequences.items():
                try:
                    mutated_seq = mutate_sequence(wt_seq, mutation_str)

                    # 가장 긴 변이 서열 추적
                    if len(mutated_seq) > max_length:
                        longest_isoform_id = isoform_id
                        longest_wt_seq = wt_seq
                        longest_mut_seq = mutated_seq
                        max_length = len(mutated_seq)
                        mutation_processed = True

                except ValueError as e:
                    print(f"Error in gene {gene} with mutation '{mutation_str}': {e}")
                    continue

            if not mutation_processed:
                print(f"Could not process mutation '{mutation_str}' for gene '{gene}'.")
                continue  # 처리할 수 없는 돌연변이는 건너뛰기

        # WT 서열 처리
        if mutation_str == 'WT':
            for isoform_id, wt_seq in isoform_sequences.items():
                if len(wt_seq) > max_length:
                    longest_isoform_id = isoform_id
                    longest_wt_seq = wt_seq
                    longest_mut_seq = None
                    max_length = len(wt_seq)

        # isoform_id, wt_seq, mut_seq 저장
        unique_mutations_df.at[row.Index, 'isoform_id'] = longest_isoform_id
        unique_mutations_df.at[row.Index, 'wt_seq'] = longest_wt_seq
        unique_mutations_df.at[row.Index, 'mut_seq'] = longest_mut_seq or longest_wt_seq

    return unique_mutations_df