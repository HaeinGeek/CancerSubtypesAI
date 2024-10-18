import re
import numpy as np
from collections import Counter


# 아미노산 정보 사전 정의
AMINO_ACID_INFO = {
    'A': 'Alanine',
    'R': 'Arginine',
    'N': 'Asparagine',
    'D': 'Aspartic Acid',
    'C': 'Cysteine',
    'Q': 'Glutamine',
    'E': 'Glutamic Acid',
    'G': 'Glycine',
    'H': 'Histidine',
    'I': 'Isoleucine',
    'L': 'Leucine',
    'K': 'Lysine',
    'M': 'Methionine',
    'F': 'Phenylalanine',
    'P': 'Proline',
    'S': 'Serine',
    'T': 'Threonine',
    'W': 'Tryptophan',
    'Y': 'Tyrosine',
    'V': 'Valine',
    '*': 'Stop Codon'
}


def get_amino_acid_info(amino_acid):
    """
    주어진 아미노산 코드에 해당하는 이름을 반환합니다.
    예외 상황에서는 'Unknown'을 반환합니다.

    Parameters:
    - amino_acid (str): 아미노산 코드 (단일 문자).

    Returns:
    - str: 아미노산 이름 또는 'Unknown'.
    """
    return AMINO_ACID_INFO.get(amino_acid, 'Unknown')



def parse_mutation(mutation):
    """
    변이 문자열을 원래 아미노산, 위치, 돌연변이 아미노산으로 파싱합니다.

    Parameters:
    - mutation (str): 변이 문자열 (예: 'A123T', '*363*', 'P34fs').

    Returns:
    - tuple: (원본 아미노산, 위치 리스트, 변이 아미노산 또는 None).

    Raises:
    - ValueError: 변이 형식이 올바르지 않은 경우.
    """
    if mutation == 'WT' or mutation == 'nan':
        return None, None, None
    
    # single deletion (e.g. 490del)
    position_del_pattern = r'^(\d+)del$'
    match = re.match(position_del_pattern, mutation)

    if match:
        pos = int(match.group(1))
        return None, pos, '-'

    # missense, nonsense, frameshift
    single_pattern = r'^([A-Z*-]+)(\d+)([A-Z*]?|fs\*\d*|fs|del)?$'
    match = re.match(single_pattern, mutation)

    if not match:
        raise ValueError(f"Invalid mutation format: {mutation}")

    orig, pos, mut = match.groups()

    # 'del'이 변이 아미노산인 경우 '-'로 반환
    if mut == 'del':
        mut = '-'

    return orig, int(pos), mut if mut else None




multiple_mutation_patterns = [
    # 연속 위치 변이 (예: 12_14AA>AG)
    (r'^(\d+)_(\d+)([A-Z*]{2,})>([A-Z*]{2,})$', 'range'),
    # 결실 변이 (예: N162_Q172del, 292_293del)
    (r'^([A-Z])(\d+)_([A-Z])(\d+)del$', 'deletion_with_aa'),
    (r'^(\d+)_(\d+)del$', 'deletion_without_aa'),
    # 아미노산 변화가 있는 프레임시프트 변이 (예: D989Tfs)
    (r'^([A-Z])(\d+)([A-Z])fs$', 'complex_fs'),
    # 삽입 변이 (예: C1479_T1480insFND)
    (r'^([A-Z])(\d+)_([A-Z])(\d+)ins([A-Z]+)$', 'insertion'),
    (r'([A-Z])(\d+)_(\d+)ins([A-Z]+)', 'long_insertion'),
    # 결실+삽입 (예: L12delinsRV, D401_L713delinsG)
    (r'([A-Z])(\d+)delins([A-Z]+)', 'single_delins'),
    (r'([A-Z])(\d+)_([A-Z])(\d+)delins([A-Z]+)', 'delins'),
    ]

def parse_multiple_mutations(mutation_str):
    """
    여러 변이 정보를 파싱하여 원본 아미노산과 변이 아미노산을 리스트로 반환합니다.

    Parameters:
    - mutation_str (str): 여러 변이가 포함된 문자열 (예: 'A123T 1499_1500HL>HL D989Tfs N162_Q172del C1479_T1480insFND').

    Returns:
    - tuple: ([원본 아미노산 리스트], [변이 위치 리스트], [변이 아미노산 리스트]).

    Examples:
    >>> parse_multiple_mutations("A123T 1499_1500HL>HL Q581*")
    (['A', 'H', 'L', 'Q'], [123, 1499, 1500, 581], ['T', 'H', 'L', '*'])

    >>> parse_multiple_mutations("A123T L123K L123K")
    Warning: Mutation 'L123K' appears 2 times. Duplicates will be removed.
    (['A', 'L'], [123, 123], ['T', 'K'])

    >>> parse_multiple_mutations("WT")
    ([], [], [])

    >>> parse_multiple_mutations("P34fs")
    (['P'], [34], ['fs'])

    >>> parse_multiple_mutations("12_14AA>AG")
    (['A', 'A', 'A'], [12, 13, 14], ['A', 'A', 'G'])

    >>> parse_multiple_mutations("292_293del")
    ([], [292, 293], ['-', '-'])

    >>> parse_multiple_mutations("L12delinsRV")
    (['L'], [12], ['-', 'R', 'V'])

    >>> parse_multiple_mutations("N162_Q172del")
    (['N', 'Q'], [162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172], ['-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-'])
    
    >>> parse_multiple_mutations("N162_Q172del")
    (['N', 'Q'], [162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172], ['-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-'])

    >>> parse_multiple_mutations("C1479_T1480insFND")
    (['C', 'T'], [1479, 1480], ['C', 'F', 'N', 'D', 'T'])
    """
    if mutation_str == 'WT' or mutation_str == 'nan':
        return [], [], []

    mutations = mutation_str.split()

    mutation_counts = Counter(mutations)
    for mutation, count in mutation_counts.items():
        if count > 1:
            print(f"Warning: Mutation '{mutation}' appears {count} times. Duplicates will be removed.")

    unique_mutations = list(mutation_counts.keys())

    origins, positions, mutants = [], [], []

    for mutation in unique_mutations:
        try:
            # 먼저 parse_mutation 시도
            orig, pos, mut = parse_mutation(mutation)
            if orig is not None or pos is not None:
                if orig: origins.append(orig)
                if pos: positions.append(pos)
                if mut: mutants.append(mut)
                continue
        except ValueError:
            # parse_mutation에서 처리되지 않은 경우, 다른 패턴 시도
            pass

        parsed = False
        for pattern, mutation_type in multiple_mutation_patterns:
            match = re.match(pattern, mutation)
            if match:
                if mutation_type == 'range':
                    start_pos, end_pos, orig, mut = match.groups()
                    orig_list = list(orig)
                    mut_list = list(mut)
                    pos_list = list(range(int(start_pos), int(end_pos) + 1))
                    if len(orig_list) != len(mut_list) or len(orig_list) != len(pos_list):
                        raise ValueError(f"Mismatch in mutation lengths: {mutation}")
                    origins.extend(orig_list)
                    positions.extend(pos_list)
                    mutants.extend(mut_list)
                elif mutation_type == 'deletion_with_aa':
                    start_aa, start_pos, end_aa, end_pos = match.groups()
                    start_pos, end_pos = int(start_pos), int(end_pos)
                    origins.extend([start_aa, end_aa])
                    positions.extend(list(range(start_pos, end_pos + 1)))
                    mutants.extend(['-'] * (end_pos - start_pos + 1))
                elif mutation_type == 'deletion_without_aa':
                    start_pos, end_pos = map(int, match.groups())
                    positions.extend(list(range(start_pos, end_pos + 1)))
                    mutants.extend(['-'] * (end_pos - start_pos + 1))
                elif mutation_type == 'insertion':
                    start_aa, start_pos, end_aa, end_pos, inserted = match.groups()
                    origins.extend([start_aa, end_aa])
                    positions.extend([int(start_pos), int(end_pos)])
                    mutants.extend([start_aa] + list(inserted) + [end_aa])
                elif mutation_type == 'delins':
                    start_aa, start_pos, end_aa, end_pos, inserted = match.groups()
                    start_pos, end_pos = int(start_pos), int(end_pos)
                    origins.extend([start_aa, end_aa])
                    positions.extend(list(range(start_pos, end_pos + 1)))
                    mutants.extend(['-'] * (end_pos - start_pos + 1) + list(inserted))
                elif mutation_type == 'single_delins':
                    orig_aa, pos, inserted = match.groups()
                    origins.append(orig_aa)
                    positions.append(int(pos))
                    mutants.extend(['-'] + list(inserted))
                elif mutation_type == 'long_insertion':
                    start_aa, start_pos, end_pos, inserted = match.groups()
                    origins.extend([start_aa, start_aa])
                    positions.extend([int(start_pos), int(end_pos)])
                    mutants.extend([start_aa] + list(inserted) + [start_aa])
                parsed = True
                break
        
        if not parsed:
            raise ValueError(f"Invalid mutation format: {mutation}")

    return origins, positions, mutants


# # 정규 표현식 패턴 컴파일
# patterns = [
#     # Nonsense 변이
#     (re.compile(r'^([A-Z])(\d+)\*$'), 'Nonsense', lambda m: [int(m.group(2))]),
#     # Missense 변이
#     (re.compile(r'^([A-Z])(\d+)([A-Z])$'), 'Missense', lambda m: [int(m.group(2))] if m.group(1) != m.group(3) else ('Silent_Missense', [int(m.group(2))])),
#     # Silent Nonsense 변이
#     (re.compile(r'^\*(\d+)\*$'), 'Silent_Nonsense', lambda m: [int(m.group(1))]),
#     # Frameshift 변이
#     (re.compile(r'^([A-Z\-]?)(\d+)fs(.*)$'), 'Frameshift', lambda m: [int(m.group(2))]),
#     # 두 개의 연속된 아미노산 변이
#     (re.compile(r'^(\d+)_(\d+)([A-Z]{2})>([A-Z\*]{2})$'), None,
#         lambda m: ('Complex_mutation', [int(m.group(1)), int(m.group(2))]) if '*' in m.group(4) and m.group(4) != '**' else
#                   ('Missense', [int(m.group(1)), int(m.group(2))]) if '*' not in m.group(4) else
#                   ('Nonsense', [int(m.group(1)), int(m.group(2))])),
#     # Deletion 변이
#     (re.compile(r'^([A-Z]?)(\d+)del$'), 'Deletion', lambda m: [int(m.group(2))]),
# ]



# def classify_single_mutation(mutation):
#     """단일 변이를 분석하고 변이 유형 및 위치를 반환하는 함수"""

#     if mutation == 'WT':
#         return ('WT', np.nan)
    
#     for pattern, mutation_type, handler in patterns:
#         m = pattern.match(mutation)
#         if m:
#             result = handler(m)
#             # handler가 튜플을 반환하는 경우 (예: 변이 유형이 변경된 경우)
#             if isinstance(result, tuple):
#                 return result
#             return (mutation_type, result)
    
#     # 'fs'를 포함하지만 위의 패턴에 매칭되지 않는 경우
#     if 'fs' in mutation:
#         return ('Frameshift', np.nan)
    
#     # 그 외의 경우
#     return ('Unknown', np.nan)



# def classify_and_aggregate_mutations(mutations):
#     """
#     여러 변이를 처리하고 변이 유형을 분류한 후, 전체 변이 유형을 결정합니다.

#     Parameters:
#     - mutations (str): 여러 변이가 포함된 문자열 (예: 'A123T Q581*').

#     Returns:
#     - str: 전체 변이 유형.
#     """
#     mutation_list = mutations.split()
    
#     # 중복 확인 및 경고 메시지 출력
#     mutation_counts = Counter(mutation_list)
#     for mutation, count in mutation_counts.items():
#         if count > 1:
#             print(f"Warning: Mutation '{mutation}' appears {count} times. Duplicates will be removed.")

#     # 중복 제거
#     unique_mutations = list(mutation_counts.keys())
    
#     # 변이 유형 분류
#     encoded_mutations = [classify_single_mutation(m) for m in unique_mutations]
    
#     # 전체 변이 유형 결정
#     types = [t for t, _ in encoded_mutations]
#     non_wt_types = [t for t in types if t != 'WT']

#     if len(set(types)) == 1:
#         mutation_type = types[0]
#     elif len(non_wt_types) == 0:
#         mutation_type = 'WT'
#     elif len(set(non_wt_types)) > 1:
#         mutation_type = 'Complex_mutation'

#     return mutation_type

import re
import numpy as np
from collections import Counter

# Updated patterns list with the new mutation patterns
patterns = [
    # Nonsense mutation (e.g., 'Q581*')
    (re.compile(r'^([A-Z])(\d+)\*$'), 'Nonsense', lambda m: [int(m.group(2))]),

    # Missense mutation (e.g., 'A123T')
    (re.compile(r'^([A-Z])(\d+)([A-Z])$'), 'Missense',
     lambda m: [int(m.group(2))] if m.group(1) != m.group(3) else ('Silent_Missense', [int(m.group(2))])),

    # Silent Nonsense mutation (e.g., '*123*')
    (re.compile(r'^\*(\d+)\*$'), 'Silent_Nonsense', lambda m: [int(m.group(1))]),

    # Frameshift mutation (e.g., 'A123fs', '-123fs')
    (re.compile(r'^([A-Z\-]?)(\d+)fs(.*)$'), 'Frameshift', lambda m: [int(m.group(2))]),

    # Deletion at a single position (e.g., 'A123del', '123del')
    (re.compile(r'^([A-Z]?)(\d+)del$'), 'Deletion', lambda m: [int(m.group(2))]),

    # Deletion over a range with amino acids (e.g., 'N162_Q172del')
    (re.compile(r'^([A-Z])(\d+)_([A-Z])(\d+)del$'), 'Deletion',
     lambda m: list(range(int(m.group(2)), int(m.group(4)) + 1))),

    # Deletion over a range without amino acids (e.g., '292_293del')
    (re.compile(r'^(\d+)_(\d+)del$'), 'Deletion',
     lambda m: list(range(int(m.group(1)), int(m.group(2)) + 1))),

    # Insertion mutation (e.g., 'C1479_T1480insFND')
    (re.compile(r'^([A-Z])(\d+)_([A-Z])(\d+)ins([A-Z]+)$'), 'Insertion',
     lambda m: list(range(int(m.group(2)), int(m.group(4)) + 1))),

    # Insertion mutation without second amino acid (e.g., 'A123_124insQ')
    (re.compile(r'^([A-Z])(\d+)_(\d+)ins([A-Z]+)$'), 'Insertion',
     lambda m: [int(m.group(2)), int(m.group(3))]),

    # Deletion-insertion at single position (e.g., 'L12delinsRV')
    (re.compile(r'^([A-Z])(\d+)delins([A-Z]+)$'), 'Delins', lambda m: [int(m.group(2))]),

    # Deletion-insertion over range (e.g., 'D401_L713delinsG')
    (re.compile(r'^([A-Z])(\d+)_([A-Z])(\d+)delins([A-Z]+)$'), 'Delins',
     lambda m: list(range(int(m.group(2)), int(m.group(4)) + 1))),

    # Multiple amino acid substitution (e.g., '12_14AA>AG')
    (re.compile(r'^(\d+)_(\d+)([A-Z\*]+)>([A-Z\*]+)$'), None,
     lambda m: ('Complex_mutation', list(range(int(m.group(1)), int(m.group(2)) + 1)))
     if '*' in m.group(4) and m.group(4) != '*' * len(m.group(4))
     else ('Missense', list(range(int(m.group(1)), int(m.group(2)) + 1)))
     if '*' not in m.group(4)
     else ('Nonsense', list(range(int(m.group(1)), int(m.group(2)) + 1)))),

    # Complex frameshift with amino acid change (e.g., 'D989Tfs')
    (re.compile(r'^([A-Z])(\d+)([A-Z])fs$'), 'Frameshift', lambda m: [int(m.group(2))]),
]

def classify_single_mutation(mutation):
    """Analyze a single mutation and return its type and positions."""

    if mutation == 'WT':
        return ('WT', np.nan)
    
    for pattern, mutation_type, handler in patterns:
        m = pattern.match(mutation)
        if m:
            result = handler(m)
            # If handler returns a tuple (e.g., when mutation type is changed)
            if isinstance(result, tuple):
                return result
            return (mutation_type, result)
    
    # If 'fs' is in mutation but no pattern matched
    if 'fs' in mutation:
        return ('Frameshift', np.nan)
    
    # Otherwise
    return ('Unknown', np.nan)

def classify_and_aggregate_mutations(mutations):
    """
    Process multiple mutations, classify them, and determine the overall mutation type.

    Parameters:
    - mutations (str): A string containing multiple mutations (e.g., 'A123T Q581*').

    Returns:
    - str: The overall mutation type.
    """
    mutation_list = mutations.split()
    
    # Check for duplicates and print a warning
    mutation_counts = Counter(mutation_list)
    for mutation, count in mutation_counts.items():
        if count > 1:
            print(f"Warning: Mutation '{mutation}' appears {count} times. Duplicates will be removed.")
    
    # Remove duplicates
    unique_mutations = list(mutation_counts.keys())
    
    # Classify mutations
    encoded_mutations = [classify_single_mutation(m) for m in unique_mutations]
    
    # Determine overall mutation type
    types = [t for t, _ in encoded_mutations]
    non_wt_types = [t for t in types if t != 'WT']

    # If all mutations are of the same non-'WT' type
    if len(set(non_wt_types)) == 1:
        mutation_type = non_wt_types[0]
    elif len(non_wt_types) == 0:
        mutation_type = 'WT'
    else:
        mutation_type = 'Complex_mutation'

    return mutation_type


def calculate_amino_acid_diff(origin, mutant, amino_acid_features):
    """아미노산 특성 차이를 계산하는 함수"""
    if origin not in amino_acid_features.index or mutant not in amino_acid_features.index:
        return None
    
    # 계산에 사용할 데이터만 선택
    numeric_columns = ['hydrophobicity', 'polarity', 'mw', 'pI', 'charge']
    
    diff = amino_acid_features.loc[mutant, numeric_columns] - amino_acid_features.loc[origin, numeric_columns]
    return diff



def process_mutation_features(row, amino_acid_features):
    """
    변이의 아미노산 특성을 처리하고, 변화를 누적하며 상태를 결정하는 함수
    
    Parameters:
    row (pd.Series): 변이 정보를 포함하는 데이터프레임의 행
    amino_acid_features (pd.DataFrame): 아미노산 특성 정보 (index: 단일문자 아미노산명)

    Returns:
    dict: 처리된 아미노산 특성 변화와 결정된 상태
    """
    mutation_type = row['type']
    origins = row['origin']
    mutants = row['mutant']
    
    feature_changes = {
        'status': 0,
        'hydrophobicity': [],
        'polarity': [],
        'mw': [],
        'pI': [],
        'charge': []
    }
    
    if mutation_type == 'WT':
        return feature_changes
    
    if isinstance(origins, list) and isinstance(mutants, list):
        for origin, mutant in zip(origins, mutants):
            diff = calculate_amino_acid_diff(origin, mutant, amino_acid_features)
            if diff is None:
                feature_changes['status'] = 1
                return feature_changes
            for key in diff.index:
                feature_changes[key].append(diff[key])
    
    elif isinstance(origins, str) and isinstance(mutants, str):
        diff = calculate_amino_acid_diff(origins, mutants, amino_acid_features)
        if diff is None:
            feature_changes['status'] = 1
            return feature_changes
        for key in diff.index:
            feature_changes[key].append(diff[key])
    
    return feature_changes
