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
    if mutation == 'WT':
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



def parse_multiple_mutations(mutation_str):
    """
    여러 변이 정보를 파싱하여 원본 아미노산과 변이 아미노산을 리스트로 반환합니다.

    Parameters:
    - mutation_str (str): 여러 변이가 포함된 문자열 (예: 'A123T 1499_1500HL>HL').

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
    """
    # 'WT' 처리
    if mutation_str == 'WT':
        return [], [], []

    # 다중 변이 처리 (공백 기준으로 분리)
    mutations = mutation_str.split()

    # 중복 확인 및 경고 메시지 출력
    mutation_counts = Counter(mutations)
    for mutation, count in mutation_counts.items():
        if count > 1:
            print(f"Warning: Mutation '{mutation}' appears {count} times. Duplicates will be removed.")

    # 중복 제거
    unique_mutations = list(mutation_counts.keys())

    origins, positions, mutants = [], [], []

    for mutation in unique_mutations:
        # 연속 위치 변이 처리 (예: '12_14AA>AG')
        if '_' in mutation and '>' in mutation:
            range_pattern = r'^(\d+)_(\d+)([A-Z*]{2,})>([A-Z*]{2,})$'
            match = re.match(range_pattern, mutation)

            if not match:
                raise ValueError(f"Invalid mutation format: {mutation}")

            start_pos, end_pos, orig, mut = match.groups()

            # 연속된 위치를 각각 리스트로 변환
            orig_list = list(orig)
            mut_list = list(mut)
            pos_list = list(range(int(start_pos), int(end_pos) + 1))
            
            # 변이 길이와 위치 길이가 일치하지 않으면 오류 발생
            if len(orig_list) != len(mut_list) or len(orig_list) != len(pos_list):
                raise ValueError(f"Mismatch in mutation lengths: {mutation}")

            # 결과 리스트에 추가
            origins.extend(orig_list)
            positions.extend(pos_list)
            mutants.extend(mut_list)

        else:
            # 단일 변이에 대해 parse_mutation 함수 적용
            orig, pos, mut = parse_mutation(mutation)
            # 유효한 변이만 리스트에 추가
            if pos and mut:
                origins.append(orig)
                positions.append(pos)
                mutants.append(mut)

    return origins, positions, mutants



# 정규 표현식 패턴 컴파일
patterns = [
    # Nonsense 변이
    (re.compile(r'^([A-Z])(\d+)\*$'), 'Nonsense', lambda m: [int(m.group(2))]),
    # Missense 변이
    (re.compile(r'^([A-Z])(\d+)([A-Z])$'), 'Missense', lambda m: [int(m.group(2))] if m.group(1) != m.group(3) else ('Silent_Missense', [int(m.group(2))])),
    # Silent Nonsense 변이
    (re.compile(r'^\*(\d+)\*$'), 'Silent_Nonsense', lambda m: [int(m.group(1))]),
    # Frameshift 변이 (삽입)
    (re.compile(r'^([A-Z]+)(\d+)fs(.*)$'), 'Frameshift_insertion', lambda m: [int(m.group(2))]),
    # Frameshift 변이 (단일 아미노산 또는 '-'로 시작)
    (re.compile(r'^([A-Z\-]?)(\d+)fs(.*)$'), None, lambda m: ('Frameshift_deletion', [int(m.group(2))]) if m.group(1) == '-' else ('Frameshift_insertion', [int(m.group(2))])),
    # 두 개의 연속된 아미노산 변이
    (re.compile(r'^(\d+)_(\d+)([A-Z]{2})>([A-Z]{2})$'), None, lambda m: ('Silent_Multiple', [int(m.group(1)), int(m.group(2))]) if m.group(3) == m.group(4) else ('Multiple_Missense', [int(m.group(1)), int(m.group(2))])),
]



def classify_single_mutation(mutation):
    """단일 변이를 분석하고 변이 유형 및 위치를 반환하는 함수"""

    if mutation == 'WT':
        return ('WT', np.nan)
    
    for pattern, mutation_type, handler in patterns:
        m = pattern.match(mutation)
        if m:
            result = handler(m)
            # handler가 튜플을 반환하는 경우 (예: 변이 유형이 변경된 경우)
            if isinstance(result, tuple):
                return result
            return (mutation_type, result)
    
    # 'fs'를 포함하지만 위의 패턴에 매칭되지 않는 경우
    if 'fs' in mutation:
        return ('Frameshift', np.nan)
    
    # 그 외의 경우
    return ('Unknown', np.nan)



def classify_and_aggregate_mutations(mutations):
    """
    여러 변이를 처리하고 변이 유형을 분류한 후, 전체 변이 유형을 결정합니다.

    Parameters:
    - mutations (str): 여러 변이가 포함된 문자열 (예: 'A123T Q581*').

    Returns:
    - str: 전체 변이 유형.
    """
    mutation_list = mutations.split()
    
    # 중복 확인 및 경고 메시지 출력
    mutation_counts = Counter(mutation_list)
    for mutation, count in mutation_counts.items():
        if count > 1:
            print(f"Warning: Mutation '{mutation}' appears {count} times. Duplicates will be removed.")

    # 중복 제거
    unique_mutations = list(mutation_counts.keys())
    
    # 변이 유형 분류
    encoded_mutations = [classify_single_mutation(m) for m in unique_mutations]
    
    # 전체 변이 유형 결정
    types = [t for t, _ in encoded_mutations]
    non_wt_types = [t for t in types if t != 'WT']

    if len(set(types)) == 1:
        mutation_type = types[0]
    elif len(non_wt_types) == 0:
        mutation_type = 'WT'
    elif len(set(non_wt_types)) > 1:
        mutation_type = 'Complex_mutation'

    return mutation_type