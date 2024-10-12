import re
import numpy as np
from collections import Counter

# 정규 표현식을 미리 컴파일하여 재사용성 및 효율성 향상
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
    # Wild Type
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

def process_mutations(mutations):
    """여러 변이를 처리하고 변이 종류와 위치를 추출하는 함수"""
    mutation_list = mutations.split()
    
    # 중복 확인 및 경고 메시지 생성
    mutation_counts = Counter(mutation_list)
    for mutation, count in mutation_counts.items():
        if count > 1:
            print(f"Warning: Mutation '{mutation}' appears {count} times. Duplicates will be removed.")
    
    # 중복 제거
    unique_mutations = list(mutation_counts.keys())
    
    # 각 변이를 처리하고 변이 종류와 위치를 추출
    encoded_mutations = [classify_single_mutation(m) for m in unique_mutations]
    
    # 전체 변이 유형 결정
    types = [t for t, _ in encoded_mutations]
    non_wt_types = [t for t in types if t != 'WT']

    if len(set(types)) == 1:
        mutation_type = types[0]
    elif len(set(non_wt_types)) == 1:
        mutation_type = f'Multiple_{non_wt_types[0]}'
    elif len(non_wt_types) == 0:
        mutation_type = 'WT'
    else:
        mutation_type = 'Complex_mutation'
    
    # 모든 위치를 수집하고 정렬 (np.nan이 아닌 경우만)
    positions = []
    for _, pos_list in encoded_mutations:
        if isinstance(pos_list, list):
            positions.extend(pos_list)
    
    # 빈 리스트로 처리
    if not positions:
        positions = []
    else:
        positions = sorted(set(positions))
    
    return mutation_type, positions
