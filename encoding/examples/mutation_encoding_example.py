from CancerSubtypesAI.encoding.mutation_encoder import process_mutations

# 사용 예시 1: 여러 변이가 있을 때 처리하는 방법
mutations = "L951V L851V L801V P34fs Q58*"
mutation_type, positions = process_mutations(mutations)

print("Mutation Type:", mutation_type)  # 결과: Complex_mutation
print("Mutation Positions:", positions)  # 결과: [34, 58, 801, 851, 951]


# 사용 예시 2: 데이터프레임 처리하는 방법
import pandas as pd
# CSV 파일 읽기
df = pd.read_csv('train.csv')  # 변이 정보가 들어 있는 데이터프레임

df_type = df.copy()
df_position = df.copy()

for column in df.columns:
    if column not in ['ID', 'SUBCLASS']:  # ID와 SUBCLASS 열 제외
        # encode_mutations 함수가 변이 종류와 변이 위치 두 개의 값을 반환하므로, 이를 각각 저장
        df_type[column], df_position[column] = zip(*df[column].apply(encode_mutations))
