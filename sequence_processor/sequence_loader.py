import ast
import pandas as pd

def load_protein_dict(db_name):
    """
    주어진 데이터베이스의 단백질 서열 파일을 불러와 딕셔너리로 변환합니다.

    Parameters:
    - db_name (str): 처리할 데이터베이스 이름 ('uniprot', 'ncbi', 'pdb', 'ensembl').

    Returns:
    - dict: 유전자 이름을 키로, 단백질 서열을 값으로 갖는 딕셔너리.
    """
    try:
        # 1. 단백질 서열 DataFrame 불러오기
        filename = f'data/processed/protein_sequences_{db_name}.csv'
        protein_df = pd.read_csv(filename)

        # 2. 문자열을 딕셔너리로 변환
        protein_df['sequence'] = protein_df['sequence'].apply(
            lambda x: ast.literal_eval(x) if isinstance(x, str) else x
        )

        # 3. 딕셔너리로 복원 (유전자 이름을 key로, 서열을 value로 저장)
        return protein_df.set_index('gene')['sequence'].to_dict()

    except FileNotFoundError:
        print(f"파일 '{filename}'을 찾을 수 없습니다. 해당 DB는 건너뜁니다.")
        return {}
    except Exception as e:
        print(f"'{db_name}' 처리 중 오류가 발생했습니다: {e}")
        return {}

class ProteinDatabase:
    """
    단백질 서열을 데이터베이스별로 동적으로 로드하는 클래스.
    필요할 때만 데이터베이스를 메모리에 로드하여 메모리 사용을 최적화합니다.
    """
    def __init__(self, db_names):
        self.db_names = db_names
        self.loaded_dicts = {}  # 이미 로드된 데이터베이스를 저장

    def get_protein_dict(self, db_name):
        """해당 데이터베이스의 단백질 서열을 동적으로 로드하거나 반환합니다."""
        if db_name not in self.db_names:
            raise ValueError(f"'{db_name}'은 지원되지 않는 데이터베이스입니다.")

        if db_name not in self.loaded_dicts:
            print(f"Loading protein dictionary for '{db_name}'...")
            self.loaded_dicts[db_name] = load_protein_dict(db_name)

        return self.loaded_dicts[db_name]