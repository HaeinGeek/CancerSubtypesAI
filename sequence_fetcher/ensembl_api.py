import requests
from concurrent.futures import ThreadPoolExecutor

class EnsemblAPI:
    """Ensembl API를 사용해 단백질 서열을 가져오는 클래스."""

    def get_ensembl_isoforms(self, gene_name, species='homo_sapiens', max_workers=5):
        server = "https://rest.ensembl.org"
        headers = {"Content-Type": "application/json"}

        # 1. 유전자 이름으로 유전자 정보 조회
        ext = f"/lookup/symbol/{species}/{gene_name}?expand=1"
        response = requests.get(server + ext, headers=headers)
        if not response.ok:
            print(f"유전자 정보를 가져오는 데 실패했습니다: {response.text}")
            return {}
        gene_data = response.json()
        gene_id = gene_data.get('id')
        if not gene_id:
            print(f"{gene_name}에 해당하는 Ensembl 유전자 ID를 찾을 수 없습니다.")
            return {}

        # 2. 전사체 목록 가져오기
        transcripts = gene_data.get('Transcript')
        if not transcripts:
            print(f"{gene_name} 유전자의 전사체를 찾을 수 없습니다.")
            return {}

        # 3. 병렬로 단백질 서열 가져오기
        isoform_sequences = {}
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = [
                executor.submit(self.fetch_protein_sequence, t)
                for t in transcripts
            ]
            for future in futures:
                result = future.result()
                if result:
                    protein_id, sequence = result
                    isoform_sequences[protein_id] = sequence

        return isoform_sequences

    @staticmethod
    def fetch_protein_sequence(transcript):
        if transcript.get('biotype') != 'protein_coding':
            return None
        translation = transcript.get('Translation')
        if not translation:
            return None

        protein_id = translation.get('id')
        if not protein_id:
            return None

        ext = f"/sequence/id/{protein_id}?type=protein"
        response = requests.get(f"https://rest.ensembl.org{ext}", headers={"Content-Type": "application/json"})
        if response.ok:
            seq_data = response.json()
            return protein_id, seq_data.get('seq')
        else:
            print(f"{protein_id}의 단백질 서열을 가져오는 데 실패했습니다.")
            return None
