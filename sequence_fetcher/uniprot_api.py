import requests

class UniProtAPI:
    """UniProt API를 사용해 단백질 서열을 가져오는 클래스."""
    
    def get_uniprot_isoforms(self, gene_name):
        """
        UniProt에서 주어진 유전자에 대해 사람의 모든 isoform 아미노산 서열을 가져옵니다.
    
        매개변수:
        - gene_name: 유전자 이름 (예: 'TP53')
    
        반환값:
        - isoform_sequences: isoform ID를 키로 하고 서열을 값으로 갖는 딕셔너리
        """
        # UniProt에서 해당 유전자의 모든 isoform 정보를 가져옵니다.
        base_url = "https://rest.uniprot.org/uniprotkb/stream"
        params = {
            "compressed": "false",
            "format": "fasta",
            "query": f"gene_exact:{gene_name} AND organism_id:9606"
        }
    
        response = requests.get(base_url, params=params, stream=True)
    
        if response.status_code != 200:
            print(f"{gene_name} 유전자의 isoform 정보를 가져오는데 실패했습니다: {response.status_code}")
            return {}
    
        isoform_sequences = {}
        current_header = None
        current_sequence_lines = []
    
        for line in response.iter_lines(decode_unicode=True):
            if line.startswith('>'):
                # 이전 레코드 처리
                if current_header is not None:
                    sequence = ''.join(current_sequence_lines)
                    # 헤더에서 isoform ID 추출
                    fields = current_header.split('|')
                    if len(fields) >= 3:
                        isoform_id = fields[1]
                        isoform_sequences[isoform_id] = sequence
                    else:
                        print(f"헤더에서 isoform ID를 파싱하는데 실패했습니다: {current_header}")
                # 새로운 레코드 시작
                current_header = line
                current_sequence_lines = []
            else:
                current_sequence_lines.append(line.strip())
    
        # 마지막 레코드 처리
        if current_header is not None:
            sequence = ''.join(current_sequence_lines)
            fields = current_header.split('|')
            if len(fields) >= 3:
                isoform_id = fields[1]
                isoform_sequences[isoform_id] = sequence
            else:
                print(f"헤더에서 isoform ID를 파싱하는데 실패했습니다: {current_header}")
    
        return isoform_sequences
