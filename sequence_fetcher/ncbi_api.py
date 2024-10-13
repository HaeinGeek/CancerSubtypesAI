from Bio import Entrez, SeqIO

class NCBIAPI:
    """NCBI API를 사용해 단백질 서열을 가져오는 클래스."""

    def get_ncbi_isoforms(self, gene_name, email):
        """
        NCBI Entrez를 사용하여 주어진 사람 유전자에 대한 모든 isoform 아미노산 서열을 가져옵니다.

        매개변수:
        - gene_name (str): 검색할 유전자 이름.
        - email (str): NCBI API 호출 시 사용할 이메일 주소.
        """
        Entrez.email = email  # NCBI Entrez API에 사용할 이메일 설정
        protein_seqs = {}
        search_term = f"{gene_name}[Gene Name] AND Homo sapiens[Organism]"

        try:
            # 단백질 ID를 검색합니다.
            with Entrez.esearch(db="protein", term=search_term, retmax=1000) as handle:
                record = Entrez.read(handle)
            id_list = record.get("IdList", [])

            if not id_list:
                print(f"{gene_name} 유전자에 해당하는 단백질 ID를 찾을 수 없습니다.")
                return {}

            # 여러 단백질 ID를 한 번에 가져옵니다.
            id_str = ','.join(id_list)
            with Entrez.efetch(db="protein", id=id_str, rettype="fasta", retmode="text") as fetch_handle:
                seq_records = list(SeqIO.parse(fetch_handle, "fasta"))

            if not seq_records:
                print(f"{gene_name} 유전자에 해당하는 단백질 서열을 찾을 수 없습니다.")
                return {}

            # 시퀀스 처리
            for seq_record in seq_records:
                accession = self.extract_accession(seq_record.id)  
                if not accession:
                    print(f"접근 번호를 추출할 수 없습니다: {seq_record.id}")
                    continue
                if accession in protein_seqs:
                    continue  # 중복 제거
                sequence = str(seq_record.seq)
                protein_seqs[accession] = sequence

            if not protein_seqs:
                print(f"{gene_name} 유전자에 대한 유효한 단백질 서열을 찾을 수 없습니다.")

        except Exception as e:
            print(f"{gene_name} 유전자의 데이터를 가져오는 중 오류가 발생했습니다: {e}")
            return {}

        return protein_seqs

    
    def extract_accession(self, seq_id):
        """
        시퀀스 ID에서 접근 번호를 추출합니다.

        매개변수:
        - seq_id: 시퀀스 ID 문자열입니다.

        반환값:
        - 접근 번호 문자열 또는 None
        """
        if '|' in seq_id:
            parts = seq_id.split('|')
            # 접근 번호 패턴에 맞는 부분을 찾습니다.
            for part in parts:
                if part.startswith(('NP_', 'XP_', 'WP_', 'YP_', 'AP_')):
                    return part
            # 일반적인 접근 번호 패턴이 없을 경우 마지막 부분을 사용합니다.
            return parts[-1] if parts[-1] else None
        else:
            return seq_id if seq_id else None
