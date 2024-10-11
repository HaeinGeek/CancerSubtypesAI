from Bio import Entrez, SeqIO

def get_entrez_isoforms(gene_name, email):
    """
    NCBI Entrez를 사용하여 주어진 유전자의 단백질 서열을 가져와 반환합니다.

    매개변수:
    - gene_name: 단백질 서열을 가져올 유전자 이름입니다.
    - email: NCBI Entrez API를 사용하기 위한 사용자 이메일입니다.

    반환값:
    - protein_seqs: 접근번호를 키로 하고 서열을 값으로 갖는 딕셔너리입니다.
    """
    Entrez.email = email  # Set email for NCBI Entrez API
    protein_seqs = {}
    print(f"\nFetching protein sequences for gene: {gene_name}")
    search_term = f"{gene_name}[Gene Name] AND Homo sapiens[Organism]"
    try:
        # 단백질 ID를 검색합니다.
        with Entrez.esearch(db="protein", term=search_term, retmax=1000) as handle:
            record = Entrez.read(handle)
        id_list = record.get("IdList", [])

        if not id_list:
            print(f"No protein IDs found for gene: {gene_name}")
            return {}

        # 배치 처리: 여러 단백질 ID를 한 번에 가져옵니다.
        id_str = ','.join(id_list)
        with Entrez.efetch(db="protein", id=id_str, rettype="fasta", retmode="text") as fetch_handle:
            seq_records = list(SeqIO.parse(fetch_handle, "fasta"))

        if not seq_records:
            print(f"No protein sequences found for gene: {gene_name}")
            return {}

        for seq_record in tqdm(seq_records, desc=f"Processing sequences for {gene_name}", leave=False):
            accession = extract_accession(seq_record.id)
            if not accession:
                print(f"Could not extract accession: {seq_record.id}")
                continue
            if accession in protein_seqs:
                continue  # 중복 제거
            sequence = str(seq_record.seq)
            protein_seqs[accession] = sequence

        if protein_seqs:
            print(f"Successfully fetched protein sequences for gene: {gene_name}")
        else:
            print(f"No valid protein sequences found for gene: {gene_name}")

    except Exception as e:
        print(f"Error fetching data for gene {gene_name}: {e}")
        return {}

    return protein_seqs

def extract_accession(seq_id):
    """
    시퀀스 ID에서 접근번호를 추출합니다.

    매개변수:
    - seq_id: 시퀀스 ID 문자열입니다.

    반환값:
    - 접근번호 문자열 또는 None
    """
    if '|' in seq_id:
        parts = seq_id.split('|')
        # 접근번호 패턴에 맞는 부분을 찾습니다.
        for part in parts:
            if part.startswith(('NP_', 'XP_', 'WP_', 'YP_', 'AP_')):
                return part
        # 일반적인 접근번호 패턴이 없을 경우 마지막 부분을 사용합니다.
        return parts[-1] if parts[-1] else None
    else:
        return seq_id if seq_id else None
