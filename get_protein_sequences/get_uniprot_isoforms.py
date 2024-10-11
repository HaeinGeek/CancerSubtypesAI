import requests

def get_uniprot_isoforms(gene_name):
    """
    Uniprot에서 해당 유전자의 모든 isoform 정보를 가져와 반환합니다.

    매개변수:
    - gene_name: 유전자 이름 문자열입니다.

    반환값:
    - protein_seqs: 접근번호를 키로 하고 서열을 값으로 갖는 딕셔너리입니다.
    """
    base_url = ("https://rest.uniprot.org/uniprotkb/stream?"
                "compressed=false&format=fasta&query=gene_exact:{}+AND+organism_id:9606")
    url = base_url.format(gene_name)
    response = requests.get(url)
    
    if response.status_code != 200:
        print(f"Failed to retrieve isoforms for gene: {gene_name}")
        return {}
    
    fasta_data = response.text
    protein_seqs = {}
    entries = fasta_data.strip().split('>')[1:]
    
    for entry in entries:
        lines = entry.strip().split('\n')
        header = lines[0]
        sequence = ''.join(lines[1:])
        # 헤더에서 접근번호 추출
        fields = header.split('|')
        if len(fields) >= 3:
            accession = fields[1]
            protein_seqs[accession] = sequence
        else:
            print(f"Failed to parse accession from header: {header}")
    
    return protein_seqs
