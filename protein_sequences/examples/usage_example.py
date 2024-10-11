from CancerSubtypesAI.get_protein_sequences import get_uniprot_isoforms, get_entrez_isoforms

# 처리할 유전자 리스트
genes = ['BRCA1', 'TP53', 'EGFR']

# 단백질 서열을 저장할 딕셔너리
protein_sequences = {}

# UniProt에서 서열 가져오기
for gene_name in genes:
    print(f"Processing gene: {gene_name}")
    protein_seqs = get_uniprot_isoforms(gene_name)
    if not protein_seqs:
        print(f"No isoforms found for gene: {gene_name}")
    else:
        for accession, sequence in protein_seqs.items():
            print(f"Accession: {accession}, Sequence Length: {len(sequence)}")
    protein_sequences[gene_name] = protein_seqs

# UniProt에서 실패한 유전자 리스트
failed_genes = [gene for gene, seqs in protein_sequences.items() if not seqs]

# NCBI Entrez에서 서열 가져오기
for gene_name in failed_genes:
    protein_seqs = get_entrez_isoforms(gene_name, email="your_email@example.com")
    if protein_seqs:
        protein_sequences[gene_name] = protein_seqs
        for accession, sequence in protein_seqs.items():
            print(f"Accession: {accession}, Sequence Length: {len(sequence)}")
    else:
        print(f"No sequences found for gene: {gene_name}")
