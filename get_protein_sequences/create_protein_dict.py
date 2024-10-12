from tqdm import tqdm

def create_protein_dict(genes, get_isoforms_func):
    protein_dict = {}
    not_found_genes = []
    
    for gene in tqdm(genes, desc="Processing genes", unit="gene"):
        isoform_sequences = get_isoforms_func(gene)
        
        if isoform_sequences:
            tqdm.write(f"Processing gene: {gene}")
            protein_dict[gene] = isoform_sequences
            for isoform_id, sequence in isoform_sequences.items():
                tqdm.write(f"Isoform ID: {isoform_id}, Sequence Length: {len(sequence)}")
        else:
            tqdm.write(f"No isoforms found for gene: {gene}")
            not_found_genes.append(gene)
    
    # Summary of Results
    tqdm.write("\nSummary:")
    tqdm.write(f"Total number of genes processed: {len(genes)}")
    tqdm.write(f"Genes with isoforms found: {len(protein_dict)}")
    tqdm.write(f"Genes without isoforms: {len(not_found_genes)}")
    if not_found_genes:
        tqdm.write("Genes without isoforms:")
        print(not_found_genes)
    
    return protein_dict, not_found_genes
