# Get Protein Sequences

This repository contains functions to fetch all isoform protein sequences for given genes from UniProt and NCBI Entrez databases. This can assist with mutation-based cancer subtype analysis by retrieving necessary protein sequence data.

## Functions

- **`create_protein_dict(genes, get_isoforms_func)`**
  - Processes a list of genes and retrieves all isoform protein sequences.
  - **Arguments**:
    - `genes`: List of gene names.
    - `get_isoforms_func`: Function to fetch isoform sequences for a given gene.
  - **Returns**:
    - `protein_dict`: Dictionary containing genes and their isoform sequences.
    - `not_found_genes`: List of genes with no isoform sequences found.

- **`get_uniprot_isoforms(gene_name)`**
  - Fetches isoform protein sequences from UniProt.

- **`get_entrez_isoforms(gene_names, email)`**
  - Fetches protein sequences from NCBI Entrez.

- **`extract_accession(seq_id)`**
  - Extracts accession numbers from sequence IDs for Entrez.

## Installation

Ensure the following Python packages are installed:

```bash
pip install requests biopython tqdm
```

## Usage

### Example Usage

```python
from your_module import create_protein_dict, get_uniprot_isoforms

# List of genes to process
genes = ["BRCA1", "TP53", "EGFR"]

# Retrieve protein sequences using the provided function
protein_dict, not_found_genes = create_protein_dict(genes, get_uniprot_isoforms)

# Display the results
print(f"Protein sequences: {protein_dict}")
print(f"Genes not found: {not_found_genes}")
```

## Notes

- You will need a valid email address to use the NCBI Entrez API.
- Use `tqdm` to efficiently track the progress of large datasets.
- Ensure proper error handling when fetching data from external sources to manage network-related issues.
- Consider caching results to avoid repeated API calls for the same gene sequences.
