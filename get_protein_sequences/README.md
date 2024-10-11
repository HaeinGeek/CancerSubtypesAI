# Get Protein Sequences

This repository contains functions to fetch all isoform protein sequences for given genes from UniProt and NCBI Entrez databases.

## Functions

- `get_uniprot_isoforms(gene_name)`
  - Fetches isoform protein sequences from UniProt.
- `get_entrez_isoforms(gene_names, protein_sequences, email)`
  - Fetches protein sequences from NCBI Entrez.
- `extract_accession(seq_id)`
  - Extracts accession numbers from sequence IDs for Entrez.

## Usage

### Installation

No special installation is required. However, you need to have the following Python packages installed:

- `requests`
- `biopython`

You can install them using `pip`:

```bash
pip install requests biopython
