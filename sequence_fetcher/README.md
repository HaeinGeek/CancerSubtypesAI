# CancerSubtypesAI Sequence Fetcher

The Sequence Fetcher module of CancerSubtypesAI is a tool designed to retrieve protein sequence information from various biological databases. This module supports databases such as UniProt, NCBI, PDB, and Ensembl.

## Features

- Retrieve protein sequences from multiple databases (UniProt, NCBI, PDB, Ensembl)
- Extract isoform sequence information based on gene names
- Save search results to CSV files

## Dependencies

- biopython
- pandas
- requests
- tqdm

Ensure these libraries are installed before using the Sequence Processor module.

## Usage

1. Initialize the `SequenceFetcher` class:

```python
from sequence_fetcher import SequenceFetcher

fetcher = SequenceFetcher()
```

2. Search for isoform sequences of genes from a specific database:

```python
genes = ['TP53', 'BRCA1', 'EGFR']
db_name = 'uniprot'  # Choose from 'uniprot', 'ncbi', 'pdb', 'ensembl'

# Example using UniProt database
results = fetcher.process_isoforms(db_name, genes)
```

3. When using the NCBI database, an email address is required:

```python
results = fetcher.process_isoforms('ncbi', genes, email='your.email@example.com')
```

## Supported Databases

- UniProt (`uniprot`)
- NCBI (`ncbi`)
- PDB (`pdb`)
- Ensembl (`ensembl`)

## API Documentation

For detailed information on the APIs used, please refer to the following documentation:

- UniProt: [https://www.uniprot.org/help/programmatic_access](https://www.uniprot.org/help/programmatic_access)
- NCBI: [https://www.ncbi.nlm.nih.gov/home/develop/api/](https://www.ncbi.nlm.nih.gov/home/develop/api/)
- PDB: [https://data.rcsb.org/#data-api](https://data.rcsb.org/#data-api)
- Ensembl: [http://asia.ensembl.org/info/docs/index.html](http://asia.ensembl.org/info/docs/index.html)

## Notes

- An email address is required when using the NCBI API.
- When processing large amounts of data, please check the usage policies of each database.
- Internet connection is required, and processing time may vary depending on the number of genes being searched and network conditions.