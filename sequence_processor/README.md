# CancerSubtypesAI Sequence Processor

The Sequence Processor module is a part of the CancerSubtypesAI project, designed to handle protein sequence data, generate mutant sequences, and manage protein databases efficiently. This module consists of two main components: the Mutant Generator and the Sequence Loader.

## Features

- Parse and validate mutation strings
- Generate mutant protein sequences based on wild-type sequences and mutation information
- Select the longest isoform from protein sequence data
- Load and manage protein sequence data from various databases
- Optimize memory usage through dynamic loading of protein databases

## Dependencies

- pandas
- numpy
- tqdm

Ensure these libraries are installed before using the Sequence Processor module.

## Components

### 1. Mutant Generator (`mutant_generator.py`)

The Mutant Generator provides functionality to:

- Parse mutation strings (e.g., 'A123T', 'Q58*')
- Select the longest isoform from a set of protein sequences
- Generate mutant sequences based on wild-type sequences and mutation information
- Add mutated sequences to a DataFrame containing gene and mutation information

Key functions:
- `parse_mutation(mutation)`
- `select_longest_isoform(isoform_sequences, db_name)`
- `mutate_sequence(wt_sequence, mutation_str)`
- `add_mutated_sequences(unique_mutations_df, protein_dict, db_name)`

### 2. Sequence Loader (`sequence_loader.py`)

The Sequence Loader is responsible for:

- Loading protein sequence data from CSV files
- Managing multiple protein databases
- Optimizing memory usage by dynamically loading databases as needed

Key components:
- `load_protein_dict(db_name)` function
- `ProteinDatabase` class

## Usage

Here's a basic example of how to use the Sequence Processor:

```python
from sequence_processor.mutant_generator import add_mutated_sequences
from sequence_processor.sequence_loader import ProteinDatabase
import pandas as pd

# Initialize the ProteinDatabase
db_names = ['uniprot', 'ncbi', 'pdb', 'ensembl']
protein_db = ProteinDatabase(db_names)

# Load mutation data
mutations_df = pd.read_csv('path/to/mutations.csv')

# Process mutations for each database
for db_name in db_names:
    protein_dict = protein_db.get_protein_dict(db_name)
    mutations_df = add_mutated_sequences(mutations_df, protein_dict, db_name)

# Save the results
mutations_df.to_csv('path/to/processed_mutations.csv', index=False)
```

## Notes

- Ensure that the protein sequence CSV files are located in the `data/processed/` directory with the naming convention `protein_sequences_{db_name}.csv`.
- The mutation string format should be in the form of 'A123T' or 'Q58*' for single mutations, or a comma-separated list for multiple mutations.
- The module handles various mutation types, including substitutions, truncations (*), and frameshifts (fs).
- Memory optimization is achieved by loading protein databases only when needed.