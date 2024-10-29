# CancerSubtypesAI Mutation Encoder

The Mutation Encoder module of CancerSubtypesAI is a comprehensive tool designed for parsing, classifying, and processing genetic mutations commonly found in cancer research. It identifies various mutation types such as missense, nonsense, frameshift, and more, while extracting mutation information from mutation strings and calculating amino acid property changes.

## Features

- Parses mutation strings to extract original amino acid, position, and mutated amino acid.
- Classifies mutations into various types (e.g., Missense, Nonsense, Frameshift, Silent).
- Handles single mutations, multiple mutations, and complex mutation patterns.
- Supports the analysis of consecutive position mutations.
- Provides amino acid information lookup.
- Handles duplicate mutations by treating them as a single mutation.
- Calculates changes in amino acid properties for mutations.

## Dependencies

- numpy

Ensure these libraries are installed before using the Sequence Processor module.

## Mutation Patterns

The module recognizes and processes the following mutation patterns:

| Type                           | Pattern (Regex)                                  | Example            |
|--------------------------------|--------------------------------------------------|--------------------|
| Single Amino Acid Mutation (Missense) | `^([A-Z])(\d+)([A-Z])$`                   | `A123T`            |
| Nonsense Mutation              | `^([A-Z])(\d+)\*$`                               | `Q581*`            |
| Frameshift Mutation            | `^([A-Z]+)(\d+)fs(.*)$`                          | `P34fs`            |
| Frameshift Deletion Mutation   | `^([A-Z\-]?)(\d+)fs(.*)$`                        | `-62fs`            |
| Silent Nonsense Mutation       | `^\*(\d+)\*$`                                    | `*363*`            |
| Consecutive Position Mutation  | `^(\d+)_(\d+)([A-Z]{2})>([A-Z]{2})$`             | `49_50WE>*K`       |
| Single Deletion                | `^(\d+)del$`                                     | `490del`           |

## Mutation Classification

The module classifies mutations into the following types:

| Classification         | Description                                           | Example            |
|------------------------|-------------------------------------------------------|-------------------|
| Missense               | Single amino acid change                              | `A123T`           |
| Nonsense               | Mutation to stop codon                                | `Q581*`           |
| Silent_Missense        | No amino acid change                                  | `L145L`           |
| Silent_Nonsense        | No change in stop codon                               | `*363*`           |
| Frameshift             | Insertion or deletion causing frameshift              | `P34fs`           |
| Deletion               | Deletion of an amino acid                             | `490del`          |
| Complex_mutation       | Combination of different mutation types               | `A123T Q581* P34fs` |
| WT                     | Wild type (no mutation)                               | `WT`              |

## Main Functions

- `get_amino_acid_info(amino_acid)`: Returns the full name of an amino acid given its single-letter code.
- `parse_mutation(mutation)`: Parses a single mutation string and returns original amino acid, position, and mutated amino acid.
- `parse_multiple_mutations(mutation_str)`: Parses a string containing multiple mutations and returns lists of original amino acids, positions, and mutated amino acids.
- `classify_single_mutation(mutation)`: Classifies a single mutation and returns the mutation type and position(s).
- `classify_and_aggregate_mutations(mutations)`: Processes multiple mutations, classifies them, and determines an overall mutation type.
- `calculate_amino_acid_diff(origin, mutant, amino_acid_features)`: Calculates the difference in amino acid properties between the original and mutated amino acids.
- `process_mutation_features(row, amino_acid_features)`: Processes mutation features and calculates cumulative changes in amino acid properties.

## Handling Duplicate Mutations

When processing mutations, this module treats duplicate mutations as a single mutation. For example, if the input string contains `A123K A123K A123K`, it will be processed as if it were a single `A123K` mutation. This approach helps to handle potential input errors where the same mutation might be accidentally repeated.

## Mutation Type Determination Logic

The `classify_and_aggregate_mutations` function uses the following logic to determine the overall mutation type:

1. If all mutations are of the same type (including WT), that type is returned.
2. If all mutations are WT, 'WT' is returned.
3. If there are multiple non-WT mutation types, 'Complex_mutation' is returned.

## Amino Acid Property Changes

The module can calculate changes in amino acid properties for mutations, including:

- Hydrophobicity
- Polarity
- Molecular Weight (MW)
- Isoelectric Point (pI)
- Charge

These changes are calculated and accumulated for multiple mutations.

## Usage

To classify and process mutations:

```python
from mutation_encoder import classify_and_aggregate_mutations

mutations = "A123T Q581*"
mutation_type = classify_and_aggregate_mutations(mutations)

print("Overall Mutation Type:", mutation_type)  # Output: Complex_mutation
```

To parse multiple mutations:

```python
from mutation_encoder import parse_multiple_mutations

mutations = "A123T Q581*"
origins, positions, mutants = parse_multiple_mutations(mutations)

print("Original Amino Acids:", origins)  # Output: ['A', 'Q']
print("Positions:", positions)  # Output: [123, 581]
print("Mutated Amino Acids:", mutants)  # Output: ['T', '*']
```

To process mutation features:

```python
import pandas as pd
from mutation_encoder import process_mutation_features

# Assuming you have a DataFrame with mutation information and amino acid features
row = pd.Series({'type': 'Missense', 'origin': 'A', 'mutant': 'T'})
amino_acid_features = pd.DataFrame(...)  # Your amino acid features data

feature_changes = process_mutation_features(row, amino_acid_features)
print("Feature Changes:", feature_changes)
```

## Note

This module includes comprehensive error handling and provides warnings for duplicate mutations. It's designed to handle a wide range of mutation formats commonly encountered in genetic research and provides detailed analysis of amino acid property changes resulting from mutations.