# Mutation Encoder

This module provides functionality for classifying and processing genetic mutations commonly found in cancer research. It identifies mutation types such as missense, nonsense, frameshift, and more, while extracting mutation positions from mutation strings.

## Key Features

- Classifies mutations into various types (e.g., Missense, Nonsense, Frameshift).
- Extracts mutation positions from mutation strings.
- Handles duplicate mutations and provides warning messages.
- Supports the analysis of multiple mutations.

## Functions

- `classify_single_mutation(mutation)`: Classifies a single mutation and returns the mutation type and positions.
- `process_mutations(mutations)`: Processes a string of mutations and returns a collective mutation type and a sorted list of positions.

## Example Usage

To classify and process mutations, you can use the `process_mutations` function:

```python
from mutation_encoder import process_mutations

mutations = "L951V L851V L801V P34fs Q58*"
mutation_type, positions = process_mutations(mutations)

print("Mutation Type:", mutation_type)  # Output: Complex_mutation
print("Mutation Positions:", positions)  # Output: [34, 58, 801, 851, 951]
