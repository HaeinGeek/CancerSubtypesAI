from sequence_processor.sequence_loader import load_protein_dict, ProteinDatabase
from sequence_processor.mutant_generator import(
    select_longest_isoform,
    parse_mutation,
    mutate_sequence,
    add_mutated_sequences    
)

__all__ = [
    'load_protein_dict',
    'ProteinDatabase',
    'select_longest_isoform',
    'parse_mutation',
    'mutate_sequence',
    'add_mutated_sequences'
]