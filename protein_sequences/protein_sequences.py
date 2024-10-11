import requests
from Bio import Entrez, SeqIO
from tqdm import tqdm

def get_uniprot_isoforms(gene_name):
    """
    Fetch all isoform protein sequences for a given gene from UniProt.

    Parameters:
    - gene_name (str): Name of the gene to search for.

    Returns:
    - dict: A dictionary where keys are accession numbers and values are protein sequences.
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
        # Extract accession from header
        fields = header.split('|')
        if len(fields) >= 3:
            accession = fields[1]
            protein_seqs[accession] = sequence
        else:
            print(f"Failed to parse accession from header: {header}")
    
    return protein_seqs

def get_entrez_isoforms(gene_names, protein_sequences, email):
    """
    Fetch protein sequences for a list of genes from NCBI Entrez and update the given dictionary.

    Parameters:
    - gene_names (list): List of gene names to fetch protein sequences for.
    - protein_sequences (dict): Dictionary to update with fetched protein sequences.
    - email (str): Email address for NCBI Entrez API.

    Returns:
    - None
    """
    Entrez.email = email  # Set email for NCBI Entrez API
    for gene_name in tqdm(gene_names, desc="Processing genes"):
        print(f"\nFetching protein sequences for gene: {gene_name}")
        search_term = f"{gene_name}[Gene Name] AND Homo sapiens[Organism]"
        try:
            # Search for protein IDs
            with Entrez.esearch(db="protein", term=search_term, retmax=1000) as handle:
                record = Entrez.read(handle)
            id_list = record.get("IdList", [])

            if not id_list:
                print(f"No protein IDs found for gene: {gene_name}")
                continue

            # Batch fetch protein sequences
            id_str = ','.join(id_list)
            with Entrez.efetch(db="protein", id=id_str, rettype="fasta", retmode="text") as fetch_handle:
                seq_records = list(SeqIO.parse(fetch_handle, "fasta"))

            if not seq_records:
                print(f"No protein sequences found for gene: {gene_name}")
                continue

            protein_seqs = {}
            for seq_record in tqdm(seq_records, desc=f"Processing sequences for {gene_name}", leave=False):
                accession = extract_accession(seq_record.id)
                if not accession:
                    print(f"Could not extract accession: {seq_record.id}")
                    continue
                if accession in protein_seqs:
                    continue  # Avoid duplicates
                sequence = str(seq_record.seq)
                protein_seqs[accession] = sequence

            if protein_seqs:
                protein_sequences[gene_name] = protein_seqs
                print(f"Successfully fetched protein sequences for gene: {gene_name}")
            else:
                print(f"No valid protein sequences found for gene: {gene_name}")

        except Exception as e:
            print(f"Error fetching data for gene {gene_name}: {e}")

def extract_accession(seq_id):
    """
    Extract the accession number from a sequence ID.

    Parameters:
    - seq_id (str): Sequence ID string.

    Returns:
    - str or None: Accession number or None if not found.
    """
    if '|' in seq_id:
        parts = seq_id.split('|')
        # Find part that matches common accession prefixes
        for part in parts:
            if part.startswith(('NP_', 'XP_', 'WP_', 'YP_', 'AP_')):
                return part
        # Use the last part if no common prefix is found
        return parts[-1] if parts[-1] else None
    else:
        return seq_id if seq_id else None
