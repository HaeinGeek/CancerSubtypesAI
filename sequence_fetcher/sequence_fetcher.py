import pandas as pd
from tqdm.auto import tqdm

from sequence_fetcher.uniprot_api import UniProtAPI
from sequence_fetcher.ncbi_api import NCBIAPI
from sequence_fetcher.pdb_api import PDBAPI
from sequence_fetcher.ensembl_api import EnsemblAPI

class SequenceFetcher:
    """여러 API를 통합해 단백질 서열을 가져오는 클래스."""

    def __init__(self):
        self.uniprot_api = UniProtAPI()
        self.ncbi_api = NCBIAPI()
        self.pdb_api = PDBAPI()
        self.ensembl_api = EnsemblAPI()

    def get_isoform_function(self, db_name):
        """
        주어진 데이터베이스 이름에 해당하는 isoform 가져오기 함수를 반환합니다.
        """
        isoform_functions = {
            'uniprot': self.uniprot_api.get_uniprot_isoforms,
            'ncbi': self.ncbi_api.get_ncbi_isoforms,
            'pdb': self.pdb_api.get_pdb_isoforms,
            'ensembl': self.ensembl_api.get_ensembl_isoforms
        }

        if db_name not in isoform_functions:
            raise ValueError(f"지원되지 않는 데이터베이스입니다: {db_name}")

        return isoform_functions[db_name]

    def create_protein_dict(self, genes, get_isoforms_func):
        """
        각 유전자에 대해 isoform 서열을 가져와 딕셔너리에 저장합니다.

        Parameters:
        - genes (list): 처리할 유전자 이름들의 리스트.
        - get_isoforms_func (function): isoform 서열을 가져오는 함수.

        Returns:
        - tuple: (protein_dict, not_found_genes)
            - protein_dict: 유전자와 해당 isoform 서열을 저장한 딕셔너리.
            - not_found_genes: isoform을 찾지 못한 유전자 리스트.
        """
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

        return protein_dict

    def process_isoforms(self, db_name, genes, **kwargs):
        """
        주어진 데이터베이스에 대해 각 유전자에 대한 isoform 서열을 가져와 처리하고 저장합니다.

        Parameters:
        - db_name (str): 데이터베이스 이름 ('uniprot', 'entrez', 'pdb', 'ensembl').
        - genes (list): 처리할 유전자 이름들의 리스트.
        - kwargs: isoform 함수에 필요한 추가 매개변수.

        Returns:
        - pd.DataFrame: 유전자와 isoform 서열 정보가 포함된 데이터프레임.
        """
        get_isoforms = self.get_isoform_function(db_name)  # 적절한 isoform 함수 선택
        protein_dict = self.create_protein_dict(genes, lambda gene: get_isoforms(gene, **kwargs))

        # protein_dict를 데이터프레임으로 변환
        protein_df = pd.DataFrame(list(protein_dict.items()), columns=['gene', 'sequence'])

        # CSV 파일로 저장
        filename = f'data/processed/protein_sequences_{db_name}.csv'
        protein_df.to_csv(filename, index=False)
        print(f"Saved {db_name} isoform sequences to {filename}")

        return protein_df