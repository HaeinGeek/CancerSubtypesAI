from processing.uniprot_api import UniProtAPI
from processing.pdb_api import PDBAPI
from processing.ncbi_api import NCBIAPI
from processing.ensembl_api import EnsemblAPI

class SequenceFetcher:
    """여러 API를 통합해 단백질 서열을 가져오는 클래스."""

    def __init__(self):
        self.uniprot_api = UniProtAPI()
        self.pdb_api = PDBAPI()
        self.ncbi_api = NCBIAPI()
        self.ensembl_api = EnsemblAPI()

    # 메서드들을 각 API 클래스의 메서드로 위임합니다.
    def get_uniprot_sequences(self, *args, **kwargs):
        return self.uniprot_api.get_uniprot_sequences(*args, **kwargs)

    def get_pdb_isoforms(self, *args, **kwargs):
        return self.pdb_api.get_pdb_isoforms(*args, **kwargs)

    def get_ncbi_sequences(self, *args, **kwargs):
        return self.ncbi_api.get_ncbi_sequences(*args, **kwargs)

    def get_ensembl_isoforms(self, *args, **kwargs):
        return self.ensembl_api.get_ensembl_isoforms(*args, **kwargs)
