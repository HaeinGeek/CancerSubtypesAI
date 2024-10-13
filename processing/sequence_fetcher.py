from processing.uniprot_api import UniProtAPI
from processing.pdb_api import PBDAPI
from processing.ncbi_api import NCBIAPI
from processing.ensembl_api import EnsemblAPI

class SequenceFetcher(UniProtAPI, PDBAPI, NCBIAPI, EnsemblAPI):
    """여러 API를 통합해 단백질 서열을 가져오는 클래스."""
    pass
