from processing.uniprot_api import UniProtAPI
from processing.pdb_api import PBDAPI
from processing.ncbi_api import NCBIAPI
from processing.ensembl_api import EnsemblAPI

class SequenceFetcher(EnsemblAPI, UniProtAPI, PBDAPI, NCBIAPI):
    pass
