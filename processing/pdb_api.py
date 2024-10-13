import requests
import json

class PDBAPI:
    """PDB API를 사용해 단백질 서열을 가져오는 클래스."""

    def search_pdb_by_gene(self, gene_name):
        """
        주어진 유전자 이름으로 PDB 엔트리를 검색합니다.
        반환값: (entry_id, entity_id)의 리스트
        """
        url = "https://search.rcsb.org/rcsbsearch/v2/query"
    
        # 검색 쿼리 작성
        query = {
            "query": {
                "type": "group",
                "logical_operator": "and",
                "nodes": [
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "rcsb_entity_source_organism.rcsb_gene_name.value",
                            "operator": "exact_match",
                            "value": gene_name
                        }
                    },
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "rcsb_entity_source_organism.ncbi_scientific_name",
                            "operator": "exact_match",
                            "value": "Homo sapiens"
                        }
                    },
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "entity_poly.rcsb_entity_polymer_type",
                            "operator": "exact_match",
                            "value": "Protein"
                        }
                    }
                ]
            },
            "return_type": "polymer_entity",
            "request_options": {
                "return_all_hits": True
            }
        }
    
        entries = []
        try:
            response = requests.post(url, json=query)
            response.raise_for_status()
        except requests.exceptions.RequestException as e:
            print(f"유전자 검색 중 오류 발생: {e}")
        else:
            try:
                results = response.json()
                if "result_set" in results and results["result_set"]:
                    # 각 엔트리의 entry_id와 entity_id 추출
                    for result in results["result_set"]:
                        identifier = result["identifier"]  # 예: '4HHB_1'
                        try:
                            entry_id, entity_id = identifier.split('_')
                            entries.append((entry_id, entity_id))
                        except ValueError:
                            print(f"Identifier {identifier}를 파싱하는 데 실패했습니다.")
                else:
                    print(f"유전자 '{gene_name}'에 대한 결과를 찾을 수 없습니다.")
            except (json.JSONDecodeError, ValueError) as e:
                print("검색 응답을 JSON으로 파싱하는 데 실패했습니다.")
                print(f"응답 내용: {response.text}")
    
        return entries

    def get_isoform_sequences(self, entry_entity_pairs):
        """
        엔트리 ID와 엔티티 ID 목록을 사용하여 isoform 서열을 가져옵니다.
        """
        isoform_sequences = {}
    
        if not entry_entity_pairs:
            return isoform_sequences
    
        url = "https://data.rcsb.org/graphql"
    
        # 엔트리 목록을 여러 개의 작은 배치로 나눕니다 (한 번에 너무 많은 요청을 보내지 않도록)
        batch_size = 100  # 필요한 경우 조정 가능
        for i in range(0, len(entry_entity_pairs), batch_size):
            batch_pairs = entry_entity_pairs[i:i+batch_size]
            identifiers = [f"{entry}_{entity}" for entry, entity in batch_pairs]
    
            # GraphQL 쿼리
            query = """
            query getSequences($entity_ids: [String!]!) {
              polymer_entities(entity_ids: $entity_ids) {
                rcsb_id
                entity_poly {
                  pdbx_seq_one_letter_code_can
                }
              }
            }
            """
    
            variables = {"entity_ids": identifiers}
    
            try:
                response = requests.post(url, json={"query": query, "variables": variables})
                response.raise_for_status()
            except requests.exceptions.RequestException as e:
                print(f"API 요청 중 오류 발생: {e}")
                continue
    
            try:
                data = response.json()
                if "errors" in data:
                    print(f"GraphQL 쿼리 실행 중 오류가 발생했습니다: {data['errors']}")
                    continue
    
                if "data" not in data or "polymer_entities" not in data["data"]:
                    print("예상한 데이터 구조가 API 응답에 없습니다.")
                    print(f"전체 응답: {json.dumps(data, indent=2)}")
                    continue
    
                for entity in data["data"]["polymer_entities"]:
                    try:
                        rcsb_id = entity["rcsb_id"]
                        sequence = entity["entity_poly"]["pdbx_seq_one_letter_code_can"]
                        if sequence:
                            isoform_sequences[rcsb_id] = sequence.replace("\n", "")  # 개행 문자 제거
                        else:
                            print(f"{rcsb_id}에 대한 서열을 찾을 수 없습니다.")
                    except (KeyError, TypeError) as e:
                        print(f"엔티티 데이터 처리 중 오류 발생: {e}")
                        print(f"문제의 엔티티: {json.dumps(entity, indent=2)}")
                        continue
            except (json.JSONDecodeError, ValueError):
                print("API 응답을 JSON으로 파싱하는 데 실패했습니다.")
                print(f"응답 내용: {response.text}")
                continue
    
        return isoform_sequences

    def get_pdb_isoforms(self, gene_name):
        """
        유전자 이름을 통해 해당 isoform 서열들을 가져옵니다.
        """
        entries = self.search_pdb_by_gene(gene_name)
        sequences = self.get_isoform_sequences(entries)
        return sequences
