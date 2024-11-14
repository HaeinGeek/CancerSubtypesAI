import os
import subprocess
import time
import pandas as pd
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIWWW, NCBIXML
import multiprocessing as mp
from tqdm.auto import tqdm
import logging

# 로깅 설정 - 로그 파일을 'my_log.log'에 저장
logging.basicConfig(
    level=logging.INFO, 
    format='%(asctime)s - %(levelname)s - %(message)s',
    filename='logs/fetch_msa.log',  # 로그 파일 이름 지정
    filemode='w'  # 파일 모드 'a'는 추가 모드 ('w'는 덮어쓰기 모드)
)

class GeneWiseBLASTMSA:
    def __init__(self, input_df, output_dir, email, api_key=None, n_jobs=1):
        self.input_df = input_df
        self.output_dir = output_dir
        self.n_jobs = n_jobs if n_jobs > 0 else 1  # 코랩 환경에서는 n_jobs=1 권장

        # NCBI Entrez 설정
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key

        os.makedirs(output_dir, exist_ok=True)

    def clean_sequence(self, seq):
        valid_chars = set("ACDEFGHIKLMNPQRSTVWY")
        return ''.join([char for char in seq.upper() if char in valid_chars])

    def run_blast(self, gene_name, sequences, isoform_ids, max_hits, e_value_thresh, txid, identity_thresh, coverage_thresh):
        logging.info(f"{gene_name} 유전자의 BLAST 수행 중...")
        gene_output_dir = os.path.join(self.output_dir, gene_name)
        os.makedirs(gene_output_dir, exist_ok=True)
        
        # 원본 서열 저장
        query_records = [SeqRecord(seq, id=isoform_id, description="") for seq, isoform_id in zip(sequences, isoform_ids)]
        query_file = os.path.join(gene_output_dir, f"{gene_name}_query.fasta")
        SeqIO.write(query_records, query_file, "fasta")
        
        # 모든 서열의 BLAST 결과 서열 저장을 위한 리스트
        all_hit_records = []
        
        for record in query_records:
            try:
                # BLASTP 검색 수행
                result_handle = NCBIWWW.qblast(
                    program="blastp",
                    database="nr",
                    sequence=record.seq,
                    hitlist_size=max_hits,
                    expect=e_value_thresh,
                    format_type="XML",
                    entrez_query=txid
                )
                blast_records = NCBIXML.read(result_handle)
                result_handle.close()  # 결과 핸들 닫기
                
                # BLAST 결과가 있는지 확인
                if len(blast_records.alignments) == 0:
                    logging.warning(f"{record.id} 서열의 BLAST 결과가 없습니다.")
                    continue
                
                # 히트된 서열 수집 및 필터링
                for alignment in blast_records.alignments:
                    for hsp in alignment.hsps:
                        identity = (hsp.identities / hsp.align_length) * 100
                        coverage = (hsp.align_length / len(record.seq)) * 100
                        if identity >= identity_thresh and coverage >= coverage_thresh:
                            hit_seq = hsp.sbjct
                            hit_id = alignment.hit_id
                            hit_def = alignment.hit_def
                            hit_record = SeqRecord(Seq(hit_seq), id=hit_id, description=hit_def)
                            all_hit_records.append(hit_record)
                # BLAST 요청 후 지연 시간 추가
                time.sleep(1)  # 1초 지연
            except Exception as e:
                logging.error(f"{record.id} 서열의 BLAST 수행 중 오류 발생: {e}")
                if 'result_handle' in locals():
                    result_handle.close()  # 예외 발생 시에도 핸들 닫기
                continue        
                
        # 중복 제거
        unique_hit_records = {rec.id: rec for rec in all_hit_records}.values()
        
        # 히트된 서열 저장
        hits_file = os.path.join(gene_output_dir, f"{gene_name}_hits.fasta")
        SeqIO.write(unique_hit_records, hits_file, "fasta")
        
        return query_file, hits_file

    def run_mafft(self, gene_name, query_file, hits_file):
        logging.info(f"{gene_name} 유전자의 MSA 수행 중...")
        gene_output_dir = os.path.join(self.output_dir, gene_name)

        # 서열 합치기
        combined_file = os.path.join(gene_output_dir, f"{gene_name}_combined.fasta")
        with open(combined_file, "w") as outfile:
            for fname in [query_file, hits_file]:
                with open(fname) as infile:
                    outfile.write(infile.read())

        # MAFFT 실행 (subprocess를 사용하여 반환 코드 확인)
        msa_output_file = os.path.join(gene_output_dir, f"{gene_name}_msa.fasta")
        mafft_command = ["mafft", "--auto", "--thread", "1", combined_file]

        try:
            with open(msa_output_file, "w") as outfile:
                process = subprocess.run(mafft_command, stdout=outfile, stderr=subprocess.PIPE, text=True)
            if process.returncode != 0:
                logging.error(f"MAFFT 오류 ({gene_name}): {process.stderr}")
                raise RuntimeError(f"MAFFT 실행 중 오류가 발생했습니다 ({gene_name}).")
        except Exception as e:
            logging.error(f"{gene_name} 유전자의 MSA 수행 중 예외 발생: {e}")
            raise e

        logging.info(f"{gene_name} 유전자의 MSA 완료. 결과가 {msa_output_file}에 저장되었습니다.")
        return msa_output_file

    def process_gene(self, gene_name, group, max_hits, e_value_thresh, txid, identity_thresh, coverage_thresh):
        sequences = [Seq(self.clean_sequence(seq)) for seq in group["wt_seq"]]
        isoform_ids = group["isoform_id"].astype(str).tolist()
        
        # BLAST 수행
        query_file, hits_file = self.run_blast(
            gene_name, sequences, isoform_ids, max_hits, e_value_thresh, txid, identity_thresh, coverage_thresh
        )
        
        # 히트된 서열이 있는 경우에만 MSA 수행
        if os.path.exists(hits_file) and os.path.getsize(hits_file) > 0:
            msa_file = self.run_mafft(gene_name, query_file, hits_file)
            return msa_file
        else:
            logging.warning(f"{gene_name} 유전자의 BLAST 결과가 없거나 히트된 서열이 없습니다.")
            return None

    def run_pipeline(self, max_hits, e_value_thresh, txid, identity_thresh, coverage_thresh):
        try:
            # 유전자별로 데이터 그룹화
            grouped = self.input_df.groupby('gene')
            gene_list = grouped.groups.keys()
            
            # 멀티프로세싱을 위한 입력 데이터 생성
            pool_input = []
            for gene_name, group in grouped:
                pool_input.append((gene_name, group, max_hits, e_value_thresh, txid, identity_thresh, coverage_thresh))
            
            # 병렬 처리
            results = []
            with mp.Pool(processes=self.n_jobs) as pool:
                for res in tqdm(pool.imap_unordered(self.process_gene_wrapper, pool_input), total=len(pool_input)):
                    results.append(res)
                    
            logging.info("모든 유전자의 BLAST 및 MSA가 완료되었습니다.")
            return results
        except Exception as e:
            logging.error(f"오류 발생: {e}")
            return None

    def process_gene_wrapper(self, args):
        return self.process_gene(*args)