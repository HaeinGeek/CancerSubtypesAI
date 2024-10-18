from embedding.protbert_embedding import (
    extract_data_for_embedding,
    load_protbert,
    process_and_save_embeddings
)

from utils.utils import(
    setup_environment,
    load_logger,
    load_full_feature_data
)

def main():
    setup_environment()
    logger = load_logger()

    # 임베딩 계산할 데이터 로딩
    mutation_df = load_full_feature_data()
    wt_prot, mut_prot = extract_data_for_embedding(mutation_df)

    print(f'전체 샘플 길이:{len(mutation_df)}',
        f'\nwt 서열 길이: {len(wt_prot)}',
        f'\n변이 서열 길이: {len(mut_prot)}'
    )

    # ProtBERT tokenizer & model 로딩
    tokenizer, model, device = load_protbert()
    

    logger.info("Starting processing of wild type sequences")
    process_and_save_embeddings(wt_prot, 
                                tokenizer, 
                                model, 
                                device,
                                'data/processed/wt_embedding_full.pt')
    
    logger.info("Starting processing of mutant sequences")
    process_and_save_embeddings(mut_prot, 
                                tokenizer, 
                                model, 
                                device, 
                                'data/processed/mut_embedding_full.pt',
                                is_mutant=True)

    logger.info("Embedding extraction process completed")

if __name__ == "__main__":
    main()