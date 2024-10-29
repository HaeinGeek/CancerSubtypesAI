from embedding.protbert_embedding import (
    extract_data_for_embedding,
    load_protbert,
    extract_embeddings,
    process_embeddings,
    save_embeddings_to_zip,
    load_compressed_embeddings
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

    # wt 임베딩 계산
    logger.info("Starting processing of wild type sequences")
    wt_embeddings = process_embeddings(wt_prot, tokenizer, model, device, batch_size=32)

    # mut 임베딩 계산
    logger.info("Starting processing of mutant sequences")
    mut_embeddings = process_embeddings(mut_prot, tokenizer, model, device, is_mutant=True, batch_size=32)

    # 임베딩 압축파일 저장
    logger.info("Saving wild type embeddings")
    save_embeddings_to_zip(wt_embeddings, 'data/processed/embeddings/wt_embedding_full.pt')

    logger.info("Saving mutant embeddings")
    save_embeddings_to_zip(mut_embeddings, 'data/processed/embeddings/mut_embedding_full.pt')

    logger.info("Embedding extraction process completed")

if __name__ == "__main__":
    main()