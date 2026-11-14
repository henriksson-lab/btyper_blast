use btyper_blast::{PreparedProteinBatch, PreparedProteinQuery};

#[test]
fn protein_batch_retains_originals_and_hard_masks_seg_ranges() {
    let varied = b"ACDEFGHIKLMNPQRSTVWY".repeat(4);
    let queries = [
        PreparedProteinQuery::new("varied", &varied).unwrap(),
        PreparedProteinQuery::new("low", &[b'A'; 48]).unwrap(),
    ];
    assert!(queries[0].seg_masks().is_empty());
    assert_eq!(queries[1].seg_masks().len(), 1);
    let batch = PreparedProteinBatch::new(&queries).unwrap();
    let info = batch.query_info();
    assert_eq!(info.contexts.len(), 2);
    assert_eq!(info.contexts[0].query_offset, 0);
    assert_eq!(info.contexts[1].query_offset, 81);
    for (i, query) in queries.iter().enumerate() {
        let start = info.contexts[i].query_offset as usize;
        assert_eq!(
            &batch.original_contexts()[start..start + query.len()],
            query.ncbistdaa()
        );
        let expected = if i == 0 {
            query.ncbistdaa().to_vec()
        } else {
            vec![21; query.len()]
        };
        assert_eq!(
            &batch.encoded_contexts()[start..start + query.len()],
            expected
        );
        assert_eq!(batch.sequence_with_sentinels()[start], 0);
    }
    assert_eq!(
        batch.original_contexts().len(),
        batch.encoded_contexts().len()
    );
}
