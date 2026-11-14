use btyper_blast::nucleotide_setup::NucleotideInitialSearch;
use btyper_blast::{
    DatabaseBuilder, NucleotideTask, PreparedNucleotideBatch, PreparedNucleotideQuery,
};

#[test]
fn standard_tasks_assemble_and_extend_exact_subject_match() {
    let mut state = 13791u32;
    let query: Vec<u8> = (0..128)
        .map(|_| {
            state = state.wrapping_mul(1664525).wrapping_add(1013904223);
            b"ACGT"[(state >> 28) as usize % 4]
        })
        .collect();
    let mut subject = b"GGCATCA".to_vec();
    subject.extend_from_slice(&query);
    subject.extend_from_slice(b"CGTACCA");
    let mut builder = DatabaseBuilder::new();
    builder.add_sequence("contig", &subject).unwrap();
    let db = builder.build().unwrap();
    for task in [NucleotideTask::Megablast, NucleotideTask::Blastn] {
        let batch =
            PreparedNucleotideBatch::new(&[
                PreparedNucleotideQuery::new("query", &query, task).unwrap()
            ])
            .unwrap();
        let mut search = NucleotideInitialSearch::new(&batch, &db).unwrap();
        let reward = search.scoring_options.reward as i32;
        let hits = search.scan_subject(db.sequence(0).unwrap()).unwrap();
        assert!(hits
            .iter()
            .any(|h| h.ungapped_data.is_some_and(|d| d.q_start == 0
                && d.s_start == 7
                && d.length == 128
                && d.score == 128 * reward)));
        assert!(search.stats().lookup_hits > 0);
        let preliminary = search.preliminary_subject(db.sequence(0).unwrap()).unwrap();
        assert!(preliminary.hsp_array[..preliminary.hspcnt as usize]
            .iter()
            .flatten()
            .any(|h| h.context == 0
                && h.query.offset == 0
                && h.query.end == 128
                && h.subject.offset == 7
                && h.subject.end == 135
                && h.score == 128 * reward));
        let final_hits = search.traceback_subject(db.sequence(0).unwrap()).unwrap();
        let hsp = final_hits.hsp_array[..final_hits.hspcnt as usize]
            .iter()
            .flatten()
            .find(|h| h.context == 0)
            .unwrap();
        assert_eq!((hsp.query.offset, hsp.query.end), (0, 128));
        assert_eq!((hsp.subject.offset, hsp.subject.end), (7, 135));
        assert_eq!(hsp.score, 128 * reward);
        assert_eq!(hsp.num_ident, 128);
        let script = hsp.gap_info.as_ref().unwrap();
        assert_eq!(script.size, 1);
        assert_eq!(script.num[0], 128);
        assert!(hsp.bit_score > 0.0 && hsp.evalue < 1e-20);
        let report = btyper_blast::report::AlignmentReport::nucleotide(
            hsp,
            &batch.blastna_contexts()[..128],
            &db.sequence(0).unwrap().blastna(),
        )
        .unwrap();
        assert_eq!(
            (report.qstart, report.qend, report.sstart, report.send),
            (1, 128, 8, 135)
        );
        assert_eq!(report.btop, "128");
        assert_eq!(report.qseq.as_bytes(), query);
        assert_eq!(report.qseq, report.sseq);
    }
}
