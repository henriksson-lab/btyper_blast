//! Two focused original-C parameter setups, recorded from BLAST+ 2.17.0+.
use btyper_blast::core::{
    blast_options::*, blast_parameters::*, blast_program::EBlastProgramType, blast_query_info::*,
    blast_stat::*, blast_stat_gapped::*,
};

#[test]
fn nucleotide_contexts_feed_original_default_word_cutoffs() {
    // C outputs: score_min, xdrop_max, cutoff, xdrop_init, reduced, round_down.
    let expected = [(13, 11, 13, 11, 10, false), (23, 22, 23, 22, 18, true)];
    for (profile, (minimum, xmax, cutoff, xinit, reduced, rounded)) in
        expected.into_iter().enumerate()
    {
        let (reward, penalty, gap_open, gap_extend) = if profile == 0 {
            (1, -2, 0, 0)
        } else {
            (2, -3, 5, 2)
        };
        let mut matrix = [[0; 16]; 16];
        assert_eq!(
            blast_score_blk_nucl_matrix_create(reward, penalty, &mut matrix),
            0
        );
        let mut score = BlastScoreBlk {
            alphabet_code: 99,
            alphabet_size: 16,
            ambiguous_res: vec![14],
            matrix: matrix.iter().map(|r| r.to_vec()).collect(),
            reward,
            penalty,
            scale_factor: 1.0,
            ..Default::default()
        };
        blast_score_blk_max_score_set(&mut score);
        let mut info = blast_query_info_new(EBlastProgramType::Blastn, 1).unwrap();
        info.contexts[0].query_length = 32;
        info.contexts[1].query_offset = 33;
        info.contexts[1].query_length = 32;
        let mut query: Vec<u8> = (0..32).map(|i| i % 4).collect();
        query.push(15);
        query.extend_from_slice(&[14; 32]);
        let mut warnings = vec![];
        assert_eq!(
            blast_score_blk_kbp_ungapped_calc(
                EBlastProgramType::Blastn,
                &mut score,
                &query,
                &mut info,
                Some(&mut warnings)
            ),
            0
        );
        assert!(info.contexts[0].is_valid);
        assert!(!info.contexts[1].is_valid);
        assert_eq!(warnings.len(), 1);
        assert_eq!(warnings[0].context, 1);
        let mut gapped = blast_karlin_blk_new();
        let mut round_down = false;
        assert_eq!(
            blast_karlin_blk_nucl_gapped_calc(
                &mut gapped,
                gap_open,
                gap_extend,
                reward,
                penalty,
                score.kbp_std[0].as_ref().unwrap(),
                &mut round_down,
                None
            ),
            0
        );
        assert_eq!(round_down, rounded);
        score.kbp_gap = Some(vec![Some(gapped), None]);
        let options = blast_initial_word_options_new(EBlastProgramType::Blastn);
        let hits = BlastHitSavingParameters {
            cutoffs: vec![
                BlastGappedCutoffs {
                    cutoff_score: 0,
                    cutoff_score_max: 23,
                },
                BlastGappedCutoffs {
                    cutoff_score: 0,
                    cutoff_score_max: i32::MAX,
                },
            ],
            ..Default::default()
        };
        let mut parameters = None;
        assert_eq!(
            blast_initial_word_parameters_new(
                EBlastProgramType::Blastn,
                &options,
                &hits,
                &score,
                &info,
                1000,
                Some(&mut parameters)
            ),
            0
        );
        let parameters = parameters.unwrap();
        assert_eq!(
            (
                parameters.cutoff_score_min,
                parameters.x_dropoff_max,
                parameters.cutoffs[0].cutoff_score,
                parameters.cutoffs[0].x_dropoff_init,
                parameters.cutoffs[0].reduced_nucl_cutoff_score
            ),
            (minimum, xmax, cutoff, xinit, reduced)
        );
        assert_eq!(parameters.cutoffs[1].cutoff_score, i32::MAX);
        assert_eq!(parameters.container_type, ESeedContainerType::DiagArray);
        assert_eq!(parameters.nucl_score_table[0], 4 * reward);
        assert_eq!(parameters.nucl_score_table[255], 4 * penalty);
    }
}
