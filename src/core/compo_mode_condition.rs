// NCBI BLAST+ 2.17.0+, compo_mode_condition.c.
// Authors: NCBI BLAST developers. See LICENSE-NCBI.
use super::compo_frequency_data::blast_get_matrix_background_freq;
use super::{blast_get_relative_entropy, ECompoAdjustModes, EMatrixAdjustRule};
pub const HALF_CIRCLE_DEGREES: i32 = 180;
pub const PI: f64 = 3.1415926543;
pub const QUERY_MATCH_DISTANCE_THRESHOLD: f64 = 0.16;
pub const LENGTH_RATIO_THRESHOLD: f64 = 3.0;
pub const ANGLE_DEGREE_THRESHOLD: f64 = 70.0;
pub const HIGH_PAIR_THRESHOLD: f64 = 0.4;
pub const LENGTH_LOWER_THRESHOLD: i32 = 50;
pub fn s_high_pair_frequencies(letter_probs: &[f64], length: i32) -> bool {
    if length <= LENGTH_LOWER_THRESHOLD {
        return false;
    }
    let mut max = 0.0;
    let mut second = 0.0;
    for i in 0..20 {
        if letter_probs[i] > second {
            second = letter_probs[i];
            if letter_probs[i] > max {
                second = max;
                max = letter_probs[i];
            }
        }
    }
    max + second > HIGH_PAIR_THRESHOLD
}
pub fn s_high_pair_either_seq(
    p_query: &[f64],
    length1: i32,
    p_match: &[f64],
    length2: i32,
) -> bool {
    let result1 = s_high_pair_frequencies(p_query, length1);
    let result2 = s_high_pair_frequencies(p_match, length2);
    result1 || result2
}
pub fn s_never_adjust_matrix(
    _len_query: i32,
    _len_match: i32,
    _p_query: &[f64],
    _p_match: &[f64],
    _matrix_name: &[u8],
) -> EMatrixAdjustRule {
    EMatrixAdjustRule::DontAdjustMatrix
}
pub fn s_just_scale_old_matrix(
    _len_query: i32,
    _len_match: i32,
    _p_query: &[f64],
    _p_match: &[f64],
    _matrix_name: &[u8],
) -> EMatrixAdjustRule {
    EMatrixAdjustRule::CompoScaleOldMatrix
}
pub fn s_test_to_apply_re_adjustment_unconditional(
    _len_query: i32,
    _len_match: i32,
    _p_query: &[f64],
    _p_match: &[f64],
    _matrix_name: &[u8],
) -> EMatrixAdjustRule {
    EMatrixAdjustRule::UserSpecifiedRelEntropy
}
pub fn s_test_to_apply_re_adjustment_conditional(
    len_query: i32,
    len_match: i32,
    p_query: &[f64],
    p_match: &[f64],
    matrix_name: &[u8],
) -> EMatrixAdjustRule {
    let p_matrix = blast_get_matrix_background_freq(matrix_name)
        .expect("default BLOSUM62 background required");
    let mut query = [0.0; 20];
    let mut subject = [0.0; 20];
    let mut _corr_factor = 0.0;
    for i in 0..20 {
        query[i] = p_query[i];
        subject[i] = p_match[i];
        _corr_factor += (query[i] - p_matrix[i]) * (subject[i] - p_matrix[i]);
    }
    let d_m_mat = blast_get_relative_entropy(&subject, p_matrix);
    let d_q_mat = blast_get_relative_entropy(&query, p_matrix);
    let d_m_q = blast_get_relative_entropy(&subject, &query);
    let mut angle =
        ((d_m_mat * d_m_mat + d_q_mat * d_q_mat - d_m_q * d_m_q) / 2.0 / d_m_mat / d_q_mat).acos();
    angle = angle * HALF_CIRCLE_DEGREES as f64 / PI;
    let len_q = 1.0 * len_query as f64;
    let len_m = 1.0 * len_match as f64;
    let (len_large, len_small) = if len_q > len_m {
        (len_q, len_m)
    } else {
        (len_m, len_q)
    };
    if s_high_pair_either_seq(p_query, len_query, p_match, len_match) {
        EMatrixAdjustRule::UserSpecifiedRelEntropy
    } else if d_m_q > QUERY_MATCH_DISTANCE_THRESHOLD
        && len_large / len_small > LENGTH_RATIO_THRESHOLD
        && angle > ANGLE_DEGREE_THRESHOLD
    {
        EMatrixAdjustRule::CompoScaleOldMatrix
    } else {
        EMatrixAdjustRule::UserSpecifiedRelEntropy
    }
}
pub fn blast_choose_matrix_adjust_rule(
    length1: i32,
    length2: i32,
    prob_array1: &[f64],
    prob_array2: &[f64],
    matrix_name: &[u8],
    mode: ECompoAdjustModes,
) -> EMatrixAdjustRule {
    let condition: fn(i32, i32, &[f64], &[f64], &[u8]) -> EMatrixAdjustRule = match mode {
        ECompoAdjustModes::NoCompositionBasedStats => s_never_adjust_matrix,
        ECompoAdjustModes::CompositionBasedStats => s_just_scale_old_matrix,
        ECompoAdjustModes::CompositionMatrixAdjust => s_test_to_apply_re_adjustment_conditional,
        ECompoAdjustModes::CompoForceFullMatrixAdjust => {
            s_test_to_apply_re_adjustment_unconditional
        }
        ECompoAdjustModes::NumCompoAdjustModes => panic!("invalid composition mode"),
    };
    condition(length1, length2, prob_array1, prob_array2, matrix_name)
}
