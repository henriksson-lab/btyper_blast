// NCBI BLAST+ 2.17.0+, default composition_adjustment.c dependency closure.
// Authors: NCBI BLAST developers. See LICENSE-NCBI.
use super::blast_stat::{blast_karlin_lambda_nr, BlastScoreFreq};
use super::compo_linear_algebra::nlm_add_vectors;
use super::compo_optimize_target_freq::blast_optimize_target_frequencies;
#[path = "compo_blosum62_data.rs"]
mod compo_blosum62_data;
#[path = "compo_frequency_data.rs"]
mod compo_frequency_data;
#[path = "compo_mode_condition.rs"]
mod compo_mode_condition;
use compo_blosum62_data::BLOSUM62_FREQRATIOS;
pub use compo_frequency_data::{
    blast_frequency_data_is_available, blast_get_joint_probs_for_matrix,
    blast_get_matrix_background_freq,
};
pub use compo_mode_condition::blast_choose_matrix_adjust_rule;

pub const COMPO_NUM_TRUE_AA: usize = 20;
pub const COMPO_LARGEST_ALPHABET: usize = 28;
pub const COMPO_SCORE_MIN: i32 = i16::MIN as i32;
pub const TRUE_CHAR_POSITIONS: [usize; 20] = [
    1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22,
];
pub const ALPHA_CONVERT: [i32; 28] = [
    -1, 0, -1, 4, 3, 6, 13, 7, 8, 9, 11, 10, 12, 2, 14, 5, 1, 15, 16, 19, 17, -1, 18, -1, -1, -1,
    -1, -1,
];
pub const K_COMPOSITION_MARGIN: i32 = 20;
pub const K_COMPO_ADJUST_ERR_TOLERANCE: f64 = 0.00000001;
pub const K_COMPO_ADJUST_ITERATION_LIMIT: i32 = 2000;
pub const K_FIXED_RE_BLOSUM62: f64 = 0.44;
pub const K_MAXIMUM_XSCORE: f64 = -1.0;
pub const LAMBDA_RATIO_LOWER_BOUND: f64 = 0.5;
pub const SCALING_FACTOR: i32 = 32;
pub const K_RE_MATRIX_ADJUSTMENT_PSEUDOCOUNTS: i32 = 20;
pub const E_GAP_CHAR: usize = 0;
pub const E_BCHAR: usize = 2;
pub const E_CCHAR: usize = 3;
pub const E_DCHAR: usize = 4;
pub const E_ECHAR: usize = 5;
pub const E_ICHAR: usize = 9;
pub const E_LCHAR: usize = 11;
pub const E_NCHAR: usize = 13;
pub const E_QCHAR: usize = 15;
pub const E_XCHAR: usize = 21;
pub const E_ZCHAR: usize = 23;
pub const E_SELENOCYSTEINE: usize = 24;
pub const E_STOP_CHAR: usize = 25;
pub const E_OCHAR: usize = 26;
pub const E_JCHAR: usize = 27;
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[repr(i32)]
pub enum ECompoAdjustModes {
    NoCompositionBasedStats = 0,
    CompositionBasedStats = 1,
    CompositionMatrixAdjust = 2,
    CompoForceFullMatrixAdjust = 3,
    NumCompoAdjustModes = 4,
}
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[repr(i32)]
pub enum EMatrixAdjustRule {
    DontAdjustMatrix = -1,
    CompoScaleOldMatrix = 0,
    UnconstrainedRelEntropy = 1,
    RelEntropyOldMatrixNewContext = 2,
    RelEntropyOldMatrixOldContext = 3,
    UserSpecifiedRelEntropy = 4,
}
#[derive(Clone, Debug)]
pub struct BlastAminoAcidComposition {
    pub prob: [f64; 28],
    pub num_true_amino_acids: i32,
}
impl Default for BlastAminoAcidComposition {
    fn default() -> Self {
        Self {
            prob: [0.0; 28],
            num_true_amino_acids: 0,
        }
    }
}
#[derive(Debug)]
pub struct BlastMatrixInfo {
    pub matrix_name: Vec<u8>,
    pub start_matrix: Vec<i32>,
    pub start_freq_ratios: Vec<f64>,
    pub rows: usize,
    pub cols: usize,
    pub position_based: bool,
    pub ungapped_lambda: f64,
}
#[derive(Debug)]
pub struct BlastCompositionWorkspace {
    pub mat_b: Vec<f64>,
    pub mat_final: Vec<f64>,
    pub first_standard_freq: Vec<f64>,
    pub second_standard_freq: Vec<f64>,
}
pub type CalcLambda = fn(&[f64], i32, i32, f64) -> f64;

pub fn blast_apply_pseudocounts(
    probs20: &mut [f64],
    number_of_observations: i32,
    background_probs20: &[f64],
    pseudocounts: i32,
) {
    let dpseudocounts = pseudocounts as f64;
    let mut sum = 0.0;
    for i in 0..20 {
        sum += probs20[i];
    }
    if sum == 0.0 {
        sum = 1.0;
    }
    let weight = dpseudocounts / (number_of_observations as f64 + dpseudocounts);
    for i in 0..20 {
        probs20[i] = (1.0 - weight) * probs20[i] / sum + weight * background_probs20[i];
    }
}
pub fn blast_get_relative_entropy(a: &[f64], b: &[f64]) -> f64 {
    let mut value = 0.0;
    for i in 0..20 {
        let temp = (a[i] + b[i]) / 2.0;
        if temp > 0.0 {
            if a[i] > 0.0 {
                value += a[i] * (a[i] / temp).ln() / 2.0;
            }
            if b[i] > 0.0 {
                value += b[i] * (b[i] / temp).ln() / 2.0;
            }
        }
    }
    if value < 0.0 {
        value = 0.0;
    }
    value.sqrt()
}
pub fn blast_calc_freq_ratios(
    ratios: &mut [f64],
    alphsize: usize,
    row_prob: &[f64],
    col_prob: &[f64],
) {
    for i in 0..alphsize {
        if row_prob[i] > 0.0 {
            for j in 0..alphsize {
                if col_prob[j] > 0.0 {
                    ratios[i * alphsize + j] /= row_prob[i] * col_prob[j];
                }
            }
        }
    }
}
pub fn blast_freq_ratio_to_score(matrix: &mut [f64], rows: usize, cols: usize, lambda: f64) {
    for i in 0..rows {
        for j in 0..cols {
            let p = i * cols + j;
            if matrix[p] == 0.0 {
                matrix[p] = COMPO_SCORE_MIN as f64;
            } else {
                matrix[p] = matrix[p].ln() / lambda;
            }
        }
    }
}
pub fn s_gather_letter_probs(output: &mut [f64], input: &[f64], alphsize: usize) {
    for c in 0..alphsize {
        if ALPHA_CONVERT[c] != -1 {
            output[ALPHA_CONVERT[c] as usize] = input[c];
        }
    }
}
pub fn s_unpack_letter_probs(std_probs: &mut [f64], alphsize: usize, probs: &[f64]) {
    for c in 0..alphsize {
        std_probs[c] = if ALPHA_CONVERT[c] != -1 {
            probs[ALPHA_CONVERT[c] as usize]
        } else {
            0.0
        };
    }
}
pub fn s_set_pair_ambig_probs_to_sum(probs: &mut [f64], alphsize: usize) {
    probs[E_BCHAR] = probs[E_DCHAR] + probs[E_NCHAR];
    probs[E_ZCHAR] = probs[E_ECHAR] + probs[E_QCHAR];
    if alphsize > E_JCHAR {
        probs[E_JCHAR] = probs[E_ICHAR] + probs[E_LCHAR];
    }
}

pub fn blast_true_aa_to_std_target_freqs(std_freq: &mut [f64], std_alphsize: usize, freq: &[f64]) {
    let mut sum = 0.0;
    for a in 0..20 {
        for b in 0..20 {
            sum += freq[a * 20 + b];
        }
    }
    for a_big in 0..std_alphsize {
        if ALPHA_CONVERT[a_big] < 0 {
            for b_big in 0..std_alphsize {
                std_freq[a_big * std_alphsize + b_big] = 0.0;
            }
        } else {
            let a = ALPHA_CONVERT[a_big] as usize;
            for b_big in 0..std_alphsize {
                std_freq[a_big * std_alphsize + b_big] = if ALPHA_CONVERT[b_big] < 0 {
                    0.0
                } else {
                    freq[a * 20 + ALPHA_CONVERT[b_big] as usize] / sum
                };
            }
            std_freq[a_big * std_alphsize + E_BCHAR] =
                std_freq[a_big * std_alphsize + E_DCHAR] + std_freq[a_big * std_alphsize + E_NCHAR];
            std_freq[a_big * std_alphsize + E_ZCHAR] =
                std_freq[a_big * std_alphsize + E_ECHAR] + std_freq[a_big * std_alphsize + E_QCHAR];
            if std_alphsize > E_JCHAR {
                std_freq[a_big * std_alphsize + E_JCHAR] = std_freq[a_big * std_alphsize + E_ICHAR]
                    + std_freq[a_big * std_alphsize + E_LCHAR];
            }
        }
    }
    std_freq.copy_within(
        E_DCHAR * std_alphsize..(E_DCHAR + 1) * std_alphsize,
        E_BCHAR * std_alphsize,
    );
    {
        let (left, right) = std_freq.split_at_mut(E_NCHAR * std_alphsize);
        nlm_add_vectors(
            &mut left[E_BCHAR * std_alphsize..(E_BCHAR + 1) * std_alphsize],
            std_alphsize,
            1.0,
            &right[..std_alphsize],
        );
    }
    std_freq.copy_within(
        E_ECHAR * std_alphsize..(E_ECHAR + 1) * std_alphsize,
        E_ZCHAR * std_alphsize,
    );
    {
        let (left, right) = std_freq.split_at_mut(E_ZCHAR * std_alphsize);
        nlm_add_vectors(
            &mut right[..std_alphsize],
            std_alphsize,
            1.0,
            &left[E_QCHAR * std_alphsize..(E_QCHAR + 1) * std_alphsize],
        );
    }
    if std_alphsize > E_JCHAR {
        std_freq.copy_within(
            E_ICHAR * std_alphsize..(E_ICHAR + 1) * std_alphsize,
            E_JCHAR * std_alphsize,
        );
        let (left, right) = std_freq.split_at_mut(E_JCHAR * std_alphsize);
        nlm_add_vectors(
            &mut right[..std_alphsize],
            std_alphsize,
            1.0,
            &left[E_LCHAR * std_alphsize..(E_LCHAR + 1) * std_alphsize],
        );
    }
}
pub fn s_calc_avg_score(m: &[f64], alphsize: usize, inc_m: usize, probs: &[f64]) -> f64 {
    let mut score = 0.0;
    for j in 0..alphsize {
        if ALPHA_CONVERT[j] >= 0 {
            score += m[j * inc_m] * probs[j];
        }
    }
    score
}
pub fn s_calc_x_score(m: &[f64], alphsize: usize, inc_m: usize, probs: &[f64]) -> f64 {
    let avg = s_calc_avg_score(m, alphsize, inc_m, probs);
    if avg < K_MAXIMUM_XSCORE {
        avg
    } else {
        K_MAXIMUM_XSCORE
    }
}
pub fn s_set_xuo_scores(m: &mut [f64], alphsize: usize, row_probs: &[f64], col_probs: &[f64]) {
    let mut score_xx = 0.0;
    for i in 0..alphsize {
        if ALPHA_CONVERT[i] >= 0 {
            let avg_ix = s_calc_avg_score(&m[i * alphsize..], alphsize, 1, col_probs);
            m[i * alphsize + E_XCHAR] = if avg_ix < K_MAXIMUM_XSCORE {
                avg_ix
            } else {
                K_MAXIMUM_XSCORE
            };
            score_xx += avg_ix * row_probs[i];
            m[E_XCHAR * alphsize + i] = s_calc_x_score(&m[i..], alphsize, alphsize, row_probs);
        }
    }
    m[E_XCHAR * alphsize + E_XCHAR] = if score_xx < K_MAXIMUM_XSCORE {
        score_xx
    } else {
        K_MAXIMUM_XSCORE
    };
    m[E_BCHAR * alphsize + E_XCHAR] =
        s_calc_x_score(&m[E_BCHAR * alphsize..], alphsize, 1, col_probs);
    m[E_XCHAR * alphsize + E_BCHAR] = s_calc_x_score(&m[E_BCHAR..], alphsize, alphsize, row_probs);
    m[E_ZCHAR * alphsize + E_XCHAR] =
        s_calc_x_score(&m[E_ZCHAR * alphsize..], alphsize, 1, col_probs);
    m[E_XCHAR * alphsize + E_ZCHAR] = s_calc_x_score(&m[E_ZCHAR..], alphsize, alphsize, row_probs);
    if alphsize > E_JCHAR {
        m[E_JCHAR * alphsize + E_XCHAR] =
            s_calc_x_score(&m[E_JCHAR * alphsize..], alphsize, 1, col_probs);
        m[E_XCHAR * alphsize + E_JCHAR] =
            s_calc_x_score(&m[E_JCHAR..], alphsize, alphsize, row_probs);
    }
    m.copy_within(
        E_CCHAR * alphsize..(E_CCHAR + 1) * alphsize,
        E_SELENOCYSTEINE * alphsize,
    );
    for i in 0..alphsize {
        m[i * alphsize + E_SELENOCYSTEINE] = m[i * alphsize + E_CCHAR];
    }
    if alphsize > E_OCHAR {
        m.copy_within(
            E_XCHAR * alphsize..(E_XCHAR + 1) * alphsize,
            E_OCHAR * alphsize,
        );
        for i in 0..alphsize {
            m[i * alphsize + E_OCHAR] = m[i * alphsize + E_XCHAR];
        }
    }
}
pub fn nint(mut x: f64) -> std::os::raw::c_long {
    x += if x >= 0.0 { 0.5 } else { -0.5 };
    x as std::os::raw::c_long
}
pub fn s_round_score_matrix(matrix: &mut [i32], rows: usize, cols: usize, float_scores: &[f64]) {
    for p in 0..rows {
        for c in 0..cols {
            let i = p * cols + c;
            matrix[i] = if float_scores[i] < i32::MIN as f64 {
                i32::MIN
            } else {
                nint(float_scores[i]) as i32
            };
        }
    }
}
pub fn s_scores_std_alphabet(
    matrix: &mut [i32],
    alphsize: usize,
    target_freq: &[f64],
    start_matrix: &[i32],
    row_prob: &[f64],
    col_prob: &[f64],
    lambda: f64,
) -> i32 {
    let mut row = [0.0; 28];
    let mut col = [0.0; 28];
    let mut scores = vec![0.0; alphsize * alphsize];
    s_unpack_letter_probs(&mut row, alphsize, row_prob);
    s_set_pair_ambig_probs_to_sum(&mut row, alphsize);
    s_unpack_letter_probs(&mut col, alphsize, col_prob);
    s_set_pair_ambig_probs_to_sum(&mut col, alphsize);
    blast_true_aa_to_std_target_freqs(&mut scores, alphsize, target_freq);
    blast_calc_freq_ratios(&mut scores, alphsize, &row, &col);
    blast_freq_ratio_to_score(&mut scores, alphsize, alphsize, lambda);
    s_set_xuo_scores(&mut scores, alphsize, &row, &col);
    s_round_score_matrix(matrix, alphsize, alphsize, &scores);
    for i in 0..alphsize {
        matrix[i * alphsize + E_STOP_CHAR] = start_matrix[i * alphsize + E_STOP_CHAR];
        matrix[E_STOP_CHAR * alphsize + i] = start_matrix[E_STOP_CHAR * alphsize + i];
    }
    0
}
pub fn s_get_score_range(
    obs_min: &mut i32,
    obs_max: &mut i32,
    matrix: &[i32],
    rows: usize,
    cols: usize,
) {
    let mut min_score = 0;
    let mut max_score = 0;
    for i in 0..rows {
        for aa in 0..20 {
            let j = TRUE_CHAR_POSITIONS[aa];
            let score = matrix[i * cols + j];
            if score < min_score && score > COMPO_SCORE_MIN {
                min_score = score;
            }
            if score > max_score {
                max_score = score;
            }
        }
    }
    *obs_min = min_score;
    *obs_max = max_score;
}
pub fn s_get_matrix_score_probs(
    score_prob: &mut Vec<f64>,
    obs_min: &mut i32,
    obs_max: &mut i32,
    matrix: &[i32],
    alphsize: usize,
    subject_prob: &[f64],
    query_prob: &[f64],
) -> i32 {
    s_get_score_range(obs_min, obs_max, matrix, alphsize, alphsize);
    let min_score = *obs_min;
    let range = *obs_max - *obs_min + 1;
    *score_prob = vec![0.0; range as usize];
    for i in 0..alphsize {
        for aa in 0..20 {
            let j = TRUE_CHAR_POSITIONS[aa];
            let score = matrix[i * alphsize + j];
            if score >= min_score {
                score_prob[(score - min_score) as usize] += query_prob[i] * subject_prob[j];
            }
        }
    }
    0
}
pub fn blast_int4_matrix_from_freq(matrix: &mut [i32], alphsize: usize, freq: &[f64], lambda: f64) {
    let mut row = [0.0; 28];
    for i in 0..alphsize {
        row[..alphsize].copy_from_slice(&freq[i * alphsize..(i + 1) * alphsize]);
        blast_freq_ratio_to_score(&mut row, 1, alphsize, lambda);
        s_round_score_matrix(&mut matrix[i * alphsize..], 1, alphsize, &row);
    }
}
pub fn blast_matrix_info_new(rows: usize, cols: usize, position_based: bool) -> BlastMatrixInfo {
    let mut info = BlastMatrixInfo {
        matrix_name: Vec::new(),
        rows,
        cols,
        position_based,
        start_matrix: vec![0; (rows + 1) * cols],
        start_freq_ratios: vec![0.0; (rows + 1) * cols],
        ungapped_lambda: 0.0,
    };
    for i in 0..cols {
        info.start_matrix[rows * cols + i] = COMPO_SCORE_MIN;
        info.start_freq_ratios[rows * cols + i] = COMPO_SCORE_MIN as f64;
    }
    info
}
pub fn s_scale_square_matrix(
    matrix: &mut [i32],
    alphsize: usize,
    freq_ratios: &[f64],
    start_matrix: &[i32],
    row_prob: &[f64],
    col_prob: &[f64],
    lambda: f64,
) -> i32 {
    let mut scores = vec![0.0; alphsize * alphsize];
    for i in 0..alphsize {
        scores[i * alphsize..(i + 1) * alphsize]
            .copy_from_slice(&freq_ratios[i * alphsize..(i + 1) * alphsize]);
    }
    blast_freq_ratio_to_score(&mut scores, alphsize, alphsize, lambda);
    s_set_xuo_scores(&mut scores, alphsize, row_prob, col_prob);
    s_round_score_matrix(matrix, alphsize, alphsize, &scores);
    for i in 0..alphsize {
        matrix[i * alphsize + E_STOP_CHAR] = start_matrix[i * alphsize + E_STOP_CHAR];
        matrix[E_STOP_CHAR * alphsize + i] = start_matrix[E_STOP_CHAR * alphsize + i];
    }
    0
}
pub fn blast_composition_based_stats(
    matrix: &mut [i32],
    lambda_ratio: &mut f64,
    ss: &BlastMatrixInfo,
    query_prob: &[f64],
    res_prob: &[f64],
    calc_lambda: CalcLambda,
    pvalue_adjustment: i32,
) -> i32 {
    assert!(
        !ss.position_based,
        "PSSM composition adjustment is outside default tblastn"
    );
    let mut obs_min = 0;
    let mut obs_max = 0;
    let mut score_array = Vec::new();
    let status = s_get_matrix_score_probs(
        &mut score_array,
        &mut obs_min,
        &mut obs_max,
        &ss.start_matrix,
        ss.cols,
        res_prob,
        query_prob,
    );
    if status != 0 {
        return -1;
    }
    let correct_ungapped_lambda = calc_lambda(&score_array, obs_min, obs_max, ss.ungapped_lambda);
    *lambda_ratio = correct_ungapped_lambda / ss.ungapped_lambda;
    if pvalue_adjustment == 0 {
        *lambda_ratio = if 1.0 < *lambda_ratio {
            1.0
        } else {
            *lambda_ratio
        };
    }
    *lambda_ratio = if *lambda_ratio > LAMBDA_RATIO_LOWER_BOUND {
        *lambda_ratio
    } else {
        LAMBDA_RATIO_LOWER_BOUND
    };
    if *lambda_ratio > 0.0 {
        let scaled_lambda = ss.ungapped_lambda / *lambda_ratio;
        s_scale_square_matrix(
            matrix,
            ss.cols,
            &ss.start_freq_ratios,
            &ss.start_matrix,
            query_prob,
            res_prob,
            scaled_lambda,
        );
    }
    0
}
pub fn blast_read_aa_composition(
    composition: &mut BlastAminoAcidComposition,
    alphsize: usize,
    sequence: &[u8],
    length: usize,
) {
    let mut num_true = 0;
    for i in 0..alphsize {
        composition.prob[i] = 0.0;
    }
    for i in 0..length {
        let code = sequence[i] as usize;
        if ALPHA_CONVERT[code] >= 0 || code == E_SELENOCYSTEINE {
            composition.prob[code] += 1.0;
            num_true += 1;
        }
    }
    if composition.prob[E_SELENOCYSTEINE] > 0.0 {
        composition.prob[E_CCHAR] += composition.prob[E_SELENOCYSTEINE];
        composition.prob[E_SELENOCYSTEINE] = 0.0;
    }
    composition.num_true_amino_acids = num_true;
    if num_true > 0 {
        for i in 0..alphsize {
            composition.prob[i] /= num_true as f64;
        }
    }
}
pub fn blast_get_composition_range(
    left: &mut i32,
    right: &mut i32,
    subject_data: &[u8],
    length: i32,
    start: i32,
    finish: i32,
) {
    let mut l = start;
    let mut i = l;
    while i > 0 {
        if subject_data[(i - 1) as usize] as usize == E_STOP_CHAR {
            if i + K_COMPOSITION_MARGIN < l {
                l = i + K_COMPOSITION_MARGIN;
            }
            break;
        }
        i -= 1;
    }
    if i == 0 {
        l = 0;
    }
    let mut r = finish;
    i = r;
    while i < length {
        if subject_data[i as usize] as usize == E_STOP_CHAR {
            if i - K_COMPOSITION_MARGIN > r {
                r = i - K_COMPOSITION_MARGIN;
            }
            break;
        }
        i += 1;
    }
    if i == length {
        r = length;
    }
    *left = l;
    *right = r;
}
pub fn blast_composition_workspace_new() -> BlastCompositionWorkspace {
    BlastCompositionWorkspace {
        first_standard_freq: vec![0.0; 20],
        second_standard_freq: vec![0.0; 20],
        mat_final: vec![0.0; 400],
        mat_b: vec![0.0; 400],
    }
}
pub fn blast_composition_workspace_init(
    workspace: &mut BlastCompositionWorkspace,
    matrix_name: &[u8],
) -> i32 {
    if blast_get_joint_probs_for_matrix(
        &mut workspace.mat_b,
        &mut workspace.first_standard_freq,
        &mut workspace.second_standard_freq,
        matrix_name,
    ) == 0
    {
        0
    } else {
        -1
    }
}
pub fn blast_composition_matrix_adj(
    matrix: &mut [i32],
    alphsize: usize,
    matrix_adjust_rule: EMatrixAdjustRule,
    length1: i32,
    length2: i32,
    stdaa_row_probs: &[f64],
    stdaa_col_probs: &[f64],
    pseudocounts: i32,
    specified_re: f64,
    workspace: &mut BlastCompositionWorkspace,
    matrix_info: &BlastMatrixInfo,
) -> i32 {
    let mut row_probs = [0.0; 20];
    let mut col_probs = [0.0; 20];
    s_gather_letter_probs(&mut row_probs, stdaa_row_probs, alphsize);
    s_gather_letter_probs(&mut col_probs, stdaa_col_probs, alphsize);
    let desired_re = match matrix_adjust_rule {
        EMatrixAdjustRule::UnconstrainedRelEntropy => 0.0,
        EMatrixAdjustRule::UserSpecifiedRelEntropy => specified_re,
        _ => panic!("non-default relative entropy rule is outside mode2 tblastn"),
    };
    blast_apply_pseudocounts(
        &mut row_probs,
        length1,
        &workspace.first_standard_freq,
        pseudocounts,
    );
    blast_apply_pseudocounts(
        &mut col_probs,
        length2,
        &workspace.second_standard_freq,
        pseudocounts,
    );
    let mut iteration_count = 0;
    let status = blast_optimize_target_frequencies(
        &mut workspace.mat_final,
        20,
        &mut iteration_count,
        &workspace.mat_b,
        &row_probs,
        &col_probs,
        desired_re > 0.0,
        desired_re,
        K_COMPO_ADJUST_ERR_TOLERANCE,
        K_COMPO_ADJUST_ITERATION_LIMIT,
    );
    if status != 0 {
        return status;
    }
    s_scores_std_alphabet(
        matrix,
        alphsize,
        &workspace.mat_final,
        &matrix_info.start_matrix,
        &row_probs,
        &col_probs,
        matrix_info.ungapped_lambda,
    )
}
pub fn blast_adjust_scores(
    matrix: &mut [i32],
    query_composition: &BlastAminoAcidComposition,
    query_length: i32,
    subject_composition: &BlastAminoAcidComposition,
    subject_length: i32,
    matrix_info: &BlastMatrixInfo,
    mode: ECompoAdjustModes,
    re_pseudocounts: i32,
    workspace: &mut BlastCompositionWorkspace,
    matrix_adjust_rule: &mut EMatrixAdjustRule,
    calc_lambda: CalcLambda,
    _pvalue_for_this_pair: &mut f64,
    composition_test_index: i32,
    ratio_to_pass_back: &mut f64,
) -> i32 {
    assert!(
        !matrix_info.position_based,
        "PSSM composition adjustment is outside default tblastn"
    );
    assert!(
        composition_test_index <= 0,
        "unified composition p-values are outside default mode2"
    );
    assert!(
        mode != ECompoAdjustModes::NoCompositionBasedStats,
        "mode0 bypasses this function in original pipeline"
    );
    let alphsize = matrix_info.cols;
    if query_composition.num_true_amino_acids == 0 || subject_composition.num_true_amino_acids == 0
    {
        return 1;
    }
    let mut query_probs = [0.0; 20];
    let mut subject_probs = [0.0; 20];
    if mode != ECompoAdjustModes::CompositionBasedStats {
        s_gather_letter_probs(&mut query_probs, &query_composition.prob, alphsize);
        s_gather_letter_probs(&mut subject_probs, &subject_composition.prob, alphsize);
    }
    *matrix_adjust_rule = if mode == ECompoAdjustModes::CompositionBasedStats {
        EMatrixAdjustRule::CompoScaleOldMatrix
    } else {
        blast_choose_matrix_adjust_rule(
            query_length,
            subject_length,
            &query_probs,
            &subject_probs,
            &matrix_info.matrix_name,
            mode,
        )
    };
    if *matrix_adjust_rule != EMatrixAdjustRule::CompoScaleOldMatrix {
        let status = blast_composition_matrix_adj(
            matrix,
            alphsize,
            *matrix_adjust_rule,
            query_composition.num_true_amino_acids,
            subject_composition.num_true_amino_acids,
            &query_composition.prob,
            &subject_composition.prob,
            re_pseudocounts,
            K_FIXED_RE_BLOSUM62,
            workspace,
            matrix_info,
        );
        *ratio_to_pass_back = 1.0;
        if status <= 0 {
            return status;
        }
    }
    *matrix_adjust_rule = EMatrixAdjustRule::CompoScaleOldMatrix;
    blast_composition_based_stats(
        matrix,
        ratio_to_pass_back,
        matrix_info,
        &query_composition.prob,
        &subject_composition.prob,
        calc_lambda,
        if composition_test_index > 0 { 1 } else { 0 },
    )
}

/// Existing score-frequency storage is owned; the original stack-only borrowed
/// callback therefore copies its short probability array into that storage.
pub fn s_calc_lambda(probs: &[f64], min_score: i32, max_score: i32, lambda0: f64) -> f64 {
    let score_range = max_score - min_score + 1;
    let mut avg = 0.0;
    for i in 0..score_range {
        avg += (min_score + i) as f64 * probs[i as usize];
    }
    let freq = BlastScoreFreq {
        score_min: min_score,
        score_max: max_score,
        obs_min: min_score,
        obs_max: max_score,
        score_avg: avg,
        sprob: probs[..score_range as usize].to_vec(),
    };
    blast_karlin_lambda_nr(&freq, lambda0)
}
#[derive(Debug)]
pub struct SFreqRatios {
    pub data: Vec<f64>,
    pub bit_scale_factor: i32,
}
pub fn _psi_matrix_frequency_ratios_new(matrix_name: &[u8]) -> Option<SFreqRatios> {
    let n = matrix_name
        .iter()
        .position(|v| *v == 0)
        .unwrap_or(matrix_name.len());
    if matrix_name[..n].eq_ignore_ascii_case(b"BLOSUM62")
        || matrix_name[..n].eq_ignore_ascii_case(b"BLOSUM62_20")
    {
        Some(SFreqRatios {
            data: BLOSUM62_FREQRATIOS.to_vec(),
            bit_scale_factor: 2,
        })
    } else {
        None
    }
}
pub fn s_get_start_freq_ratios(return_ratios: &mut [f64], matrix_name: &[u8]) -> i32 {
    let Some(ratios) = _psi_matrix_frequency_ratios_new(matrix_name) else {
        return -1;
    };
    for i in 0..28 {
        for j in 0..28 {
            return_ratios[i * 28 + j] = ratios.data[i * 28 + j];
        }
    }
    0
}
pub fn s_matrix_info_init(
    info: &mut BlastMatrixInfo,
    ideal_ungapped_lambda: f64,
    scale_factor: f64,
    matrix_name: &[u8],
) -> i32 {
    assert!(
        !info.position_based,
        "PSSM matrix initialization is outside default tblastn"
    );
    let len = matrix_name
        .iter()
        .position(|v| *v == 0)
        .unwrap_or(matrix_name.len());
    info.matrix_name = matrix_name[..len].to_vec();
    info.ungapped_lambda = ideal_ungapped_lambda / scale_factor;
    let status = s_get_start_freq_ratios(&mut info.start_freq_ratios, matrix_name);
    if status == 0 {
        blast_int4_matrix_from_freq(
            &mut info.start_matrix,
            info.cols,
            &info.start_freq_ratios,
            info.ungapped_lambda,
        );
    }
    status
}
