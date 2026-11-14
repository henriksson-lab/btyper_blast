// NCBI BLAST+ 2.17.0+, selected blast_stat.c function.
// Authors: NCBI BLAST developers. See LICENSE-NCBI.
use super::blast_encoding::{BLASTNA_SIZE, BLASTNA_TO_NCBI4NA};
use super::ncbi_math::blast_nint;

pub type NucleotideScoreMatrix = [[i32; BLASTNA_SIZE]; BLASTNA_SIZE];

/// Matches C `BlastScoreBlkNuclMatrixCreate`. Only the fields read/written by
/// that function are passed: reward, penalty and a fixed 16-by-16 matrix.
pub fn blast_score_blk_nucl_matrix_create(
    reward: i32,
    penalty: i32,
    matrix: &mut NucleotideScoreMatrix,
) -> i16 {
    let mut degeneracy = [0_i16; BLASTNA_SIZE + 1];
    let k_number_non_ambig_bp = 4;
    for index1 in 0..BLASTNA_SIZE {
        for index2 in 0..BLASTNA_SIZE {
            matrix[index1][index2] = 0;
        }
    }
    for index1 in 0..k_number_non_ambig_bp {
        degeneracy[index1] = 1;
    }
    for index1 in k_number_non_ambig_bp..BLASTNA_SIZE {
        let mut degen = 0;
        for index2 in 0..k_number_non_ambig_bp {
            if BLASTNA_TO_NCBI4NA[index1] & BLASTNA_TO_NCBI4NA[index2] != 0 {
                degen += 1;
            }
        }
        degeneracy[index1] = degen;
    }
    for index1 in 0..BLASTNA_SIZE {
        for index2 in index1..BLASTNA_SIZE {
            if BLASTNA_TO_NCBI4NA[index1] & BLASTNA_TO_NCBI4NA[index2] != 0 {
                matrix[index1][index2] = blast_nint(
                    f64::from((i32::from(degeneracy[index2]) - 1) * penalty + reward)
                        / f64::from(degeneracy[index2]),
                ) as i32;
                if index1 != index2 {
                    matrix[index2][index1] = matrix[index1][index2];
                }
            } else {
                matrix[index1][index2] = penalty;
                matrix[index2][index1] = penalty;
            }
        }
    }
    for index1 in 0..BLASTNA_SIZE {
        matrix[BLASTNA_SIZE - 1][index1] = i32::MIN / 2;
    }
    for index1 in 0..BLASTNA_SIZE {
        matrix[index1][BLASTNA_SIZE - 1] = i32::MIN / 2;
    }
    0
}

use super::ncbi_math::{blast_expm1, blast_gcd, blast_powi};

pub const BLAST_SCORE_MIN: i32 = i16::MIN as i32;
pub const BLAST_SCORE_MAX: i32 = i16::MAX as i32;
pub const BLAST_SCORE_RANGE_MAX: i32 = BLAST_SCORE_MAX - BLAST_SCORE_MIN;
pub const BLAST_KARLIN_K_SUMLIMIT_DEFAULT: f64 = 0.0001;
pub const BLAST_KARLIN_K_ITER_MAX: i32 = 100;
pub const BLAST_KARLIN_LAMBDA_ACCURACY_DEFAULT: f64 = 1.0e-5;
pub const BLAST_KARLIN_LAMBDA_ITER_DEFAULT: i32 = 17;
pub const BLAST_KARLIN_LAMBDA0_DEFAULT: f64 = 0.5;

#[derive(Clone, Copy, Debug, Default, PartialEq)]
#[repr(C)]
pub struct BlastKarlinBlk {
    pub lambda: f64,
    pub k: f64,
    pub log_k: f64,
    pub h: f64,
    pub param_c: f64,
}
#[derive(Clone, Debug)]
pub struct BlastScoreFreq {
    pub score_min: i32,
    pub score_max: i32,
    pub obs_min: i32,
    pub obs_max: i32,
    pub score_avg: f64,
    /// Original sprob0 storage; score i is stored at i-score_min.
    pub sprob: Vec<f64>,
}
impl BlastScoreFreq {
    pub(crate) fn prob(&self, score: i32) -> f64 {
        self.sprob[(score - self.score_min) as usize]
    }
}

pub fn blast_score_chk(lo: i32, hi: i32) -> i16 {
    if lo >= 0 || hi <= 0 || lo < BLAST_SCORE_MIN || hi > BLAST_SCORE_MAX {
        return 1;
    }
    if hi - lo > BLAST_SCORE_RANGE_MAX {
        return 1;
    }
    0
}
pub fn blast_score_freq_new(score_min: i32, score_max: i32) -> Option<BlastScoreFreq> {
    if blast_score_chk(score_min, score_max) != 0 {
        return None;
    }
    Some(BlastScoreFreq {
        score_min,
        score_max,
        obs_min: 0,
        obs_max: 0,
        score_avg: 0.0,
        sprob: vec![0.0; (score_max - score_min + 1) as usize],
    })
}
pub fn blast_karlin_blk_new() -> BlastKarlinBlk {
    BlastKarlinBlk::default()
}
pub fn blast_karlin_blk_copy(
    to: Option<&mut BlastKarlinBlk>,
    from: Option<&BlastKarlinBlk>,
) -> i16 {
    let (Some(to), Some(from)) = (to, from) else {
        return -1;
    };
    to.lambda = from.lambda;
    to.k = from.k;
    to.log_k = from.log_k;
    to.h = from.h;
    to.param_c = from.param_c;
    0
}

pub fn nlm_karlin_lambda_nr(
    sfp: &BlastScoreFreq,
    d: i32,
    low: i32,
    high: i32,
    lambda0: f64,
    tolx: f64,
    itmax: i32,
    max_newton: i32,
    itn: &mut i32,
) -> f64 {
    assert!(d > 0);
    let x0 = (-lambda0).exp();
    let mut x = if 0.0 < x0 && x0 < 1.0 { x0 } else { 0.5 };
    let (mut a, mut b) = (0.0, 1.0);
    let mut f = 4.0;
    let mut is_newton = false;
    let mut k = 0;
    while k < itmax {
        let fold = f;
        let was_newton = is_newton;
        is_newton = false;
        let mut g = 0.0;
        f = sfp.prob(low);
        let mut i = low + d;
        while i < 0 {
            g = x * g + f;
            f = f * x + sfp.prob(i);
            i += d;
        }
        g = x * g + f;
        f = f * x + sfp.prob(0) - 1.0;
        i = d;
        while i <= high {
            g = x * g + f;
            f = f * x + sfp.prob(i);
            i += d;
        }
        if f > 0.0 {
            a = x;
        } else if f < 0.0 {
            b = x;
        } else {
            break;
        }
        if b - a < 2.0 * a * (1.0 - b) * tolx {
            x = (a + b) / 2.0;
            break;
        }
        if k >= max_newton || (was_newton && f.abs() > 0.9 * fold.abs()) || g >= 0.0 {
            x = (a + b) / 2.0;
        } else {
            let p = -f / g;
            let y = x + p;
            if y <= a || y >= b {
                x = (a + b) / 2.0;
            } else {
                is_newton = true;
                x = y;
                if p.abs() < tolx * x * (1.0 - x) {
                    break;
                }
            }
        }
        k += 1;
    }
    *itn = k;
    -x.ln() / d as f64
}
pub fn blast_karlin_lambda_nr(sfp: &BlastScoreFreq, initial_lambda_guess: f64) -> f64 {
    let low = sfp.obs_min;
    let high = sfp.obs_max;
    if sfp.score_avg >= 0.0 {
        return -1.0;
    }
    if blast_score_chk(low, high) != 0 {
        return -1.0;
    }
    let mut d = -low;
    let mut i = 1;
    while i <= high - low && d > 1 {
        if sfp.prob(i + low) != 0.0 {
            d = blast_gcd(d, i);
        }
        i += 1;
    }
    let mut itn = 0;
    nlm_karlin_lambda_nr(
        sfp,
        d,
        low,
        high,
        initial_lambda_guess,
        BLAST_KARLIN_LAMBDA_ACCURACY_DEFAULT,
        20,
        20 + BLAST_KARLIN_LAMBDA_ITER_DEFAULT,
        &mut itn,
    )
}
pub fn blast_karlin_lto_h(sfp: &BlastScoreFreq, lambda: f64) -> f64 {
    let low = sfp.obs_min;
    let high = sfp.obs_max;
    if lambda < 0.0 {
        return -1.0;
    }
    if blast_score_chk(low, high) != 0 {
        return -1.0;
    }
    let etonlam = (-lambda).exp();
    let mut sum = low as f64 * sfp.prob(low);
    for score in low + 1..=high {
        sum = score as f64 * sfp.prob(score) + etonlam * sum;
    }
    let scale = blast_powi(etonlam, high);
    if scale > 0.0 {
        lambda * sum / scale
    } else {
        lambda * (lambda * high as f64 + sum.ln()).exp()
    }
}
pub fn blast_karlin_l_hto_k(sfp: &BlastScoreFreq, mut lambda: f64, h: f64) -> f64 {
    if lambda <= 0.0 || h <= 0.0 {
        return -1.0;
    }
    if sfp.score_avg >= 0.0 {
        return -1.0;
    }
    let mut low = sfp.obs_min;
    let mut high = sfp.obs_max;
    let mut range = high - low;
    let prob_start_low = low;
    let mut divisor = -low;
    let mut i = 1;
    while i <= range && divisor > 1 {
        if sfp.prob(prob_start_low + i) != 0.0 {
            divisor = blast_gcd(divisor, i);
        }
        i += 1;
    }
    high /= divisor;
    low /= divisor;
    lambda *= divisor as f64;
    range = high - low;
    let mut first_term_closed_form = h / lambda;
    let exp_minus_lambda = (-lambda).exp();
    if low == -1 && high == 1 {
        return (sfp.prob(low * divisor) - sfp.prob(high * divisor))
            * (sfp.prob(low * divisor) - sfp.prob(high * divisor))
            / sfp.prob(low * divisor);
    }
    if low == -1 || high == 1 {
        if high != 1 {
            let score_avg = sfp.score_avg / divisor as f64;
            first_term_closed_form = (score_avg * score_avg) / first_term_closed_form;
        }
        return first_term_closed_form * (1.0 - exp_minus_lambda);
    }
    let mut probabilities = vec![0.0; (BLAST_KARLIN_K_ITER_MAX * range + 1) as usize];
    let mut outer_sum = 0.0;
    let (mut low_alignment_score, mut high_alignment_score) = (0, 0);
    probabilities[0] = 1.0;
    let mut inner_sum = 1.0;
    let mut iter_counter = 0;
    while iter_counter < BLAST_KARLIN_K_ITER_MAX && inner_sum > BLAST_KARLIN_K_SUMLIMIT_DEFAULT {
        let (mut first, mut last) = (range, range);
        low_alignment_score += low;
        high_alignment_score += high;
        let mut p = high_alignment_score - low_alignment_score;
        while p >= 0 {
            let mut ptr1 = p - first;
            let ptr1e = p - last;
            let mut ptr2 = first;
            inner_sum = 0.0;
            while ptr1 >= ptr1e {
                inner_sum += probabilities[ptr1 as usize] * sfp.prob(prob_start_low + ptr2);
                ptr1 -= 1;
                ptr2 += 1;
            }
            if first != 0 {
                first -= 1;
            }
            if p <= range {
                last -= 1;
            }
            probabilities[p as usize] = inner_sum;
            p -= 1;
        }
        p += 1;
        inner_sum = probabilities[p as usize];
        let mut i = low_alignment_score + 1;
        while i < 0 {
            p += 1;
            inner_sum = probabilities[p as usize] + inner_sum * exp_minus_lambda;
            i += 1;
        }
        inner_sum *= exp_minus_lambda;
        while i <= high_alignment_score {
            p += 1;
            inner_sum += probabilities[p as usize];
            i += 1;
        }
        iter_counter += 1;
        inner_sum /= iter_counter as f64;
        outer_sum += inner_sum;
    }
    -(-2.0 * outer_sum).exp() / (first_term_closed_form * blast_expm1(-lambda))
}
pub fn blast_karlin_blk_ungapped_calc(
    kbp: Option<&mut BlastKarlinBlk>,
    sfp: Option<&BlastScoreFreq>,
) -> i16 {
    let (Some(kbp), Some(sfp)) = (kbp, sfp) else {
        return 1;
    };
    kbp.lambda = blast_karlin_lambda_nr(sfp, BLAST_KARLIN_LAMBDA0_DEFAULT);
    if !(kbp.lambda < 0.0) {
        kbp.h = blast_karlin_lto_h(sfp, kbp.lambda);
        if !(kbp.h < 0.0) {
            kbp.k = blast_karlin_l_hto_k(sfp, kbp.lambda, kbp.h);
            if !(kbp.k < 0.0) {
                kbp.log_k = kbp.k.ln();
                return 0;
            }
        }
    }
    kbp.lambda = -1.0;
    kbp.h = -1.0;
    kbp.k = -1.0;
    kbp.log_k = f64::INFINITY;
    1
}

use super::blast_encoding::{AMINOACID_TO_NCBISTDAA, BLASTAA_SEQ_CODE};

/// Score-block fields consumed by composition and ideal statistics. Matrix
/// indexes retain absolute alphabet codes; frequency arrays start at
/// alphabet_start, replacing the original shifted C pointers.
#[derive(Clone, Debug, Default)]
pub struct BlastScoreBlk {
    pub alphabet_code: u8,
    pub alphabet_size: i16,
    pub alphabet_start: i16,
    pub protein_alphabet: bool,
    pub ambiguous_res: Vec<u8>,
    pub loscore: i32,
    pub hiscore: i32,
    pub matrix: Vec<Vec<i32>>,
    pub kbp_ideal: Option<BlastKarlinBlk>,
    pub sfp: Vec<Option<BlastScoreFreq>>,
    pub kbp_std: Vec<Option<BlastKarlinBlk>>,
    pub kbp_psi: Vec<Option<BlastKarlinBlk>>,
    pub kbp_is_psi: bool,
    pub kbp_gap: Option<Vec<Option<BlastKarlinBlk>>>,
    pub reward: i32,
    pub penalty: i32,
    pub scale_factor: f64,
    pub matrix_only_scoring: bool,
    pub round_down: bool,
    pub gbp: Option<super::blast_stat_spouge::BlastGumbelBlk>,
}
#[derive(Clone, Debug)]
pub struct BlastResFreq {
    pub alphabet_code: u8,
    pub prob: Vec<f64>,
}
#[derive(Clone, Debug)]
pub struct BlastResComp {
    pub alphabet_code: u8,
    pub comp: Vec<i32>,
}

pub const ROBINSON_PROB: [(u8, f64); 20] = [
    (b'A', 78.05),
    (b'C', 19.25),
    (b'D', 53.64),
    (b'E', 62.95),
    (b'F', 38.56),
    (b'G', 73.77),
    (b'H', 21.99),
    (b'I', 51.42),
    (b'K', 57.44),
    (b'L', 90.19),
    (b'M', 22.43),
    (b'N', 44.87),
    (b'P', 52.03),
    (b'Q', 42.64),
    (b'R', 51.29),
    (b'S', 71.20),
    (b'T', 58.41),
    (b'V', 64.41),
    (b'W', 13.30),
    (b'Y', 32.16),
];
pub const NT_PROB: [(u8, f64); 4] = [(b'A', 25.0), (b'C', 25.0), (b'G', 25.0), (b'T', 25.0)];

pub fn blast_res_freq_new(sbp: Option<&BlastScoreBlk>) -> Option<BlastResFreq> {
    let sbp = sbp?;
    Some(BlastResFreq {
        alphabet_code: sbp.alphabet_code,
        prob: vec![0.0; sbp.alphabet_size as usize],
    })
}
pub fn blast_res_freq_normalize(sbp: &BlastScoreBlk, rfp: &mut BlastResFreq, norm: f64) -> i16 {
    if norm == 0.0 {
        return 1;
    }
    let mut sum = 0.0;
    let stop = sbp.alphabet_start + sbp.alphabet_size;
    for index in sbp.alphabet_start..stop {
        let p = rfp.prob[(index - sbp.alphabet_start) as usize];
        if p < 0.0 {
            return 1;
        }
        sum += p;
    }
    if sum <= 0.0 {
        return 0;
    }
    for index in sbp.alphabet_start..stop {
        rfp.prob[(index - sbp.alphabet_start) as usize] /= sum;
        rfp.prob[(index - sbp.alphabet_start) as usize] *= norm;
    }
    0
}
pub fn blast_get_std_alphabet(alphabet_code: u8, residues: &mut [u8]) -> i16 {
    if residues.len() < ROBINSON_PROB.len() {
        return -2;
    }
    for index in 0..ROBINSON_PROB.len() {
        residues[index] = if alphabet_code as usize == BLASTAA_SEQ_CODE {
            AMINOACID_TO_NCBISTDAA[ROBINSON_PROB[index].0.to_ascii_uppercase() as usize]
        } else {
            ROBINSON_PROB[index].0
        };
    }
    ROBINSON_PROB.len() as i16
}
pub fn blast_res_freq_std_comp(sbp: &BlastScoreBlk, rfp: &mut BlastResFreq) -> i16 {
    if sbp.protein_alphabet {
        let mut residues = [0; 20];
        let retval = blast_get_std_alphabet(sbp.alphabet_code, &mut residues);
        if retval < 1 {
            return retval;
        }
        for index in 0..ROBINSON_PROB.len() {
            rfp.prob[(residues[index] as i16 - sbp.alphabet_start) as usize] =
                ROBINSON_PROB[index].1;
        }
    } else {
        for index in 0..NT_PROB.len() {
            rfp.prob[(index as i16 - sbp.alphabet_start) as usize] = NT_PROB[index].1;
        }
    }
    blast_res_freq_normalize(sbp, rfp, 1.0);
    0
}
pub fn blast_res_comp_new(sbp: &BlastScoreBlk) -> BlastResComp {
    BlastResComp {
        alphabet_code: sbp.alphabet_code,
        comp: vec![0; sbp.alphabet_size as usize],
    }
}
pub fn blast_res_comp_str(
    sbp: Option<&BlastScoreBlk>,
    rcp: Option<&mut BlastResComp>,
    sequence: Option<&[u8]>,
) -> i16 {
    let (Some(sbp), Some(rcp), Some(sequence)) = (sbp, rcp, sequence) else {
        return 1;
    };
    if rcp.alphabet_code != sbp.alphabet_code {
        return 1;
    }
    let mask = if sbp.protein_alphabet { 0xff } else { 0x0f };
    for index in 0..sbp.alphabet_size {
        rcp.comp[index as usize] = 0;
    }
    for residue in sequence {
        let slot = ((residue & mask) as i16 - sbp.alphabet_start) as usize;
        rcp.comp[slot] = rcp.comp[slot].wrapping_add(1);
    }
    for residue in &sbp.ambiguous_res {
        rcp.comp[(*residue as i16 - sbp.alphabet_start) as usize] = 0;
    }
    0
}
pub fn blast_res_freq_clr(sbp: &BlastScoreBlk, rfp: &mut BlastResFreq) -> i16 {
    let stop = sbp.alphabet_start + sbp.alphabet_size;
    for index in sbp.alphabet_start..stop {
        rfp.prob[(index - sbp.alphabet_start) as usize] = 0.0;
    }
    0
}
pub fn blast_res_freq_res_comp(
    sbp: &BlastScoreBlk,
    rfp: Option<&mut BlastResFreq>,
    rcp: Option<&BlastResComp>,
) -> i16 {
    let (Some(rfp), Some(rcp)) = (rfp, rcp) else {
        return 1;
    };
    if rfp.alphabet_code != rcp.alphabet_code {
        return 1;
    }
    let stop = sbp.alphabet_start + sbp.alphabet_size;
    let mut sum = 0.0;
    for index in sbp.alphabet_start..stop {
        sum += rcp.comp[(index - sbp.alphabet_start) as usize] as f64;
    }
    if sum == 0.0 {
        blast_res_freq_clr(sbp, rfp);
        return 0;
    }
    for index in sbp.alphabet_start..stop {
        rfp.prob[(index - sbp.alphabet_start) as usize] =
            rcp.comp[(index - sbp.alphabet_start) as usize] as f64 / sum;
    }
    0
}
pub fn blast_res_freq_string(sbp: &BlastScoreBlk, rfp: &mut BlastResFreq, sequence: &[u8]) -> i16 {
    let mut rcp = blast_res_comp_new(sbp);
    blast_res_comp_str(Some(sbp), Some(&mut rcp), Some(sequence));
    blast_res_freq_res_comp(sbp, Some(rfp), Some(&rcp));
    0
}
pub fn blast_score_freq_calc(
    sbp: Option<&BlastScoreBlk>,
    sfp: Option<&mut BlastScoreFreq>,
    rfp1: &BlastResFreq,
    rfp2: &BlastResFreq,
) -> i16 {
    let (Some(sbp), Some(sfp)) = (sbp, sfp) else {
        return 1;
    };
    if sbp.loscore < sfp.score_min || sbp.hiscore > sfp.score_max {
        return 1;
    }
    for score in sfp.score_min..=sfp.score_max {
        sfp.sprob[(score - sfp.score_min) as usize] = 0.0;
    }
    let start = sbp.alphabet_start;
    let end = start + sbp.alphabet_size;
    for index1 in start..end {
        for index2 in start..end {
            let score = sbp.matrix[index1 as usize][index2 as usize];
            if score >= sbp.loscore {
                sfp.sprob[(score - sfp.score_min) as usize] +=
                    rfp1.prob[(index1 - start) as usize] * rfp2.prob[(index2 - start) as usize];
            }
        }
    }
    let mut score_sum = 0.0;
    let mut obs_min = BLAST_SCORE_MIN;
    let mut obs_max = BLAST_SCORE_MIN;
    for score in sfp.score_min..=sfp.score_max {
        if sfp.prob(score) > 0.0 {
            score_sum += sfp.prob(score);
            obs_max = score;
            if obs_min == BLAST_SCORE_MIN {
                obs_min = score;
            }
        }
    }
    sfp.obs_min = obs_min;
    sfp.obs_max = obs_max;
    let mut score_avg = 0.0;
    if score_sum > 0.0001 || score_sum < -0.0001 {
        for score in obs_min..=obs_max {
            sfp.sprob[(score - sfp.score_min) as usize] /= score_sum;
            score_avg += score as f64 * sfp.prob(score);
        }
    }
    sfp.score_avg = score_avg;
    0
}
pub fn blast_score_blk_kbp_ideal_calc(sbp: Option<&mut BlastScoreBlk>) -> i16 {
    let Some(sbp) = sbp else {
        return 1;
    };
    let mut stdrfp = blast_res_freq_new(Some(sbp)).unwrap();
    blast_res_freq_std_comp(sbp, &mut stdrfp);
    let mut sfp = blast_score_freq_new(sbp.loscore, sbp.hiscore);
    blast_score_freq_calc(Some(sbp), sfp.as_mut(), &stdrfp, &stdrfp);
    let mut kbp = blast_karlin_blk_new();
    blast_karlin_blk_ungapped_calc(Some(&mut kbp), sfp.as_ref());
    sbp.kbp_ideal = Some(kbp);
    0
}

pub fn blast_score_blk_max_score_set(sbp: &mut BlastScoreBlk) -> i16 {
    sbp.loscore = BLAST_SCORE_MAX;
    sbp.hiscore = BLAST_SCORE_MIN;
    for index1 in 0..sbp.alphabet_size {
        for index2 in 0..sbp.alphabet_size {
            let score = sbp.matrix[index1 as usize][index2 as usize];
            if score <= BLAST_SCORE_MIN || score >= BLAST_SCORE_MAX {
                continue;
            }
            if sbp.loscore > score {
                sbp.loscore = score;
            }
            if sbp.hiscore < score {
                sbp.hiscore = score;
            }
        }
    }
    if sbp.loscore < BLAST_SCORE_MIN {
        sbp.loscore = BLAST_SCORE_MIN;
    }
    if sbp.hiscore > BLAST_SCORE_MAX {
        sbp.hiscore = BLAST_SCORE_MAX;
    }
    0
}

impl BlastScoreBlk {
    pub(crate) fn ungapped_kbp(&self) -> &[Option<BlastKarlinBlk>] {
        if self.kbp_is_psi {
            &self.kbp_psi
        } else {
            &self.kbp_std
        }
    }
}

use super::blast_program::{blast_query_is_pssm, blast_query_is_translated, EBlastProgramType};
use super::blast_query_info::BlastQueryInfo;

pub const K_BLAST_ERR_MSG_CANT_CALCULATE_UNGAPPED_KA_PARAMS: &str = "Could not calculate ungapped Karlin-Altschul parameters due to an invalid query sequence or its translation. Please verify the query sequence(s) and/or filtering options";
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StatisticalWarning {
    pub context: i32,
    pub message: &'static str,
}

pub fn blast_score_blk_kbp_ungapped_calc(
    program: EBlastProgramType,
    sbp: &mut BlastScoreBlk,
    query: &[u8],
    query_info: &mut BlastQueryInfo,
    mut warnings: Option<&mut Vec<StatisticalWarning>>,
) -> i16 {
    let mut status = blast_score_blk_kbp_ideal_calc(Some(sbp));
    if status != 0 {
        return status;
    }
    let mut stdrfp = blast_res_freq_new(Some(sbp)).unwrap();
    blast_res_freq_std_comp(sbp, &mut stdrfp);
    let mut rfp = blast_res_freq_new(Some(sbp)).unwrap();
    let check_ideal = matches!(
        program,
        EBlastProgramType::Blastx | EBlastProgramType::Tblastx | EBlastProgramType::RpsTblastn
    );
    let mut valid_context = false;
    let count = (query_info.last_context + 1) as usize;
    if sbp.sfp.len() < count {
        sbp.sfp.resize_with(count, || None);
    }
    if sbp.kbp_std.len() < count {
        sbp.kbp_std.resize(count, None);
    }
    if sbp.kbp_psi.len() < count {
        sbp.kbp_psi.resize(count, None);
    }
    for context in query_info.first_context..=query_info.last_context {
        let c = context as usize;
        if !query_info.contexts[c].is_valid {
            continue;
        }
        let query_length = query_info.contexts[c].query_length as usize;
        let context_offset = query_info.contexts[c].query_offset as usize;
        let buffer = &query[context_offset..context_offset + query_length];
        blast_res_freq_string(sbp, &mut rfp, buffer);
        let mut sfp = blast_score_freq_new(sbp.loscore, sbp.hiscore);
        blast_score_freq_calc(Some(sbp), sfp.as_mut(), &rfp, &stdrfp);
        let mut kbp = blast_karlin_blk_new();
        let loop_status = blast_karlin_blk_ungapped_calc(Some(&mut kbp), sfp.as_ref());
        if loop_status != 0 {
            query_info.contexts[c].is_valid = false;
            sbp.sfp[c] = None;
            sbp.kbp_std[c] = None;
            if !blast_query_is_translated(program) {
                if let Some(out) = warnings.as_deref_mut() {
                    out.push(StatisticalWarning {
                        context,
                        message: K_BLAST_ERR_MSG_CANT_CALCULATE_UNGAPPED_KA_PARAMS,
                    });
                }
            }
            continue;
        }
        if check_ideal && kbp.lambda >= sbp.kbp_ideal.as_ref().unwrap().lambda {
            blast_karlin_blk_copy(Some(&mut kbp), sbp.kbp_ideal.as_ref());
        }
        sbp.sfp[c] = sfp;
        sbp.kbp_std[c] = Some(kbp);
        let mut psi = blast_karlin_blk_new();
        let loop_status = blast_karlin_blk_ungapped_calc(Some(&mut psi), sbp.sfp[c].as_ref());
        if loop_status != 0 {
            query_info.contexts[c].is_valid = false;
            sbp.sfp[c] = None;
            sbp.kbp_std[c] = None;
            sbp.kbp_psi[c] = None;
            continue;
        }
        sbp.kbp_psi[c] = Some(psi);
        valid_context = true;
    }
    if !valid_context {
        if blast_query_is_translated(program) {
            if let Some(out) = warnings {
                out.push(StatisticalWarning {
                    context: -1,
                    message: K_BLAST_ERR_MSG_CANT_CALCULATE_UNGAPPED_KA_PARAMS,
                });
            }
        }
        status = 1;
    }
    sbp.kbp_is_psi = blast_query_is_pssm(program);
    status
}

pub fn blast_karlin_eto_s_simple(mut e: f64, kbp: &BlastKarlinBlk, searchsp: i64) -> i32 {
    if kbp.lambda < 0.0 || kbp.k < 0.0 || kbp.h < 0.0 {
        return BLAST_SCORE_MIN;
    }
    let small = 1.0e-297;
    e = if e > small { e } else { small };
    ((kbp.k * searchsp as f64 / e).ln() / kbp.lambda).ceil() as i32
}
pub fn blast_gap_decay_divisor(decayrate: f64, nsegs: u32) -> f64 {
    (1.0 - decayrate) * blast_powi(decayrate, nsegs.wrapping_sub(1) as i32)
}
pub fn blast_karlin_sto_e_simple(s: i32, kbp: &BlastKarlinBlk, searchsp: i64) -> f64 {
    if kbp.lambda < 0.0 || kbp.k < 0.0 || kbp.h < 0.0 {
        return -1.0;
    }
    searchsp as f64 * (-kbp.lambda * s as f64 + kbp.log_k).exp()
}
pub fn blast_cutoffs(
    s: &mut i32,
    e: &mut f64,
    kbp: &BlastKarlinBlk,
    searchsp: i64,
    dodecay: bool,
    gap_decay_rate: f64,
) -> i16 {
    let mut score = *s;
    let mut expect = *e;
    let mut s_changed = false;
    if kbp.lambda == -1.0 || kbp.k == -1.0 || kbp.h == -1.0 {
        return 1;
    }
    let mut es = 1;
    let esave = expect;
    if expect > 0.0 {
        if dodecay && gap_decay_rate > 0.0 && gap_decay_rate < 1.0 {
            expect *= blast_gap_decay_divisor(gap_decay_rate, 1);
        }
        es = blast_karlin_eto_s_simple(expect, kbp, searchsp);
    }
    if es > score {
        s_changed = true;
        *s = es;
        score = es;
    }
    if esave <= 0.0 || !s_changed {
        expect = blast_karlin_sto_e_simple(score, kbp, searchsp);
        if dodecay && gap_decay_rate > 0.0 && gap_decay_rate < 1.0 {
            expect /= blast_gap_decay_divisor(gap_decay_rate, 1);
        }
        *e = expect;
    }
    0
}
