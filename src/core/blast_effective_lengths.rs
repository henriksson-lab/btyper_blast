// NCBI BLAST+ 2.17.0+, selected blast_stat.c/blast_setup.c length setup.
// Authors: NCBI BLAST developers. See LICENSE-NCBI.
use super::blast_program::{
    blast_program_is_mapping, blast_program_is_phi_blast, blast_subject_is_translated,
    EBlastProgramType,
};
use super::blast_protein_matrix::blast_get_alpha_beta;
use super::blast_query_info::BlastQueryInfo;
use super::blast_stat::BlastScoreBlk;
use super::blast_stat_gapped::blast_get_nucl_alpha_beta;

pub const BLAST_REWARD: i32 = 1;
pub const BLAST_PENALTY: i32 = -3;

#[derive(Clone, Debug, Default)]
pub struct BlastEffectiveLengthsOptions {
    pub db_length: i64,
    pub dbseq_num: i32,
    pub num_searchspaces: i32,
    pub searchsp_eff: Option<Vec<i64>>,
}
#[derive(Clone, Debug, Default)]
pub struct BlastEffectiveLengthsParameters {
    pub options: BlastEffectiveLengthsOptions,
    pub real_db_length: i64,
    pub real_num_seqs: i32,
}
/// Only the scoring fields read by effective-length setup. Standard search
/// constructors must supply these values; Default is the C calloc state.
#[derive(Clone, Copy, Debug, Default)]
pub struct EffectiveLengthScoringOptions {
    pub gapped_calculation: bool,
    pub reward: i32,
    pub penalty: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
}
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct BlastEffectiveLengthWarning {
    pub context_index: i32,
    pub message: String,
}

pub fn blast_effective_lengths_options_new() -> BlastEffectiveLengthsOptions {
    BlastEffectiveLengthsOptions::default()
}
pub fn blast_effective_lengths_options_is_search_space_set(
    options: Option<&BlastEffectiveLengthsOptions>,
) -> bool {
    let Some(options) = options else {
        return false;
    };
    let Some(spaces) = &options.searchsp_eff else {
        return false;
    };
    for i in 0..options.num_searchspaces {
        if spaces[i as usize] != 0 {
            return true;
        }
    }
    false
}
pub fn blast_effective_lengths_parameters_new(
    options: &BlastEffectiveLengthsOptions,
    db_length: i64,
    num_seqs: i32,
) -> BlastEffectiveLengthsParameters {
    BlastEffectiveLengthsParameters {
        options: options.clone(),
        real_db_length: db_length,
        real_num_seqs: num_seqs,
    }
}
pub fn s_get_effective_search_space_for_context(
    options: &BlastEffectiveLengthsOptions,
    context_index: i32,
    blast_message: Option<&mut Vec<BlastEffectiveLengthWarning>>,
) -> i64 {
    let retval;
    if options.num_searchspaces == 0 {
        retval = 0;
    } else if options.num_searchspaces == 1 {
        if context_index != 0 {
            if let Some(messages) = blast_message {
                messages.push(BlastEffectiveLengthWarning {
                    context_index,
                    message: "One search space is being used for multiple sequences".to_owned(),
                });
            }
        }
        retval = options
            .searchsp_eff
            .as_ref()
            .expect("search-space array required")[0];
    } else if options.num_searchspaces > 1 {
        assert!(context_index < options.num_searchspaces);
        retval = options
            .searchsp_eff
            .as_ref()
            .expect("search-space array required")[context_index as usize];
    } else {
        panic!("negative effective search-space count");
    }
    retval
}

pub fn blast_compute_length_adjustment(
    k: f64,
    log_k: f64,
    alpha_d_lambda: f64,
    beta: f64,
    query_length: i32,
    db_length: i64,
    db_num_seqs: i32,
    length_adjustment: &mut i32,
) -> i32 {
    let k_max_iterations = 20;
    let m = query_length as f64;
    let n = db_length as f64;
    let number = db_num_seqs as f64;
    let mut ell_min = 0.0;
    let mut ell_max;
    let mut converged = false;
    let mut ell_next = 0.0;
    {
        let a = number;
        let mb = m * number + n;
        let c = n * m - (if m > n { m } else { n }) / k;
        if c < 0.0 {
            *length_adjustment = 0;
            return 1;
        } else {
            ell_max = 2.0 * c / (mb + (mb * mb - 4.0 * a * c).sqrt());
        }
    }
    for i in 1..=k_max_iterations {
        let ell = ell_next;
        let ss = (m - ell) * (n - number * ell);
        let ell_bar = alpha_d_lambda * (log_k + ss.ln()) + beta;
        if ell_bar >= ell {
            ell_min = ell;
            if ell_bar - ell_min <= 1.0 {
                converged = true;
                break;
            }
            if ell_min == ell_max {
                break;
            }
        } else {
            ell_max = ell;
        }
        if ell_min <= ell_bar && ell_bar <= ell_max {
            ell_next = ell_bar;
        } else {
            ell_next = if i == 1 {
                ell_max
            } else {
                (ell_min + ell_max) / 2.0
            };
        }
    }
    if converged {
        *length_adjustment = ell_min as i32;
        let ell = ell_min.ceil();
        if ell <= ell_max {
            let ss = (m - ell) * (n - number * ell);
            if alpha_d_lambda * (log_k + ss.ln()) + beta >= ell {
                *length_adjustment = ell as i32;
            }
        }
    } else {
        *length_adjustment = ell_min as i32;
    }
    if converged {
        0
    } else {
        1
    }
}

/// The standard non-PSSM gapped alias is supplied by sbp.kbp_gap. The matrix
/// name is explicit in place of the original sbp->name member. Nucleotide
/// setup uses its original alpha/beta tables; protein setup uses the scoped
/// original BLOSUM62 loader supplied by blast_protein_matrix.
pub fn blast_calc_eff_lengths(
    program_number: EBlastProgramType,
    scoring_options: &EffectiveLengthScoringOptions,
    eff_len_params: &BlastEffectiveLengthsParameters,
    sbp: Option<&BlastScoreBlk>,
    query_info: Option<&mut BlastQueryInfo>,
    matrix_name: Option<&[u8]>,
    blast_message: Option<&mut Vec<BlastEffectiveLengthWarning>>,
) -> i16 {
    let (Some(query_info), Some(sbp)) = (query_info, sbp) else {
        return -1;
    };
    let eff_len_options = &eff_len_params.options;
    let (mut alpha, mut beta) = (0.0, 0.0);
    let mut db_length = if eff_len_options.db_length > 0 {
        eff_len_options.db_length
    } else {
        eff_len_params.real_db_length
    };
    if db_length == 0 && !blast_effective_lengths_options_is_search_space_set(Some(eff_len_options))
    {
        return 0;
    }
    if blast_subject_is_translated(program_number) {
        db_length /= 3;
    }
    let db_num_seqs = if eff_len_options.dbseq_num > 0 {
        eff_len_options.dbseq_num
    } else {
        eff_len_params.real_num_seqs
    };
    if blast_program_is_mapping(program_number) {
        for index in query_info.first_context..=query_info.last_context {
            query_info.contexts[index as usize].eff_searchsp = db_length;
        }
        return 0;
    }
    if blast_program_is_phi_blast(program_number) {
        for index in query_info.first_context..=query_info.last_context {
            let adjustment = query_info.contexts[index as usize].length_adjustment;
            let effective_search_space =
                db_length.wrapping_sub(db_num_seqs.wrapping_mul(adjustment) as i64);
            query_info.contexts[index as usize].eff_searchsp = effective_search_space;
        }
        return 0;
    }
    // In the supported non-PSSM programs kbp_gap is the kbp_gap_std array.
    let kbp_ptr = if scoring_options.gapped_calculation {
        sbp.kbp_gap
            .as_deref()
            .expect("gapped statistical blocks required")
    } else {
        sbp.ungapped_kbp()
    };
    let mut messages = blast_message;
    for index in query_info.first_context..=query_info.last_context {
        let mut length_adjustment = 0;
        let mut effective_search_space = s_get_effective_search_space_for_context(
            eff_len_options,
            index,
            messages.as_deref_mut(),
        );
        let kbp = kbp_ptr[index as usize].as_ref();
        let context = &mut query_info.contexts[index as usize];
        let query_length = context.query_length;
        if context.is_valid && query_length > 0 {
            let ungapped = sbp.kbp_std[index as usize]
                .as_ref()
                .expect("valid ungapped context required");
            if program_number == EBlastProgramType::Blastn {
                if scoring_options.reward == 0 && scoring_options.penalty == 0 {
                    blast_get_nucl_alpha_beta(
                        BLAST_REWARD,
                        BLAST_PENALTY,
                        scoring_options.gap_open,
                        scoring_options.gap_extend,
                        ungapped,
                        scoring_options.gapped_calculation,
                        &mut alpha,
                        &mut beta,
                    );
                } else {
                    blast_get_nucl_alpha_beta(
                        scoring_options.reward,
                        scoring_options.penalty,
                        scoring_options.gap_open,
                        scoring_options.gap_extend,
                        ungapped,
                        scoring_options.gapped_calculation,
                        &mut alpha,
                        &mut beta,
                    );
                }
            } else {
                blast_get_alpha_beta(
                    matrix_name,
                    &mut alpha,
                    &mut beta,
                    scoring_options.gapped_calculation,
                    scoring_options.gap_open,
                    scoring_options.gap_extend,
                    ungapped,
                );
            }
            let kbp = kbp.expect("valid context statistical block required");
            blast_compute_length_adjustment(
                kbp.k,
                kbp.log_k,
                alpha / kbp.lambda,
                beta,
                query_length,
                db_length,
                db_num_seqs,
                &mut length_adjustment,
            );
            if effective_search_space == 0 {
                let mut effective_db_length = db_length
                    .wrapping_sub((db_num_seqs as i64).wrapping_mul(length_adjustment as i64));
                if effective_db_length <= 0 {
                    effective_db_length = 1;
                }
                effective_search_space = effective_db_length
                    .wrapping_mul(query_length.wrapping_sub(length_adjustment) as i64);
            }
        }
        context.eff_searchsp = effective_search_space;
        context.length_adjustment = length_adjustment;
    }
    0
}
