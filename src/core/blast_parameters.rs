// NCBI BLAST+ 2.17.0+, selected blast_parameters.h data definitions.
// Authors: NCBI BLAST developers. See LICENSE-NCBI.
use super::blast_options::BlastInitialWordOptions;

#[derive(Clone, Copy, Debug, Default)]
#[repr(C)]
pub struct BlastUngappedCutoffs {
    pub x_dropoff_init: i32,
    pub x_dropoff: i32,
    pub cutoff_score: i32,
    pub reduced_nucl_cutoff_score: i32,
}

/// Fields consumed by the supported nucleotide/protein initial extension and
/// statistical parameter constructors.
#[derive(Debug)]
pub struct BlastInitialWordParameters {
    pub options: BlastInitialWordOptions,
    pub cutoffs: Vec<BlastUngappedCutoffs>,
    pub nucl_score_table: [i32; 256],
    pub matrix_only_scoring: bool,
    pub ungapped_extension: bool,
    pub x_dropoff_max: i32,
    pub cutoff_score_min: i32,
    pub container_type: ESeedContainerType,
}

use super::blast_program::{
    blast_program_is_nucleotide, blast_program_is_phi_blast, blast_query_is_pattern,
    EBlastProgramType,
};
use super::blast_query_info::BlastQueryInfo;
use super::blast_stat::{blast_cutoffs, BlastKarlinBlk, BlastScoreBlk};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[repr(i32)]
pub enum ESeedContainerType {
    #[default]
    DiagArray = 0,
    DiagHash = 1,
    MaxContainerType = 2,
}
#[derive(Clone, Copy, Debug, Default)]
pub struct BlastGappedCutoffs {
    pub cutoff_score: i32,
    pub cutoff_score_max: i32,
}
#[derive(Debug, Default)]
pub struct BlastLinkHSPParameters {
    pub gap_prob: f64,
    pub gap_size: i32,
    pub overlap_size: i32,
    pub cutoff_small_gap: i32,
    pub cutoff_big_gap: i32,
    pub longest_intron: i32,
    pub gap_decay_rate: f64,
}
#[derive(Debug, Default)]
pub struct BlastHitSavingParameters {
    pub options: Option<super::blast_options::BlastHitSavingOptions>,
    pub cutoff_score_min: i32,
    pub restricted_align: bool,
    pub do_sum_stats: bool,
    pub mask_level: i32,
    pub low_score: Option<Vec<i32>>,
    pub prelim_evalue: f64,
    pub cutoffs: Vec<BlastGappedCutoffs>,
    pub link_hsp_params: Option<BlastLinkHSPParameters>,
}

pub const BLASTERR_NOVALIDKARLINALTSCHUL: i16 = 104;
pub const CUTOFF_E_BLASTN: f64 = 0.05;
pub const CUTOFF_E_BLASTP: f64 = 1e-300;
pub const CUTOFF_E_BLASTX: f64 = 1.0;
pub const CUTOFF_E_TBLASTN: f64 = 1.0;
pub const CUTOFF_E_TBLASTX: f64 = 1e-300;

pub fn s_blast_karlin_blk_is_valid(kbp: Option<&BlastKarlinBlk>) -> bool {
    match kbp {
        None => false,
        Some(kbp) => kbp.lambda > 0.0 && kbp.k > 0.0 && kbp.h > 0.0,
    }
}
pub fn s_blast_find_valid_karlin_blk<'a>(
    kbp_in: &'a [Option<BlastKarlinBlk>],
    query_info: &BlastQueryInfo,
    kbp_ret: &mut Option<&'a BlastKarlinBlk>,
) -> i16 {
    let mut status = BLASTERR_NOVALIDKARLINALTSCHUL;
    for i in query_info.first_context..=query_info.last_context {
        let kbp = kbp_in[i as usize].as_ref();
        assert_eq!(
            s_blast_karlin_blk_is_valid(kbp),
            query_info.contexts[i as usize].is_valid
        );
        if s_blast_karlin_blk_is_valid(kbp) {
            *kbp_ret = kbp;
            status = 0;
            break;
        }
    }
    status
}
pub fn s_blast_find_smallest_lambda<'a>(
    kbp_in: &'a [Option<BlastKarlinBlk>],
    query_info: &BlastQueryInfo,
    mut kbp_out: Option<&mut Option<&'a BlastKarlinBlk>>,
) -> f64 {
    let mut min_lambda = i32::MAX as f64;
    for i in query_info.first_context..=query_info.last_context {
        let kbp = kbp_in[i as usize].as_ref();
        assert_eq!(
            s_blast_karlin_blk_is_valid(kbp),
            query_info.contexts[i as usize].is_valid
        );
        if s_blast_karlin_blk_is_valid(kbp) {
            let kbp = kbp.unwrap();
            if min_lambda > kbp.lambda {
                min_lambda = kbp.lambda;
                if let Some(out) = kbp_out.as_deref_mut() {
                    *out = Some(kbp);
                }
            }
        }
    }
    assert!(min_lambda > 0.0);
    min_lambda
}
pub fn s_get_cutoff_evalue(program: EBlastProgramType) -> f64 {
    match program {
        EBlastProgramType::Mapping | EBlastProgramType::Blastn => CUTOFF_E_BLASTN,
        EBlastProgramType::Blastp | EBlastProgramType::RpsBlast => CUTOFF_E_BLASTP,
        EBlastProgramType::Blastx => CUTOFF_E_BLASTX,
        EBlastProgramType::Tblastn
        | EBlastProgramType::PsiTblastn
        | EBlastProgramType::RpsTblastn => CUTOFF_E_TBLASTN,
        EBlastProgramType::Tblastx => CUTOFF_E_TBLASTX,
        _ => panic!("unsupported cutoff program"),
    }
}

pub fn blast_initial_word_parameters_update(
    program_number: EBlastProgramType,
    hit_params: &BlastHitSavingParameters,
    sbp: &BlastScoreBlk,
    query_info: &BlastQueryInfo,
    subj_length: u32,
    parameters: &mut BlastInitialWordParameters,
) -> i16 {
    let (kbp_array, gapped_calculation) = if let Some(gap) = &sbp.kbp_gap {
        (gap.as_slice(), true)
    } else if !sbp.kbp_std.is_empty() {
        (sbp.kbp_std.as_slice(), false)
    } else {
        return -1;
    };
    let gap_decay_rate = hit_params
        .link_hsp_params
        .as_ref()
        .map_or(0.0, |p| p.gap_decay_rate);
    let mut cutoff_min = i32::MAX;
    let mut xdrop_max = 0;
    for context in query_info.first_context..=query_info.last_context {
        let c = context as usize;
        let mut gap_trigger = i32::MAX;
        let mut new_cutoff = 1;
        let curr_cutoffs = &mut parameters.cutoffs[c];
        if !query_info.contexts[c].is_valid {
            curr_cutoffs.cutoff_score = i32::MAX;
            continue;
        }
        if !sbp.kbp_std.is_empty() {
            let kbp = sbp.kbp_std[c].as_ref();
            if s_blast_karlin_blk_is_valid(kbp) {
                let kbp = kbp.unwrap();
                gap_trigger = ((parameters.options.gap_trigger * std::f64::consts::LN_2
                    + kbp.log_k)
                    / kbp.lambda) as i32;
            }
        }
        if !gapped_calculation || sbp.matrix_only_scoring {
            let mut cutoff_e = s_get_cutoff_evalue(program_number);
            let mut query_length = query_info.contexts[c].query_length;
            assert!(query_length > 0);
            if matches!(
                program_number,
                EBlastProgramType::Blastn | EBlastProgramType::Mapping
            ) {
                query_length = query_length.wrapping_mul(2);
            }
            let kbp = kbp_array[c].as_ref().unwrap();
            assert!(s_blast_karlin_blk_is_valid(Some(kbp)));
            let searchsp = (subj_length as u64)
                .min(query_length as u64)
                .wrapping_mul(subj_length as u64) as i64;
            blast_cutoffs(
                &mut new_cutoff,
                &mut cutoff_e,
                kbp,
                searchsp,
                true,
                gap_decay_rate,
            );
            if program_number != EBlastProgramType::Blastn {
                new_cutoff = new_cutoff.min(gap_trigger);
            }
        } else {
            new_cutoff = gap_trigger;
        }
        new_cutoff = new_cutoff.wrapping_mul(sbp.scale_factor as i32);
        new_cutoff = new_cutoff.min(hit_params.cutoffs[c].cutoff_score_max);
        curr_cutoffs.cutoff_score = new_cutoff;
        curr_cutoffs.x_dropoff = if curr_cutoffs.x_dropoff_init == 0 {
            new_cutoff
        } else {
            curr_cutoffs.x_dropoff_init
        };
        if program_number == EBlastProgramType::Blastn && sbp.matrix_only_scoring {
            curr_cutoffs.x_dropoff = curr_cutoffs.x_dropoff_init;
        }
        if new_cutoff < cutoff_min {
            cutoff_min = new_cutoff;
        }
        if xdrop_max < curr_cutoffs.x_dropoff {
            xdrop_max = curr_cutoffs.x_dropoff;
        }
        if matches!(
            program_number,
            EBlastProgramType::Blastn | EBlastProgramType::Mapping
        ) {
            curr_cutoffs.reduced_nucl_cutoff_score = (0.8 * new_cutoff as f64) as i32;
        }
    }
    parameters.cutoff_score_min = cutoff_min;
    parameters.x_dropoff_max = xdrop_max;
    0
}
pub fn blast_initial_word_parameters_new(
    program_number: EBlastProgramType,
    word_options: &BlastInitialWordOptions,
    hit_params: &BlastHitSavingParameters,
    sbp: &BlastScoreBlk,
    query_info: &BlastQueryInfo,
    subject_length: u32,
    parameters: Option<&mut Option<BlastInitialWordParameters>>,
) -> i16 {
    let Some(parameters) = parameters else {
        return 0;
    };
    let mut kbp = None;
    let status = s_blast_find_valid_karlin_blk(sbp.ungapped_kbp(), query_info, &mut kbp);
    if status != 0 {
        return status;
    }
    let mut p = BlastInitialWordParameters {
        options: word_options.clone(),
        cutoffs: vec![BlastUngappedCutoffs::default(); (query_info.last_context + 1) as usize],
        nucl_score_table: [0; 256],
        matrix_only_scoring: false,
        ungapped_extension: !blast_program_is_phi_blast(program_number),
        x_dropoff_max: 0,
        cutoff_score_min: 0,
        container_type: ESeedContainerType::DiagArray,
    };
    for context in query_info.first_context..=query_info.last_context {
        let c = context as usize;
        if !query_info.contexts[c].is_valid {
            continue;
        }
        let kbp = sbp.ungapped_kbp()[c].as_ref().unwrap();
        assert!(s_blast_karlin_blk_is_valid(Some(kbp)));
        p.cutoffs[c].x_dropoff_init =
            if program_number == EBlastProgramType::Blastn && sbp.matrix_only_scoring {
                word_options.x_dropoff as i32
            } else {
                (sbp.scale_factor
                    * (word_options.x_dropoff * std::f64::consts::LN_2 / kbp.lambda).ceil())
                    as i32
            };
    }
    if blast_program_is_nucleotide(program_number)
        && !blast_query_is_pattern(program_number)
        && (query_info.contexts[query_info.last_context as usize].query_offset
            + query_info.contexts[query_info.last_context as usize].query_length)
            > 8000
    {
        p.container_type = ESeedContainerType::DiagHash;
    }
    let status = blast_initial_word_parameters_update(
        program_number,
        hit_params,
        sbp,
        query_info,
        subject_length,
        &mut p,
    );
    if matches!(
        program_number,
        EBlastProgramType::Blastn | EBlastProgramType::Mapping
    ) {
        for i in 0..256 {
            let mut score = 0i32;
            score = score.wrapping_add(if i & 3 != 0 { sbp.penalty } else { sbp.reward });
            score = score.wrapping_add(if (i >> 2) & 3 != 0 {
                sbp.penalty
            } else {
                sbp.reward
            });
            score = score.wrapping_add(if (i >> 4) & 3 != 0 {
                sbp.penalty
            } else {
                sbp.reward
            });
            score = score.wrapping_add(if i >> 6 != 0 { sbp.penalty } else { sbp.reward });
            p.nucl_score_table[i] = score;
        }
    }
    p.matrix_only_scoring = program_number == EBlastProgramType::Blastn && sbp.matrix_only_scoring;
    *parameters = Some(p);
    status
}

use super::blast_options::{BlastExtensionOptions, BlastScoringOptions};
#[derive(Debug)]
pub struct BlastExtensionParameters {
    pub options: BlastExtensionOptions,
    pub gap_x_dropoff: i32,
    pub gap_x_dropoff_final: i32,
}
#[derive(Debug)]
pub struct BlastScoringParameters {
    pub options: BlastScoringOptions,
    pub scale_factor: f64,
    pub reward: i32,
    pub penalty: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub shift_pen: i32,
}

pub fn blast_scoring_parameters_new(
    score_options: Option<&BlastScoringOptions>,
    sbp: &BlastScoreBlk,
    parameters: &mut Option<BlastScoringParameters>,
) -> i16 {
    let Some(options) = score_options else {
        return 1;
    };
    let scale_factor = sbp.scale_factor;
    *parameters = Some(BlastScoringParameters {
        options: options.clone(),
        scale_factor,
        reward: options.reward as i32,
        penalty: options.penalty as i32,
        gap_open: options.gap_open.wrapping_mul(scale_factor as i32),
        gap_extend: options.gap_extend.wrapping_mul(scale_factor as i32),
        shift_pen: options.shift_pen.wrapping_mul(scale_factor as i32),
    });
    0
}
pub fn blast_extension_parameters_new(
    program_number: EBlastProgramType,
    options: &BlastExtensionOptions,
    sbp: &BlastScoreBlk,
    query_info: &BlastQueryInfo,
    parameters: Option<&mut Option<BlastExtensionParameters>>,
) -> i16 {
    let Some(parameters) = parameters else {
        return 0;
    };
    if !sbp.ungapped_kbp().is_empty() {
        let mut kbp = None;
        let status = s_blast_find_valid_karlin_blk(sbp.ungapped_kbp(), query_info, &mut kbp);
        if status != 0 {
            return status;
        }
    } else {
        *parameters = None;
        return -1;
    }
    let mut params = BlastExtensionParameters {
        options: options.clone(),
        gap_x_dropoff: 0,
        gap_x_dropoff_final: 0,
    };
    if let Some(gap) = &sbp.kbp_gap {
        let min_lambda = s_blast_find_smallest_lambda(gap, query_info, None);
        params.gap_x_dropoff = (options.gap_x_dropoff * std::f64::consts::LN_2 / min_lambda) as i32;
        let drop_final = options.gap_x_dropoff_final * std::f64::consts::LN_2 / min_lambda;
        params.gap_x_dropoff_final = if drop_final > params.gap_x_dropoff as f64 {
            drop_final as i32
        } else {
            params.gap_x_dropoff
        };
    }
    if sbp.scale_factor > 1.0 {
        assert!(super::blast_program::blast_program_is_rps_blast(
            program_number
        ));
        params.gap_x_dropoff = params.gap_x_dropoff.wrapping_mul(sbp.scale_factor as i32);
        params.gap_x_dropoff_final = params
            .gap_x_dropoff_final
            .wrapping_mul(sbp.scale_factor as i32);
    }
    if program_number == EBlastProgramType::Blastn && sbp.matrix_only_scoring {
        params.gap_x_dropoff = options.gap_x_dropoff as i32;
        params.gap_x_dropoff_final = options.gap_x_dropoff_final as i32;
    }
    if program_number == EBlastProgramType::Mapping {
        params.gap_x_dropoff = options.gap_x_dropoff as i32;
    }
    *parameters = Some(params);
    0
}

use super::blast_options::BlastHitSavingOptions;
pub const BLAST_GAP_PROB: f64 = 0.5;
pub const BLAST_GAP_PROB_GAPPED: f64 = 1.0;
pub const BLAST_GAP_DECAY_RATE: f64 = 0.5;
pub const BLAST_GAP_DECAY_RATE_GAPPED: f64 = 0.1;
pub const BLAST_GAP_SIZE: i32 = 40;
pub const BLAST_OVERLAP_SIZE: i32 = 9;
pub const DEFAULT_LONGEST_INTRON: i32 = 122;
pub const RESTRICTED_ALIGNMENT_WORST_EVALUE: f64 = 10.0;

pub fn blast_link_hsp_parameters_new(
    program_number: EBlastProgramType,
    gapped_calculation: bool,
    link_hsp_params: Option<&mut Option<BlastLinkHSPParameters>>,
) -> i16 {
    let Some(out) = link_hsp_params else {
        return -1;
    };
    let mut params = BlastLinkHSPParameters::default();
    if program_number == EBlastProgramType::Blastn || !gapped_calculation {
        params.gap_prob = BLAST_GAP_PROB;
        params.gap_decay_rate = BLAST_GAP_DECAY_RATE;
    } else {
        params.gap_prob = BLAST_GAP_PROB_GAPPED;
        params.gap_decay_rate = BLAST_GAP_DECAY_RATE_GAPPED;
    }
    params.gap_size = BLAST_GAP_SIZE;
    params.overlap_size = BLAST_OVERLAP_SIZE;
    *out = Some(params);
    0
}
pub fn blast_link_hsp_parameters_update(
    word_params: Option<&BlastInitialWordParameters>,
    hit_params: Option<&mut BlastHitSavingParameters>,
    _gapped_calculation: bool,
) -> i16 {
    let (Some(word), Some(hit)) = (word_params, hit_params) else {
        return -1;
    };
    let Some(link) = &mut hit.link_hsp_params else {
        return 0;
    };
    link.cutoff_small_gap = word.cutoff_score_min;
    0
}

/// Default nucleotide and tblastn cutoff paths, including Gumbel/Spouge.
/// PHI and RPS-only branches remain outside the supported programs.
pub fn blast_hit_saving_parameters_update(
    program_number: EBlastProgramType,
    sbp: &BlastScoreBlk,
    query_info: &BlastQueryInfo,
    avg_subject_length: i32,
    composition_based_stats: i32,
    params: &mut BlastHitSavingParameters,
) -> i16 {
    assert!(
        matches!(
            program_number,
            EBlastProgramType::Blastn | EBlastProgramType::Mapping | EBlastProgramType::Tblastn
        ),
        "PHI/RPS hit-saving setup is outside the supported default programs"
    );
    let options = params
        .options
        .as_ref()
        .expect("hit-saving options required");
    let (kbp_array, gapped_calculation) = if let Some(gap) = &sbp.kbp_gap {
        (gap.as_slice(), true)
    } else if !sbp.ungapped_kbp().is_empty() {
        (sbp.ungapped_kbp(), false)
    } else {
        return -1;
    };
    params.prelim_evalue = options.expect_value;
    if matches!(
        program_number,
        EBlastProgramType::Blastn | EBlastProgramType::Mapping
    ) && options.mask_level >= 0
    {
        params.mask_level = options.mask_level;
    }
    if options.cutoff_score > 0 {
        let new_cutoff = options.cutoff_score.wrapping_mul(sbp.scale_factor as i32);
        for context in query_info.first_context..=query_info.last_context {
            let cut = &mut params.cutoffs[context as usize];
            cut.cutoff_score = new_cutoff;
            cut.cutoff_score_max = new_cutoff;
            if program_number == EBlastProgramType::Blastn && sbp.matrix_only_scoring {
                cut.cutoff_score = options.cutoff_score;
                cut.cutoff_score_max = options.cutoff_score / 2;
            }
        }
        params.cutoff_score_min = new_cutoff;
    } else {
        let mut cutoff_min = i32::MAX;
        for context in query_info.first_context..=query_info.last_context {
            let c = context as usize;
            let mut new_cutoff = 1;
            let mut evalue = options.expect_value;
            if !query_info.contexts[c].is_valid {
                params.cutoffs[c].cutoff_score = i32::MAX;
                continue;
            }
            let kbp = kbp_array[c].as_ref().unwrap();
            assert!(s_blast_karlin_blk_is_valid(Some(kbp)));
            params.prelim_evalue = evalue;
            if let Some(gbp) = sbp.gbp.as_ref().filter(|g| g.filled) {
                let stretch = if composition_based_stats > 1 { 5 } else { 1 };
                params.prelim_evalue = stretch as f64 * evalue;
                new_cutoff = super::blast_stat_spouge::blast_spouge_eto_s(
                    stretch as f64 * evalue,
                    kbp,
                    gbp,
                    query_info.contexts[c].query_length,
                    avg_subject_length,
                );
            } else {
                blast_cutoffs(
                    &mut new_cutoff,
                    &mut evalue,
                    kbp,
                    query_info.contexts[c].eff_searchsp,
                    false,
                    0.0,
                );
            }
            params.cutoffs[c].cutoff_score = new_cutoff;
            params.cutoffs[c].cutoff_score_max = new_cutoff;
        }
        if let Some(link) = &params.link_hsp_params {
            if gapped_calculation {
                let mut evalue_hsp = 1.0;
                let last = &query_info.contexts[query_info.last_context as usize];
                let avg_qlen =
                    (last.query_offset + last.query_length) / (query_info.last_context + 1);
                let searchsp = (avg_qlen.min(avg_subject_length) as i64)
                    .wrapping_mul(avg_subject_length as i64);
                for context in query_info.first_context..=query_info.last_context {
                    let c = context as usize;
                    if !query_info.contexts[c].is_valid {
                        continue;
                    }
                    let mut new_cutoff = 1;
                    let kbp = kbp_array[c].as_ref().unwrap();
                    assert!(s_blast_karlin_blk_is_valid(Some(kbp)));
                    blast_cutoffs(
                        &mut new_cutoff,
                        &mut evalue_hsp,
                        kbp,
                        searchsp,
                        true,
                        link.gap_decay_rate,
                    );
                    params.cutoffs[c].cutoff_score = new_cutoff.min(params.cutoffs[c].cutoff_score);
                }
            }
        }
        for context in query_info.first_context..=query_info.last_context {
            let c = context as usize;
            if query_info.contexts[c].is_valid {
                params.cutoffs[c].cutoff_score = params.cutoffs[c]
                    .cutoff_score
                    .wrapping_mul(sbp.scale_factor as i32);
                params.cutoffs[c].cutoff_score_max = params.cutoffs[c]
                    .cutoff_score_max
                    .wrapping_mul(sbp.scale_factor as i32);
                cutoff_min = cutoff_min.min(params.cutoffs[c].cutoff_score);
            }
        }
        params.cutoff_score_min = cutoff_min;
    }
    0
}
pub fn blast_hit_saving_parameters_new(
    program_number: EBlastProgramType,
    options: &BlastHitSavingOptions,
    sbp: &BlastScoreBlk,
    query_info: &BlastQueryInfo,
    avg_subj_length: i32,
    composition_based_stats: i32,
    parameters: Option<&mut Option<BlastHitSavingParameters>>,
) -> i16 {
    let Some(parameters) = parameters else {
        return 0;
    };
    *parameters = None;
    let gapped_calculation = sbp.kbp_gap.is_some();
    if options.do_sum_stats && gapped_calculation && avg_subj_length <= 0 {
        return 1;
    }
    let mut params = BlastHitSavingParameters {
        options: Some(options.clone()),
        mask_level: 101,
        do_sum_stats: options.do_sum_stats,
        cutoffs: vec![BlastGappedCutoffs::default(); (query_info.last_context + 1) as usize],
        ..Default::default()
    };
    if params.do_sum_stats {
        blast_link_hsp_parameters_new(
            program_number,
            gapped_calculation,
            Some(&mut params.link_hsp_params),
        );
        if (super::blast_program::blast_query_is_translated(program_number)
            || super::blast_program::blast_subject_is_translated(program_number))
            && program_number != EBlastProgramType::Tblastx
        {
            let max_protein_gap = (options.longest_intron - 2) / 3;
            if gapped_calculation {
                if options.longest_intron == 0 {
                    params.link_hsp_params.as_mut().unwrap().longest_intron =
                        (DEFAULT_LONGEST_INTRON - 2) / 3;
                } else if max_protein_gap <= 0 {
                    params.link_hsp_params = None;
                    params.do_sum_stats = false;
                } else {
                    params.link_hsp_params.as_mut().unwrap().longest_intron = max_protein_gap;
                }
            } else {
                params.link_hsp_params.as_mut().unwrap().longest_intron = max_protein_gap.max(0);
            }
        }
    }
    if options.low_score_perc > 0.00001 {
        params.low_score = Some(vec![0; query_info.num_queries as usize]);
    }
    let status = blast_hit_saving_parameters_update(
        program_number,
        sbp,
        query_info,
        avg_subj_length,
        composition_based_stats,
        &mut params,
    );
    *parameters = Some(params);
    status
}
