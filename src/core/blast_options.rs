// NCBI BLAST+ 2.17.0+, selected blast_options.c/h definitions.
// Authors: NCBI BLAST developers. See LICENSE-NCBI.

/// Fields read by the translated mask-at-hash predicate. Other filtering
/// options will be added when their algorithms require them.
#[derive(Debug, Default)]
pub struct SBlastFilterOptions {
    pub mask_at_hash: bool,
}

/// Fields read by `s_HasMaskAtHashEnabled`. The deprecated filter string is
/// stored as bytes; its first NUL terminates it, as with the original C string.
#[derive(Debug, Default)]
pub struct QuerySetUpOptions {
    pub filtering_options: Option<SBlastFilterOptions>,
    pub filter_string: Option<Vec<u8>>,
}

/// Matches C `SBlastFilterOptionsMaskAtHash`.
pub fn s_blast_filter_options_mask_at_hash(filter_options: Option<&SBlastFilterOptions>) -> bool {
    let Some(filter_options) = filter_options else {
        return false;
    };
    filter_options.mask_at_hash
}

/// Original initial-word options. Values are supplied explicitly until the
/// upstream defaults and statistical parameter construction are translated.
#[derive(Clone, Debug)]
pub struct BlastInitialWordOptions {
    pub gap_trigger: f64,
    pub window_size: i32,
    pub scan_range: i32,
    pub x_dropoff: f64,
    pub program_number: super::blast_program::EBlastProgramType,
}

pub const BLAST_WINDOW_SIZE_PROT: i32 = 40;
pub const BLAST_WINDOW_SIZE_NUCL: i32 = 0;
pub const BLAST_SCAN_RANGE_NUCL: i32 = 0;
pub const BLAST_UNGAPPED_X_DROPOFF_PROT: i32 = 7;
pub const BLAST_UNGAPPED_X_DROPOFF_NUCL: i32 = 20;
pub const BLAST_GAP_TRIGGER_PROT: f64 = 22.0;
pub const BLAST_GAP_TRIGGER_NUCL: f64 = 27.0;

pub fn blast_initial_word_options_new(
    program: super::blast_program::EBlastProgramType,
) -> BlastInitialWordOptions {
    let mut options = BlastInitialWordOptions {
        gap_trigger: 0.0,
        window_size: 0,
        scan_range: 0,
        x_dropoff: 0.0,
        program_number: program,
    };
    if !super::blast_program::blast_program_is_nucleotide(program) {
        options.window_size = BLAST_WINDOW_SIZE_PROT;
        options.x_dropoff = BLAST_UNGAPPED_X_DROPOFF_PROT as f64;
        options.gap_trigger = BLAST_GAP_TRIGGER_PROT;
    } else {
        options.window_size = BLAST_WINDOW_SIZE_NUCL;
        options.scan_range = BLAST_SCAN_RANGE_NUCL;
        options.gap_trigger = BLAST_GAP_TRIGGER_NUCL;
        options.x_dropoff = BLAST_UNGAPPED_X_DROPOFF_NUCL as f64;
    }
    options.program_number = program;
    options
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(i32)]
pub enum EBlastPrelimGapExt {
    DynProgScoreOnly = 0,
    GreedyScoreOnly = 1,
    JumperWithTraceback = 2,
    SmithWatermanScoreOnly = 3,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(i32)]
pub enum EBlastTbackExt {
    DynProgTbck = 0,
    GreedyTbck = 1,
    SmithWatermanTbck = 2,
    SmithWatermanTbckFull = 3,
}
#[derive(Clone, Debug)]
pub struct BlastExtensionOptions {
    pub gap_x_dropoff: f64,
    pub gap_x_dropoff_final: f64,
    pub e_prelim_gap_ext: EBlastPrelimGapExt,
    pub e_tback_ext: EBlastTbackExt,
    pub composition_based_stats: i32,
    pub unified_p: i32,
    pub max_mismatches: i32,
    pub mismatch_window: i32,
    pub chaining: bool,
    pub program_number: super::blast_program::EBlastProgramType,
}
#[derive(Clone, Debug)]
pub struct BlastScoringOptions {
    pub matrix: Option<String>,
    pub matrix_path: Option<String>,
    pub reward: i16,
    pub penalty: i16,
    pub gapped_calculation: bool,
    pub complexity_adjusted_scoring: bool,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub is_ooframe: bool,
    pub shift_pen: i32,
    pub program_number: super::blast_program::EBlastProgramType,
}
pub const BLAST_GAP_X_DROPOFF_PROT: i32 = 15;
pub const BLAST_GAP_X_DROPOFF_NUCL: i32 = 30;
pub const BLAST_GAP_X_DROPOFF_GREEDY: i32 = 25;
pub const BLAST_GAP_X_DROPOFF_FINAL_PROT: i32 = 25;
pub const BLAST_GAP_X_DROPOFF_FINAL_NUCL: i32 = 100;
pub const BLAST_GAP_OPEN_PROT: i32 = 11;
pub const BLAST_GAP_OPEN_NUCL: i32 = 5;
pub const BLAST_GAP_EXTN_PROT: i32 = 1;
pub const BLAST_GAP_EXTN_NUCL: i32 = 2;
pub const BLAST_PENALTY: i16 = -3;
pub const BLAST_REWARD: i16 = 1;
pub const BLAST_DEFAULT_MATRIX: &str = "BLOSUM62";

pub fn blast_extension_options_new(
    program: super::blast_program::EBlastProgramType,
    gapped: bool,
) -> BlastExtensionOptions {
    use super::blast_program::*;
    let mut options = BlastExtensionOptions {
        gap_x_dropoff: 0.0,
        gap_x_dropoff_final: 0.0,
        e_prelim_gap_ext: EBlastPrelimGapExt::DynProgScoreOnly,
        e_tback_ext: EBlastTbackExt::DynProgTbck,
        composition_based_stats: 0,
        unified_p: 0,
        max_mismatches: 0,
        mismatch_window: 0,
        chaining: false,
        program_number: program,
    };
    if !blast_program_is_nucleotide(program) {
        options.gap_x_dropoff = BLAST_GAP_X_DROPOFF_PROT as f64;
        options.gap_x_dropoff_final = BLAST_GAP_X_DROPOFF_FINAL_PROT as f64;
    } else {
        options.gap_x_dropoff = BLAST_GAP_X_DROPOFF_NUCL as f64;
        options.gap_x_dropoff_final = BLAST_GAP_X_DROPOFF_FINAL_NUCL as f64;
    }
    if gapped && blast_query_is_pssm(program) && !blast_subject_is_translated(program) {
        options.composition_based_stats = 1;
    }
    options.max_mismatches = 5;
    options.mismatch_window = 10;
    options.program_number = program;
    options
}
pub fn blast_scoring_options_new(
    program_number: super::blast_program::EBlastProgramType,
) -> BlastScoringOptions {
    use super::blast_program::*;
    let mut options = BlastScoringOptions {
        matrix: None,
        matrix_path: None,
        reward: 0,
        penalty: 0,
        gapped_calculation: false,
        complexity_adjusted_scoring: false,
        gap_open: 0,
        gap_extend: 0,
        is_ooframe: false,
        shift_pen: 0,
        program_number,
    };
    if !blast_program_is_nucleotide(program_number) {
        options.shift_pen = i16::MAX as i32;
        options.is_ooframe = false;
        options.gap_open = BLAST_GAP_OPEN_PROT;
        options.gap_extend = BLAST_GAP_EXTN_PROT;
        options.matrix = Some(BLAST_DEFAULT_MATRIX.to_owned());
    } else {
        options.penalty = BLAST_PENALTY;
        options.reward = BLAST_REWARD;
        options.gap_open = BLAST_GAP_OPEN_NUCL;
        options.gap_extend = BLAST_GAP_EXTN_NUCL;
    }
    if program_number != EBlastProgramType::Tblastx {
        options.gapped_calculation = true;
    }
    options.program_number = program_number;
    options.complexity_adjusted_scoring = false;
    options
}

#[derive(Clone, Debug)]
pub struct BlastHitSavingOptions {
    pub expect_value: f64,
    pub cutoff_score: i32,
    pub cutoff_score_fun: [i32; 2],
    pub percent_identity: f64,
    pub max_edit_distance: i32,
    pub hitlist_size: i32,
    pub hsp_num_max: i32,
    pub total_hsp_limit: i32,
    pub culling_limit: i32,
    pub mask_level: i32,
    pub do_sum_stats: bool,
    pub longest_intron: i32,
    pub min_hit_length: i32,
    pub min_diag_separation: i32,
    pub program_number: super::blast_program::EBlastProgramType,
    pub low_score_perc: f64,
    pub query_cov_hsp_perc: f64,
    pub max_hsps_per_subject: i32,
    pub paired: bool,
    pub splice: bool,
}
pub const BLAST_HITLIST_SIZE: i32 = 500;
pub const BLAST_EXPECT_VALUE: f64 = 10.0;
pub fn blast_hit_saving_options_new(
    program_number: super::blast_program::EBlastProgramType,
    gapped_calculation: bool,
) -> BlastHitSavingOptions {
    use super::blast_program::*;
    BlastHitSavingOptions {
        expect_value: BLAST_EXPECT_VALUE,
        cutoff_score: 0,
        cutoff_score_fun: [0; 2],
        percent_identity: 0.0,
        max_edit_distance: i32::MAX,
        hitlist_size: BLAST_HITLIST_SIZE,
        hsp_num_max: 0,
        total_hsp_limit: 0,
        culling_limit: 0,
        mask_level: 101,
        do_sum_stats: program_number != EBlastProgramType::RpsTblastn
            && (!gapped_calculation
                || blast_query_is_translated(program_number)
                || blast_subject_is_translated(program_number)),
        longest_intron: 0,
        min_hit_length: 0,
        min_diag_separation: 0,
        program_number,
        low_score_perc: 0.0,
        query_cov_hsp_perc: 0.0,
        max_hsps_per_subject: 0,
        paired: false,
        splice: false,
    }
}
