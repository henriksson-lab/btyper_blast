//! Standard nucleotide setup, preliminary extension and per-subject traceback.
//! Database-wide candidate/result selection is assembled separately.
use crate::core::{
    blast_diagnostics::*,
    blast_dynamic_align::BlastDynamicAlign,
    blast_effective_lengths::*,
    blast_extend::*,
    blast_gapalign::{blast_gap_align_struct_new, BlastGapAlignStruct},
    blast_lookup::SSeqRange,
    blast_nascan::*,
    blast_options::*,
    blast_parameters::*,
    blast_program::EBlastProgramType,
    blast_query_info::BlastQueryInfo,
    blast_stat::*,
    blast_stat_gapped::*,
    na_ungapped::*,
};
use crate::seed::{NucleotideSeedIndex, SeedError};
use crate::{NucleotideDatabase, NucleotideSequence, NucleotideTask, PreparedNucleotideBatch};
use std::fmt;

#[derive(Debug)]
pub enum NucleotideSetupError {
    Seed(SeedError),
    InvalidSubjectRange { offset: i32, length: i32 },
    Upstream { stage: &'static str, status: i16 },
}
impl fmt::Display for NucleotideSetupError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Seed(e) => e.fmt(f),
            Self::InvalidSubjectRange { offset, length } => {
                write!(
                    f,
                    "invalid packed subject range: offset {offset}, length {length}"
                )
            }
            Self::Upstream { stage, status } => {
                write!(f, "BLAST {stage} failed with status {status}")
            }
        }
    }
}
impl std::error::Error for NucleotideSetupError {}
impl From<SeedError> for NucleotideSetupError {
    fn from(e: SeedError) -> Self {
        Self::Seed(e)
    }
}

/// Reusable standard nucleotide word-finder state, with original score,
/// effective-length and cutoff calculations. Coordinates retain concatenated
/// context offsets. The returned hits still need gapped extension and final
/// processing before they are BLAST search results.
pub struct NucleotideInitialSearch<'a> {
    batch: &'a PreparedNucleotideBatch,
    index: Option<NucleotideSeedIndex>,
    query_info: BlastQueryInfo,
    pub score_block: BlastScoreBlk,
    pub scoring_options: BlastScoringOptions,
    pub extension_parameters: BlastExtensionParameters,
    pub word_parameters: BlastInitialWordParameters,
    pub hit_parameters: BlastHitSavingParameters,
    pub round_down: bool,
    pub warnings: Vec<StatisticalWarning>,
    matrix: NucleotideScoreMatrix,
    // Independent constant rows permit statistics to mutably borrow ScoreBlk.
    traceback_matrix: Vec<Vec<i32>>,
    state: BlastExtendWord,
    pairs: Vec<BlastOffsetPair>,
    initial_hits: BlastInitHitList,
    stats: BlastUngappedStats,
    // Source GapAlignSetUp owns one scratch allocation throughout the search.
    // The stored view is empty; methods bind a short matrix borrow and return
    // the same owned buffers here without constructing a self-reference.
    gap_align: Option<BlastGapAlignStruct<'static>>,
}
impl<'a> NucleotideInitialSearch<'a> {
    pub fn new(
        batch: &'a PreparedNucleotideBatch,
        database: &NucleotideDatabase,
    ) -> Result<Self, NucleotideSetupError> {
        Self::new_with_effective_searchspaces(batch, database, &[])
    }
    /// Split-query setup supplies the full query's effective spaces before
    /// deriving chunk hit and word thresholds, as in the original splitter.
    pub(crate) fn new_with_effective_searchspaces(
        batch: &'a PreparedNucleotideBatch,
        database: &NucleotideDatabase,
        effective_searchspaces: &[i64],
    ) -> Result<Self, NucleotideSetupError> {
        Self::new_internal(batch, database, effective_searchspaces, true)
    }
    /// Original full-query traceback setup does not construct a full lookup
    /// when the preliminary search uses split query factories.
    pub(crate) fn new_for_traceback(
        batch: &'a PreparedNucleotideBatch,
        database: &NucleotideDatabase,
    ) -> Result<Self, NucleotideSetupError> {
        Self::new_internal(batch, database, &[], false)
    }
    fn new_internal(
        batch: &'a PreparedNucleotideBatch,
        database: &NucleotideDatabase,
        effective_searchspaces: &[i64],
        build_lookup: bool,
    ) -> Result<Self, NucleotideSetupError> {
        let program = EBlastProgramType::Blastn;
        let mut scoring_options = blast_scoring_options_new(program);
        // Task overrides from CBlastNucleotideOptionsHandle's scoring defaults.
        match batch.task() {
            NucleotideTask::Megablast => {
                scoring_options.reward = 1;
                scoring_options.penalty = -2;
                scoring_options.gap_open = 0;
                scoring_options.gap_extend = 0;
            }
            NucleotideTask::Blastn => {
                scoring_options.reward = 2;
                scoring_options.penalty = -3;
                scoring_options.gap_open = 5;
                scoring_options.gap_extend = 2;
            }
        }
        scoring_options.shift_pen = i16::MAX as i32;
        let mut matrix = [[0; 16]; 16];
        blast_score_blk_nucl_matrix_create(
            scoring_options.reward as i32,
            scoring_options.penalty as i32,
            &mut matrix,
        );
        let mut score_block = BlastScoreBlk {
            alphabet_code: crate::core::blast_encoding::BLASTNA_SEQ_CODE as u8,
            alphabet_size: 16,
            ambiguous_res: vec![14],
            reward: scoring_options.reward as i32,
            penalty: scoring_options.penalty as i32,
            scale_factor: 1.0,
            matrix: matrix.iter().map(|r| r.to_vec()).collect(),
            ..Default::default()
        };
        blast_score_blk_max_score_set(&mut score_block);
        let mut query_info = batch.query_info().clone();
        let mut warnings = vec![];
        let status = blast_score_blk_kbp_ungapped_calc(
            program,
            &mut score_block,
            batch.blastna_contexts(),
            &mut query_info,
            Some(&mut warnings),
        );
        if status != 0 {
            return Err(NucleotideSetupError::Upstream {
                stage: "ungapped statistics",
                status,
            });
        }
        let mut gapped = vec![None; query_info.contexts.len()];
        let mut round_down = false;
        for c in query_info.first_context..=query_info.last_context {
            if !query_info.contexts[c as usize].is_valid {
                continue;
            }
            let mut kbp = blast_karlin_blk_new();
            let status = blast_karlin_blk_nucl_gapped_calc(
                &mut kbp,
                scoring_options.gap_open,
                scoring_options.gap_extend,
                score_block.reward,
                score_block.penalty,
                score_block.kbp_std[c as usize].as_ref().unwrap(),
                &mut round_down,
                None,
            );
            if status != 0 {
                return Err(NucleotideSetupError::Upstream {
                    stage: "gapped statistics",
                    status,
                });
            }
            gapped[c as usize] = Some(kbp);
        }
        score_block.kbp_gap = Some(gapped);
        score_block.round_down = round_down;
        let mut length_options = blast_effective_lengths_options_new();
        if !effective_searchspaces.is_empty() {
            // Original chunk setup indexes this full vector by LOCAL context.
            length_options.num_searchspaces = effective_searchspaces.len() as i32;
            length_options.searchsp_eff = Some(effective_searchspaces.to_vec());
        }
        let length_parameters = blast_effective_lengths_parameters_new(
            &length_options,
            database.total_length() as i64,
            database.len() as i32,
        );
        let status = blast_calc_eff_lengths(
            program,
            &EffectiveLengthScoringOptions {
                gapped_calculation: true,
                reward: score_block.reward,
                penalty: score_block.penalty,
                gap_open: scoring_options.gap_open,
                gap_extend: scoring_options.gap_extend,
            },
            &length_parameters,
            Some(&score_block),
            Some(&mut query_info),
            None,
            None,
        );
        if status != 0 {
            return Err(NucleotideSetupError::Upstream {
                stage: "effective lengths",
                status,
            });
        }
        let avg_length = (database.total_length() / database.len() as u64) as i32;
        let mut hit_options = blast_hit_saving_options_new(program, true);
        // CBlastNucleotideOptionsHandle's task presets override the core
        // constructor's zero: SetHitSavingOptionsDefaults uses50, whereas
        // SetMBHitSavingOptionsDefaults uses6. Containment must retain hits
        // separated by more diagonals, including alternate traceback seeds.
        hit_options.min_diag_separation = match batch.task() {
            NucleotideTask::Blastn => 50,
            NucleotideTask::Megablast => 6,
        };
        let mut hit_parameters = None;
        let status = blast_hit_saving_parameters_new(
            program,
            &hit_options,
            &score_block,
            &query_info,
            avg_length,
            0,
            Some(&mut hit_parameters),
        );
        if status != 0 {
            return Err(NucleotideSetupError::Upstream {
                stage: "hit-saving parameters",
                status,
            });
        }
        let hit_parameters = hit_parameters.unwrap();
        let mut word_parameters = None;
        let status = blast_initial_word_parameters_new(
            program,
            &blast_initial_word_options_new(program),
            &hit_parameters,
            &score_block,
            &query_info,
            avg_length as u32,
            Some(&mut word_parameters),
        );
        if status != 0 {
            return Err(NucleotideSetupError::Upstream {
                stage: "initial-word parameters",
                status,
            });
        }
        let word_parameters = word_parameters.unwrap();
        let mut extension_options = blast_extension_options_new(program, true);
        if batch.task() == NucleotideTask::Megablast {
            extension_options.gap_x_dropoff = BLAST_GAP_X_DROPOFF_GREEDY as f64;
            extension_options.e_prelim_gap_ext = EBlastPrelimGapExt::GreedyScoreOnly;
            extension_options.e_tback_ext = EBlastTbackExt::GreedyTbck;
        }
        let mut extension_parameters = None;
        let status = blast_extension_parameters_new(
            program,
            &extension_options,
            &score_block,
            &query_info,
            Some(&mut extension_parameters),
        );
        if status != 0 {
            return Err(NucleotideSetupError::Upstream {
                stage: "extension parameters",
                status,
            });
        }
        let (index, state, pairs) = if build_lookup {
            let index = batch.seed_index()?;
            let state = blast_extend_word_new(
                batch.blastna_contexts().len() as u32,
                word_parameters.container_type == ESeedContainerType::DiagHash,
                word_parameters.options.window_size,
            );
            let pairs = vec![BlastOffsetPair::default(); index.recommended_capacity()];
            (Some(index), state, pairs)
        } else {
            (None, BlastExtendWord::default(), Vec::new())
        };
        let extension_parameters = extension_parameters.unwrap();
        let mut scoring = None;
        let status =
            blast_scoring_parameters_new(Some(&scoring_options), &score_block, &mut scoring);
        if status != 0 {
            return Err(NucleotideSetupError::Upstream {
                stage: "scoring parameters",
                status,
            });
        }
        let mut gap_align = None;
        let status = blast_gap_align_struct_new(
            scoring.as_ref(),
            Some(&extension_parameters),
            database.max_length() as u32,
            Some(&score_block),
            Some(&mut gap_align),
        );
        if status != 0 {
            return Err(NucleotideSetupError::Upstream {
                stage: "gapped scratch",
                status,
            });
        }
        let gap_align = Some(bind_gap_align(gap_align.unwrap(), &[]));
        let traceback_matrix = score_block.matrix.clone();
        Ok(Self {
            batch,
            index,
            query_info,
            score_block,
            scoring_options,
            extension_parameters,
            word_parameters,
            hit_parameters,
            round_down,
            warnings,
            matrix,
            traceback_matrix,
            state,
            pairs,
            initial_hits: blast_init_hit_list_new(),
            stats: BlastUngappedStats::default(),
            gap_align,
        })
    }
    pub fn query_info(&self) -> &BlastQueryInfo {
        &self.query_info
    }
    pub fn stats(&self) -> &BlastUngappedStats {
        &self.stats
    }
    pub fn scan_subject(
        &mut self,
        subject: &NucleotideSequence,
    ) -> Result<&[BlastInitHSP], NucleotideSetupError> {
        self.scan_subject_chunk(subject, 0, subject.len() as i32)
    }
    pub(crate) fn scan_subject_chunk(
        &mut self,
        subject: &NucleotideSequence,
        offset: i32,
        length: i32,
    ) -> Result<&[BlastInitHSP], NucleotideSetupError> {
        if offset < 0
            || offset % 4 != 0
            || length < 0
            || offset as usize + length as usize > subject.len()
        {
            return Err(NucleotideSetupError::InvalidSubjectRange { offset, length });
        }
        blast_init_hit_list_reset(&mut self.initial_hits);
        let index = self.index.as_ref().expect("preliminary lookup owner");
        let mut lookup = index.lookup_wrap();
        let extend = blast_choose_na_extend(&mut lookup);
        let query = NucleotideExtensionQuery {
            sequence_start: self.batch.sequence_with_sentinels(),
            compressed_nuc_seq_start: index.compressed_query(),
        };
        let capacity = self.pairs.len() as i32;
        let status = blast_na_word_finder(
            &PackedNucleotideSubject {
                sequence: &subject.packed_for_scanning()[offset as usize / 4..],
                mask_type: ESubjectMaskingType::NoSubjMasking,
            },
            length,
            &[SSeqRange {
                left: 0,
                right: length,
            }],
            &query,
            &self.query_info,
            &lookup,
            &self.matrix,
            &self.word_parameters,
            &mut self.state,
            &mut self.pairs,
            capacity,
            &mut self.initial_hits,
            Some(&mut self.stats),
            index.scan_callback(),
            extend,
        );
        if status != 0 {
            return Err(NucleotideSetupError::Upstream {
                stage: "word finder",
                status,
            });
        }
        Ok(&self.initial_hits.init_hsp_array)
    }
    /// Original preliminary gapped extension/containment over this subject.
    /// Traceback, ambiguity reevaluation and final HSP filtering are subsequent stages.
    pub fn preliminary_subject(
        &mut self,
        subject: &NucleotideSequence,
    ) -> Result<Box<crate::core::blast_hits::BlastHSPList>, NucleotideSetupError> {
        self.preliminary_subject_chunk(subject, 0, subject.len() as i32)
    }
    pub(crate) fn preliminary_subject_chunk(
        &mut self,
        subject: &NucleotideSequence,
        offset: i32,
        length: i32,
    ) -> Result<Box<crate::core::blast_hits::BlastHSPList>, NucleotideSetupError> {
        use crate::core::{
            blast_dynamic_align::DynamicSequence,
            blast_hits::blast_hsp_list_new,
            blast_prelim::{blast_get_gapped_score, GapSequence, GapSubject},
        };
        self.scan_subject_chunk(subject, offset, length)?;
        let mut scoring = None;
        let status = blast_scoring_parameters_new(
            Some(&self.scoring_options),
            &self.score_block,
            &mut scoring,
        );
        if status != 0 {
            return Err(NucleotideSetupError::Upstream {
                stage: "scoring parameters",
                status,
            });
        }
        let scoring = scoring.unwrap();
        let mut gap = bind_gap_align(
            self.gap_align.take().expect("search-owned gapped scratch"),
            &self.score_block.matrix,
        );
        gap.dynamic.gap_x_dropoff = self.extension_parameters.gap_x_dropoff;
        let mut list = None;
        let query = GapSequence {
            sequence: Some(DynamicSequence {
                data: self.batch.sequence_with_sentinels(),
                offset: 1,
            }),
            oof_sequence: None,
            length: self.batch.blastna_contexts().len() as i32,
        };
        let target = GapSubject {
            sequence: DynamicSequence {
                data: &subject.packed_for_scanning()[offset as usize / 4..],
                offset: 0,
            },
            length,
            frame: 1,
        };
        let status = blast_get_gapped_score(
            EBlastProgramType::Blastn,
            &query,
            &self.query_info,
            &target,
            &mut gap,
            &scoring,
            &self.extension_parameters,
            &self.hit_parameters,
            &self.word_parameters,
            &self.initial_hits,
            Some(&mut list),
            None,
            None,
        );
        self.gap_align = Some(bind_gap_align(gap, &[]));
        if status != 0 {
            return Err(NucleotideSetupError::Upstream {
                stage: "preliminary gapped alignment",
                status,
            });
        }
        Ok(list.unwrap_or_else(|| blast_hsp_list_new(i32::MAX)))
    }

    /// Extend one subject through traceback, ambiguity reevaluation and final
    /// per-HSP statistics. Coordinates are zero-based, context-local and half-open.
    /// Database-wide candidate selection and result limits are a separate stage.
    pub fn traceback_subject(
        &mut self,
        subject: &NucleotideSequence,
    ) -> Result<Box<crate::core::blast_hits::BlastHSPList>, NucleotideSetupError> {
        let list = self.preliminary_subject(subject)?;
        self.traceback_preliminary_subject(subject, list)
    }

    /// Trace back the actual HSP list retained by preliminary candidate selection.
    pub fn traceback_preliminary_subject(
        &mut self,
        subject: &NucleotideSequence,
        list: Box<crate::core::blast_hits::BlastHSPList>,
    ) -> Result<Box<crate::core::blast_hits::BlastHSPList>, NucleotideSetupError> {
        let mut expanded = Vec::new();
        expand_nucleotide_subject(subject, &mut expanded);
        self.traceback_preliminary_subject_expanded(&expanded, list)
    }

    /// Source OID-batch binding: caller retains the restored subject across all
    /// query lists in the batch. Buffer includes the two original sentinels.
    pub(crate) fn traceback_preliminary_subject_expanded(
        &mut self,
        expanded: &[u8],
        mut list: Box<crate::core::blast_hits::BlastHSPList>,
    ) -> Result<Box<crate::core::blast_hits::BlastHSPList>, NucleotideSetupError> {
        use crate::core::{
            blast_dynamic_align::DynamicSequence,
            blast_hits::blast_hsp_list_sort_by_score,
            blast_prelim::{GapSequence, GapSubject},
            blast_traceback::blast_traceback_from_hsp_list,
        };
        blast_hsp_list_sort_by_score(Some(&mut list));
        let mut scoring = None;
        let status = blast_scoring_parameters_new(
            Some(&self.scoring_options),
            &self.score_block,
            &mut scoring,
        );
        if status != 0 {
            return Err(NucleotideSetupError::Upstream {
                stage: "scoring parameters",
                status,
            });
        }
        let scoring = scoring.unwrap();
        // Standard-task scoring rows are constant for the search; keep their
        // independent owner instead of cloning rows for each traceback list.
        let mut gap = bind_gap_align(
            self.gap_align.take().expect("search-owned gapped scratch"),
            &self.traceback_matrix,
        );
        gap.dynamic.gap_x_dropoff = self.extension_parameters.gap_x_dropoff_final;
        let query_data = DynamicSequence {
            data: self.batch.sequence_with_sentinels(),
            offset: 1,
        };
        let query = GapSequence {
            sequence: Some(query_data),
            oof_sequence: None,
            length: self.batch.blastna_contexts().len() as i32,
        };
        assert!(expanded.len() >= 2);
        let target = GapSubject {
            sequence: DynamicSequence {
                data: expanded,
                offset: 1,
            },
            length: (expanded.len() - 2) as i32,
            frame: 1,
        };
        let status = blast_traceback_from_hsp_list(
            EBlastProgramType::Blastn,
            &mut list,
            &query,
            &query_data,
            &target,
            &self.query_info,
            &mut gap,
            &mut self.score_block,
            &scoring,
            &self.extension_parameters.options,
            &self.hit_parameters,
        );
        self.gap_align = Some(bind_gap_align(gap, &[]));
        if status != 0 {
            return Err(NucleotideSetupError::Upstream {
                stage: "traceback",
                status,
            });
        }
        Ok(list)
    }
}

/// Restore one NCBI4na subject directly into the original padded BLASTna
/// layout. Retaining this vector across an OID batch avoids an intermediate
/// whole-subject allocation and bounds capacity by the largest visited subject.
pub(crate) fn expand_nucleotide_subject(subject: &NucleotideSequence, expanded: &mut Vec<u8>) {
    use crate::core::blast_encoding::{K_NUCL_SENTINEL, NCBI4NA_TO_BLASTNA};
    expanded.clear();
    expanded.reserve(subject.len() + 2);
    expanded.push(K_NUCL_SENTINEL);
    expanded.extend(
        subject
            .ncbi4na()
            .iter()
            .map(|base| NCBI4NA_TO_BLASTNA[*base as usize]),
    );
    expanded.push(K_NUCL_SENTINEL);
}

// Same safe ownership reconstruction as the Kappa scratch binder: replacing
// only the borrowed matrix changes the lifetime without altering owned data.
fn bind_gap_align<'matrix>(
    gap: BlastGapAlignStruct<'_>,
    matrix: &'matrix [Vec<i32>],
) -> BlastGapAlignStruct<'matrix> {
    let dynamic = gap.dynamic;
    BlastGapAlignStruct {
        dynamic: BlastDynamicAlign {
            matrix,
            pssm: None,
            position_based: false,
            gap_x_dropoff: dynamic.gap_x_dropoff,
            dp_mem: dynamic.dp_mem,
            dp_mem_alloc: dynamic.dp_mem_alloc,
            state_struct: dynamic.state_struct,
            query_start: dynamic.query_start,
            query_stop: dynamic.query_stop,
            subject_start: dynamic.subject_start,
            subject_stop: dynamic.subject_stop,
            score: dynamic.score,
        },
        greedy_align_mem: gap.greedy_align_mem,
        edit_script: gap.edit_script,
        fwd_prelim_tback: gap.fwd_prelim_tback,
        rev_prelim_tback: gap.rev_prelim_tback,
        greedy_query_seed_start: gap.greedy_query_seed_start,
        greedy_subject_seed_start: gap.greedy_subject_seed_start,
    }
}
