//! Native standard tblastn CBS2 final redo. The core numerical and alignment
//! functions retain their pinned BLAST+ 2.17.0+ correspondence; this owner binds
//! those functions to borrowed in-memory subjects and preliminary HSP owners.
use crate::core::{
    blast_composition_adjustment::*,
    blast_composition_redo::*,
    blast_gapalign::BlastGapAlignStruct,
    blast_hits::*,
    blast_hsp_evaluate::*,
    blast_hsp_link::blast_link_hsps,
    blast_kappa::*,
    blast_parameters::{blast_scoring_parameters_new, BlastScoringParameters},
    blast_program::EBlastProgramType,
    blast_stat::BlastScoreBlk,
    blast_target_translation::{blast_hsp_get_target_translation, blast_target_translation_new},
    genetic_code::STANDARD_GENETIC_CODE,
};
use crate::{
    protein_setup::{ProteinInitialSearch, ProteinPreliminarySubject},
    NucleotideSequence,
};
use std::fmt;
const PROGRAM: EBlastProgramType = EBlastProgramType::Tblastn;
#[derive(Debug)]
pub enum ProteinRedoError {
    InvalidSubject,
    MissingState(&'static str),
    UnsupportedMode(i32),
    Upstream { stage: &'static str, status: i32 },
}
impl fmt::Display for ProteinRedoError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidSubject => write!(f, "invalid tblastn subject length or ordinal"),
            Self::MissingState(stage) => write!(f, "missing tblastn redo state: {stage}"),
            Self::UnsupportedMode(mode) => write!(
                f,
                "tblastn redo requires standard composition mode2, got {mode}"
            ),
            Self::Upstream { stage, status } => {
                write!(f, "tblastn {stage} failed with status {status}")
            }
        }
    }
}
impl std::error::Error for ProteinRedoError {}
fn checked(status: i32, stage: &'static str) -> Result<(), ProteinRedoError> {
    if status == 0 {
        Ok(())
    } else {
        Err(ProteinRedoError::Upstream { stage, status })
    }
}
#[derive(Debug)]
pub struct ProteinRedoSubject {
    pub query_lists: Vec<Option<Box<BlastHSPList>>>,
    /// Original scaled score used as the composition-heap tie key.
    pub best_scores: Vec<i32>,
}
impl ProteinRedoSubject {
    pub fn hsp_count(&self) -> usize {
        self.query_lists
            .iter()
            .flatten()
            .map(|list| list.hspcnt as usize)
            .sum()
    }
}

/// Source blast_kappa.c s_HitlistReapContained. The significance order is
/// retained; equal-or-lower-scoring HSPs contained in an earlier same-frame
/// HSP are removed, then the original live prefix is compacted.
pub fn s_hitlist_reap_contained(array: &mut [Option<Box<BlastHSP>>], count: &mut i32) {
    let old = *count as usize;
    for read in 1..old {
        for back in 0..read {
            let Some(one) = array[read].as_ref() else {
                break;
            };
            let Some(two) = array[back].as_ref() else {
                continue;
            };
            if two.query.frame == one.query.frame
                && two.subject.frame == one.subject.frame
                && two.query.offset <= one.query.offset
                && two.query.end >= one.query.offset
                && two.subject.offset <= one.subject.offset
                && two.subject.end >= one.subject.offset
                && two.query.offset <= one.query.end
                && two.query.end >= one.query.end
                && two.subject.offset <= one.subject.end
                && two.subject.end >= one.subject.end
                && one.score <= two.score
            {
                array[read] = None;
            }
        }
    }
    let mut write = 0;
    for read in 0..old {
        if array[read].is_some() {
            if write != read {
                array[write] = array[read].take();
            }
            write += 1;
        }
    }
    *count = write as i32;
    for slot in &mut array[write..old] {
        *slot = None;
    }
}
/// Source s_ComputeNumIdentities, specialized to ordinary tblastn. Original
/// partial-target translation caches serve each HSP; only the identity-helper
/// call sees local coordinates, which are restored before returning.
pub fn s_compute_num_identities(
    query: &[u8],
    query_info: &crate::core::blast_query_info::BlastQueryInfo,
    subject: &[u8],
    list: Option<&mut BlastHSPList>,
    options: &crate::core::blast_options::BlastScoringOptions,
    sbp: &BlastScoreBlk,
) -> i16 {
    let Some(list) = list else {
        return 0;
    };
    assert!(!options.is_ooframe);
    let mut target =
        blast_target_translation_new(subject, &STANDARD_GENETIC_CODE, PROGRAM, options.is_ooframe);
    for slot in &mut list.hsp_array[..list.hspcnt as usize] {
        let hsp = slot.as_mut().expect("live HSP prefix");
        let translation = blast_hsp_get_target_translation(&mut target, Some(hsp))
            .expect("source target translation for an ordinary HSP");
        hsp.subject.offset -= translation.start_shift;
        hsp.subject.end -= translation.start_shift;
        let offset = query_info.contexts[hsp.context as usize].query_offset as usize;
        let status = blast_hsp_get_num_identities_and_positives(
            Some(&query[offset..]),
            Some(translation.residues),
            hsp,
            options,
            None,
            Some(sbp),
        );
        hsp.subject.offset += translation.start_shift;
        hsp.subject.end += translation.start_shift;
        assert_eq!(status, 0);
    }
    0
}
/// Source s_HitlistEvaluateAndPurge, restricted to tblastn. Protein-only
/// composition P-value combination and diagnostic sequence-source IDs are
/// omitted because the source does not apply those branches to tblastn.
pub fn s_hitlist_evaluate_and_purge(
    best_score: &mut i32,
    best_evalue: &mut f64,
    list: &mut BlastHSPList,
    subject_length: i32,
    query_info: &crate::core::blast_query_info::BlastQueryInfo,
    sbp: &mut BlastScoreBlk,
    hit: &crate::core::blast_parameters::BlastHitSavingParameters,
) -> i32 {
    *best_evalue = f64::MAX;
    *best_score = 0;
    let status = if hit.do_sum_stats {
        blast_link_hsps(
            PROGRAM,
            Some(list),
            query_info,
            subject_length,
            sbp,
            hit.link_hsp_params
                .as_ref()
                .expect("source sum-statistics parameters"),
            true,
        )
    } else {
        blast_hsp_list_get_evalues(
            PROGRAM,
            query_info,
            subject_length,
            Some(list),
            true,
            false,
            sbp,
            0.0,
            1.0,
        )
    };
    if status == 0 {
        blast_hsp_list_reap_by_evalue(
            Some(list),
            hit.options
                .as_ref()
                .expect("source hit options")
                .expect_value,
        );
        if list.hspcnt > 0 {
            *best_evalue = list.best_evalue;
            *best_score = list.hsp_array[0].as_ref().unwrap().score;
        }
    }
    if status == 0 {
        0
    } else {
        -1
    }
}
/// Original per-search/per-thread CBS owners. Matrix, composition workspace,
/// query words and traceback scratch survive across subject calls, as they do
/// in Blast_RedoAlignmentCore. Borrowing the immutable initial-search owner
/// keeps shared query/hit-option lifetimes explicit without self references.
pub struct ProteinRedoContext<'search, 'batch> {
    search: &'search ProteinInitialSearch<'batch>,
    scaled: BlastScoreBlk,
    scoring: BlastScoringParameters,
    params: BlastRedoAlignParams,
    compo_queries: Vec<BlastCompoQueryInfo>,
    workspace: BlastCompositionWorkspace,
    matrix: Vec<i32>,
    callback_matrix_rows: Vec<Vec<i32>>,
    alignments: Vec<Vec<BlastCompoAlignment<BlastKappaAlignmentContext>>>,
    gap_align: Option<BlastGapAlignStruct<'static>>,
}
impl<'search, 'batch> ProteinRedoContext<'search, 'batch> {
    pub fn new(search: &'search ProteinInitialSearch<'batch>) -> Result<Self, ProteinRedoError> {
        let mode = search.extension_parameters.options.composition_based_stats;
        if mode != 2 {
            return Err(ProteinRedoError::UnsupportedMode(mode));
        }
        let query_info = search.query_info();
        let num_contexts = (query_info.last_context + 1) as usize;
        let mut scaled = search.score_block.clone();
        let mut scoring = None;
        checked(
            blast_scoring_parameters_new(Some(&search.scoring_options), &scaled, &mut scoring)
                as i32,
            "scoring parameters",
        )?;
        let mut scoring: BlastScoringParameters =
            scoring.ok_or(ProteinRedoError::MissingState("scoring parameters"))?;
        let scale = SCALING_FACTOR as f64;
        s_rescale_search(&mut scaled, &mut scoring, num_contexts as i32, scale);
        let params = s_get_align_params(
            &scaled,
            &scoring,
            query_info,
            &search.hit_parameters,
            &search.extension_parameters,
            scale,
        )
        .ok_or(ProteinRedoError::MissingState(
            "composition redo parameters",
        ))?;
        let compo_queries = s_get_query_info(search.batch().encoded_contexts(), query_info, false);
        let mut workspace = blast_composition_workspace_new();
        checked(
            blast_composition_workspace_init(&mut workspace, b"BLOSUM62"),
            "composition workspace",
        )?;
        let alignments = (0..num_contexts).map(|_| Vec::new()).collect();
        let matrix = scaled.matrix.iter().flatten().copied().collect();
        Ok(Self {
            search,
            scaled,
            scoring,
            params,
            compo_queries,
            workspace,
            matrix,
            callback_matrix_rows: vec![vec![0; 28]; 28],
            alignments,
            gap_align: None,
        })
    }
    /// Run the original scaled-score per-match stages using this thread's
    /// existing CBS state. Returned HSP scores are normalized after E-values.
    pub fn redo_subject(
        &mut self,
        subject: &NucleotideSequence,
        oid: i32,
        mut preliminary: ProteinPreliminarySubject,
    ) -> Result<ProteinRedoSubject, ProteinRedoError> {
        if oid < 0 || subject.len() > i32::MAX as usize || subject.len() < 3 {
            return Err(ProteinRedoError::InvalidSubject);
        }
        let Self {
            search,
            scaled,
            scoring,
            params,
            compo_queries,
            workspace,
            matrix,
            callback_matrix_rows,
            alignments,
            gap_align,
        } = self;
        let query_info = search.query_info();
        let num_contexts = (query_info.last_context + 1) as usize;
        let scale = SCALING_FACTOR as f64;
        let subjects = [BlastKappaSubject {
            sequence_start: subject.ncbi4na(),
            length: subject.len() as i32,
            genetic_code: &STANDARD_GENETIC_CODE,
        }];
        let mut callbacks = BlastKappaCallbacks {
            subjects: &subjects,
            scoring,
            gap_align: gap_align.take(),
            matrix_rows: callback_matrix_rows,
        };
        let result = (|| {
            let matching = BlastCompoMatchingSequence {
                index: 0,
                length: subject.len() as i32,
            };
            let options = search
                .hit_parameters
                .options
                .as_ref()
                .ok_or(ProteinRedoError::MissingState("hit saving options"))?;
            let mut incoming: Vec<Vec<Option<Box<BlastHSP>>>> =
                (0..num_contexts).map(|_| Vec::new()).collect();
            for frame in &mut preliminary.frames {
                for slot in &mut frame.hsp_list.hsp_array[..frame.hsp_list.hspcnt as usize] {
                    let hsp = slot.take().expect("live preliminary HSP");
                    let context = hsp.context as usize;
                    if context >= num_contexts {
                        return Err(ProteinRedoError::MissingState("preliminary context"));
                    }
                    incoming[context].push(Some(hsp));
                }
                frame.hsp_list.hspcnt = 0;
            }
            let mut output: Vec<Option<Box<BlastHSPList>>> =
                (0..num_contexts).map(|_| None).collect();
            let mut best_scores = vec![0; num_contexts];
            for (context, mut hsps) in incoming.into_iter().enumerate() {
                if hsps.is_empty() {
                    continue;
                }
                // Merged preliminary frame lists have the same original score comparator.
                hsps.sort_by(|a, b| score_compare_hs_ps(a.as_deref(), b.as_deref()));
                let count = hsps.len() as i32;
                let mut groups =
                    s_result_hsp_to_distinct_align(&mut hsps, count, context as i32, scale);
                let incoming = std::mem::take(&mut groups[0]);
                let karlin = scaled
                    .kbp_gap
                    .as_ref()
                    .and_then(|blocks| blocks[context])
                    .ok_or(ProteinRedoError::MissingState("scaled gapped statistics"))?;
                let (mut pvalue, mut lambda_ratio) = (-1.0, 0.0);
                checked(
                    blast_redo_one_match(
                        alignments,
                        params,
                        &incoming,
                        count,
                        karlin.lambda,
                        &matching,
                        -1,
                        compo_queries,
                        num_contexts as i32,
                        matrix,
                        28,
                        workspace,
                        &mut pvalue,
                        search.extension_parameters.options.unified_p,
                        &mut lambda_ratio,
                        &mut callbacks,
                    ),
                    "composition redo",
                )?;
                let mut list = blast_hsp_list_new(0);
                checked(
                    s_hsp_list_from_distinct_alignments(
                        Some(&mut list),
                        &mut alignments[context],
                        oid,
                        0,
                    ),
                    "redo HSP conversion",
                )?;
                if list.hspcnt > 1 {
                    s_hitlist_reap_contained(&mut list.hsp_array, &mut list.hspcnt);
                }
                let (mut best_score, mut best_evalue) = (0, f64::MAX);
                checked(
                    s_hitlist_evaluate_and_purge(
                        &mut best_score,
                        &mut best_evalue,
                        &mut list,
                        subject.len() as i32,
                        query_info,
                        scaled,
                        &search.hit_parameters,
                    ),
                    "final E-values and reaping",
                )?;
                if best_evalue > options.expect_value || best_evalue.is_nan() {
                    continue;
                }
                best_scores[context] = best_score;
                s_hsp_list_normalize_scores(&mut list, karlin.lambda, karlin.log_k, scale);
                // The original score block matrix points at the last adjusted matrix.
                for (row, values) in scaled.matrix.iter_mut().zip(matrix.chunks_exact(28)) {
                    row.copy_from_slice(values);
                }
                checked(
                    s_compute_num_identities(
                        // The tblastn branch in s_ComputeNumIdentities uses
                        // query_blk->sequence (SEG-masked), whereas other
                        // programs use query_nomask. Display metrics are
                        // computed separately from the original residues.
                        search.batch().encoded_contexts(),
                        query_info,
                        subject.ncbi4na(),
                        Some(&mut list),
                        &search.scoring_options,
                        scaled,
                    ) as i32,
                    "identities and positives",
                )?;
                output[context] = Some(list);
            }
            Ok(ProteinRedoSubject {
                query_lists: output,
                best_scores,
            })
        })();
        // The callback binds only local subject and scoring borrows. Its matrix
        // borrow has already been removed by the core callback adapter, so all
        // traceback arenas safely return to this reusable owner even on failure.
        *gap_align = callbacks.gap_align.take();
        result
    }
}
/// Convenience API for a single match. Database-wide callers should construct
/// ProteinRedoContext once and reuse redo_subject for each selected subject.
pub fn redo_protein_subject(
    search: &ProteinInitialSearch<'_>,
    subject: &NucleotideSequence,
    oid: i32,
    preliminary: ProteinPreliminarySubject,
) -> Result<ProteinRedoSubject, ProteinRedoError> {
    ProteinRedoContext::new(search)?.redo_subject(subject, oid, preliminary)
}
