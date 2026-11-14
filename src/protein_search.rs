//! Native database-wide standard tblastn search. Preliminary selection and
//! composition redo retain the original engine, collector and heap algorithms.
use crate::core::{
    blast_composition_redo::blast_compo_early_termination,
    blast_diagnostics::{BlastGappedStats, BlastUngappedStats},
    blast_hits::*,
    blast_parameters::blast_scoring_parameters_new,
    blast_prelim_evaluate::{s_blast_search_engine_core, PreliminaryContextResult},
    blast_program::EBlastProgramType,
    blast_query_info::BlastQueryInfo,
    compo_heap::*,
    hspfilter_collector::*,
};
use crate::protein_redo::{ProteinRedoContext, ProteinRedoError};
use crate::protein_setup::{
    ProteinInitialSearch, ProteinPreliminaryFrame, ProteinPreliminarySubject, ProteinSetupError,
};
use crate::{NucleotideDatabase, PreparedProteinBatch};
use std::fmt;
const PROGRAM: EBlastProgramType = EBlastProgramType::Tblastn;
/// Default source PSI_INCLUSION_ETHRESH used by redo even for ordinary tblastn.
pub const PSI_INCLUSION_ETHRESH: f64 = 0.002;
#[derive(Debug)]
pub enum ProteinSearchError {
    Setup(ProteinSetupError),
    Redo(ProteinRedoError),
    QuerySplit(crate::query_split::QuerySplitError),
    EmptyDatabase,
    DatabaseTooLarge,
    UnsupportedSubjectLength {
        ordinal: usize,
        length: usize,
        maximum: usize,
    },
    InvalidSubjectOrdinal(i32),
    MissingState(&'static str),
    Upstream {
        stage: &'static str,
        status: i32,
    },
}
impl fmt::Display for ProteinSearchError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Setup(e) => e.fmt(f),
            Self::Redo(e) => e.fmt(f),
            Self::QuerySplit(e) => e.fmt(f),
            Self::EmptyDatabase => write!(f, "tblastn requires a nonempty nucleotide database"),
            Self::DatabaseTooLarge => {
                write!(f, "database dimensions exceed the BLAST integer domain")
            }
            Self::UnsupportedSubjectLength {
                ordinal,
                length,
                maximum,
            } => write!(
                f,
                "subject {ordinal} length {length} exceeds the current unchunked limit {maximum}"
            ),
            Self::InvalidSubjectOrdinal(oid) => write!(f, "missing tblastn subject ordinal {oid}"),
            Self::MissingState(stage) => write!(f, "missing tblastn state: {stage}"),
            Self::Upstream { stage, status } => {
                write!(f, "tblastn {stage} failed with status {status}")
            }
        }
    }
}
impl std::error::Error for ProteinSearchError {}
impl From<ProteinSetupError> for ProteinSearchError {
    fn from(e: ProteinSetupError) -> Self {
        Self::Setup(e)
    }
}
impl From<ProteinRedoError> for ProteinSearchError {
    fn from(e: ProteinRedoError) -> Self {
        Self::Redo(e)
    }
}
impl From<crate::query_split::QuerySplitError> for ProteinSearchError {
    fn from(error: crate::query_split::QuerySplitError) -> Self {
        Self::QuerySplit(error)
    }
}
fn checked(status: i32, stage: &'static str) -> Result<(), ProteinSearchError> {
    if status == 0 {
        Ok(())
    } else {
        Err(ProteinSearchError::Upstream { stage, status })
    }
}
/// Native HSP coordinates remain context-local protein/translated-frame
/// coordinates. OIDs index subject_ids; contexts index query_info.
#[derive(Debug)]
pub struct ProteinSearchResults {
    pub results: Box<BlastHSPResults>,
    pub query_ids: Vec<String>,
    pub subject_ids: Vec<String>,
    pub query_info: BlastQueryInfo,
    pub ungapped_stats: BlastUngappedStats,
    pub gapped_stats: BlastGappedStats,
    pub preliminary_subject_lists: usize,
    pub traceback_subject_lists: usize,
    pub early_terminated_subject_lists: usize,
}
/// Source s_FillResultsFromCompoHeaps, preserving worst-first heap pops,
/// bounded HitListUpdate and the final source reverse-order operation.
pub fn s_fill_results_from_compo_heaps(
    results: &mut BlastHSPResults,
    heaps: &mut [BlastCompoHeap<Box<BlastHSPList>>],
    hitlist_size: i32,
) {
    for query in 0..results.num_queries as usize {
        results.hitlist_array[query] = Some(blast_hit_list_new(hitlist_size));
        let list = results.hitlist_array[query].as_mut().unwrap();
        while let Some(hsps) = blast_compo_heap_pop(&mut heaps[query]) {
            blast_hit_list_update(list, hsps);
        }
    }
    blast_hsp_results_reverse_order(results);
}
fn collect_protein_preliminary(
    search: &mut ProteinInitialSearch<'_>,
    database: &NucleotideDatabase,
) -> Result<Box<BlastHSPResults>, ProteinSearchError> {
    let query_info = search.query_info().clone();
    let options = search
        .hit_parameters
        .options
        .as_ref()
        .ok_or(ProteinSearchError::MissingState("hit-saving options"))?
        .clone();
    let collector = blast_hsp_collector_params_new(Some(&options), 2, true)
        .ok_or(ProteinSearchError::MissingState("preliminary collector"))?;
    let mut scoring = None;
    checked(
        blast_scoring_parameters_new(
            Some(&search.scoring_options),
            &search.score_block,
            &mut scoring,
        ) as i32,
        "scoring parameters",
    )?;
    let scoring = scoring.ok_or(ProteinSearchError::MissingState("scoring parameters"))?;
    let mut preliminary = blast_hsp_results_new(query_info.num_queries);
    for (oid, subject) in database.sequences().iter().enumerate() {
        let chunked = subject.len() / 3 > crate::core::blast_gapalign::MAX_DBSEQ_LEN as usize;
        let raw = if chunked {
            crate::protein_chunks::preliminary_chunked_protein_subject(search, subject, oid as i32)?
        } else {
            search.preliminary_subject(subject)?
        };
        let mut frames: Vec<Option<Box<BlastHSPList>>> = raw
            .frames
            .into_iter()
            .map(|frame| Some(frame.hsp_list))
            .collect();
        let lengths: Vec<i32> = [1i32, 2, 3, -1, -2, -3]
            .into_iter()
            .map(|frame| (subject.len() as i32 - frame.abs() + 1) / 3)
            .collect();
        let mut evaluated = None;
        checked(
            s_blast_search_engine_core(
                PROGRAM,
                &lengths,
                subject.len() as i32,
                oid as i32,
                &query_info,
                &scoring,
                &search.hit_parameters,
                &mut search.score_block,
                None,
                |frame| {
                    Ok(if chunked {
                        PreliminaryContextResult::Finished(frames[frame].take())
                    } else {
                        PreliminaryContextResult::GappedScore(frames[frame].take())
                    })
                },
                &mut evaluated,
            ) as i32,
            "preliminary evaluation",
        )?;
        if let Some(list) = evaluated {
            checked(
                s_blast_hsp_collector_run(Some(&mut preliminary), Some(&collector), Some(list)),
                "preliminary collector",
            )?;
            if let Some(low) = search.hit_parameters.low_score.as_mut() {
                for query in 0..preliminary.num_queries as usize {
                    if let Some(list) = preliminary.hitlist_array[query].as_ref() {
                        if list.heapified {
                            low[query] = low[query]
                                .max((options.low_score_perc * list.low_score as f64) as i32);
                        }
                    }
                }
            }
        }
    }
    Ok(preliminary)
}
// Native diagnostic aggregation across source query chunks. Counters describe
// executed chunk work; final traceback still runs once against the full query.
fn add_preliminary_stats(
    u: &mut BlastUngappedStats,
    g: &mut BlastGappedStats,
    incoming_u: &BlastUngappedStats,
    incoming_g: &BlastGappedStats,
) {
    u.lookup_hits = u.lookup_hits.wrapping_add(incoming_u.lookup_hits);
    u.num_seqs_lookup_hits = u
        .num_seqs_lookup_hits
        .wrapping_add(incoming_u.num_seqs_lookup_hits);
    u.init_extends = u.init_extends.wrapping_add(incoming_u.init_extends);
    u.good_init_extends = u
        .good_init_extends
        .wrapping_add(incoming_u.good_init_extends);
    u.num_seqs_passed = u.num_seqs_passed.wrapping_add(incoming_u.num_seqs_passed);
    g.seqs_ungapped_passed = g
        .seqs_ungapped_passed
        .wrapping_add(incoming_g.seqs_ungapped_passed);
    g.extensions = g.extensions.wrapping_add(incoming_g.extensions);
    g.good_extensions = g.good_extensions.wrapping_add(incoming_g.good_extensions);
    g.num_seqs_passed = g.num_seqs_passed.wrapping_add(incoming_g.num_seqs_passed);
}
/// Searches with standard BLOSUM62/11/1, word3/two-hit40, SEG and CBS2.
/// Whole-search parity and performance remain separate verification gates.
pub fn search_protein_database(
    batch: &PreparedProteinBatch,
    database: &NucleotideDatabase,
) -> Result<ProteinSearchResults, ProteinSearchError> {
    if database.is_empty() {
        return Err(ProteinSearchError::EmptyDatabase);
    }
    if database.len() > i32::MAX as usize || database.total_length() > i64::MAX as u64 {
        return Err(ProteinSearchError::DatabaseTooLarge);
    }
    let lengths: Vec<_> = (0..batch.query_info().num_queries as usize)
        .map(|i| batch.query_info().contexts[i].query_length as usize)
        .collect();
    let plan = crate::query_split::QuerySplitPlan::protein(&lengths)?;
    let mut search = if plan.is_split() {
        ProteinInitialSearch::new_for_traceback(batch, database)?
    } else {
        ProteinInitialSearch::new(batch, database)?
    };
    let query_info = search.query_info().clone();
    let options = search
        .hit_parameters
        .options
        .as_ref()
        .ok_or(ProteinSearchError::MissingState("hit-saving options"))?
        .clone();
    let (mut preliminary, ungapped_stats, gapped_stats) = if plan.is_split() {
        // CPrelimSearch::Run sets full-query effective spaces before setting up
        // any chunk. The unchanged vector intentionally retains local-prefix
        // indexing in each chunk; length adjustment is recomputed locally.
        let spaces: Vec<_> = query_info
            .contexts
            .iter()
            .map(|ctx| ctx.eff_searchsp)
            .collect();
        let mut combined = blast_hsp_results_new(query_info.num_queries);
        let mut ungapped = BlastUngappedStats::default();
        let mut gapped = BlastGappedStats::default();
        for chunk in 0..plan.num_chunks() {
            let chunk_batch = plan.prepare_protein_batch_chunk(chunk, batch)?;
            let mut chunk_search = match ProteinInitialSearch::new_with_effective_searchspaces(
                &chunk_batch,
                database,
                &spaces,
            ) {
                Ok(search) => search,
                // Original wrapper ignores only the all-invalid-context KA
                // error for one chunk, not arbitrary setup/search errors.
                Err(ProteinSetupError::Upstream {
                    stage: "ungapped statistics",
                    status: 1,
                }) => continue,
                Err(error) => return Err(error.into()),
            };
            let mut incoming = collect_protein_preliminary(&mut chunk_search, database)?;
            checked(
                crate::query_split::blast_hsp_stream_merge(
                    plan.split_block(),
                    chunk,
                    &mut incoming,
                    &mut combined,
                    PROGRAM,
                ),
                "split-query preliminary merge",
            )?;
            add_preliminary_stats(
                &mut ungapped,
                &mut gapped,
                chunk_search.stats(),
                chunk_search.gapped_stats(),
            );
        }
        (combined, ungapped, gapped)
    } else {
        let preliminary = collect_protein_preliminary(&mut search, database)?;
        (preliminary, *search.stats(), *search.gapped_stats())
    };
    // CBS streams enable sort_on_read: ReverseSort, then consume each query's
    // lists from the end. This differs from ordinary traceback's OID batching.
    checked(
        blast_hsp_results_reverse_sort(&mut preliminary) as i32,
        "CBS stream order",
    )?;
    let preliminary_subject_lists = preliminary
        .hitlist_array
        .iter()
        .flatten()
        .map(|list| list.hsplist_count as usize)
        .sum();
    let mut heaps: Vec<BlastCompoHeap<Box<BlastHSPList>>> = (0..query_info.num_queries)
        .map(|_| BlastCompoHeap::default())
        .collect();
    for heap in &mut heaps {
        checked(
            blast_compo_heap_initialize(heap, options.hitlist_size, PSI_INCLUSION_ETHRESH),
            "composition heap initialization",
        )?;
    }
    let (mut traceback_subject_lists, mut early_terminated_subject_lists) = (0, 0);
    let mut redo = ProteinRedoContext::new(&search)?;
    for query in 0..query_info.num_queries as usize {
        let Some(hitlist) = preliminary.hitlist_array[query].as_mut() else {
            continue;
        };
        while hitlist.hsplist_count > 0 {
            let index = hitlist.hsplist_count as usize - 1;
            let mut incoming = hitlist.hsplist_array[index]
                .take()
                .expect("selected HSP list");
            hitlist.hsplist_count -= 1;
            incoming.query_index = query as i32;
            if blast_compo_early_termination(incoming.best_evalue, &heaps, query_info.num_queries) {
                early_terminated_subject_lists += 1;
                continue;
            }
            let oid = incoming.oid;
            let subject = usize::try_from(oid)
                .ok()
                .and_then(|oid| database.sequence(oid))
                .ok_or(ProteinSearchError::InvalidSubjectOrdinal(oid))?;
            // EngineCore has already combined/sorted subject frames; redo reads the
            // frame stored in each HSP. The wrapper's frame field is not consulted.
            let selected = ProteinPreliminarySubject {
                frames: vec![ProteinPreliminaryFrame {
                    frame: 0,
                    hsp_list: incoming,
                }],
            };
            let mut redone = redo.redo_subject(subject, oid, selected)?;
            traceback_subject_lists += 1;
            if let Some(list) = redone.query_lists[query].take() {
                let best_evalue = list.best_evalue;
                let best_score = redone.best_scores[query];
                if best_evalue <= options.expect_value
                    && blast_compo_heap_would_insert(
                        &mut heaps[query],
                        best_evalue,
                        best_score,
                        oid,
                    )
                {
                    let mut discarded = None;
                    checked(
                        blast_compo_heap_insert(
                            &mut heaps[query],
                            list,
                            best_evalue,
                            best_score,
                            oid,
                            &mut discarded,
                        ),
                        "composition heap insertion",
                    )?;
                }
            }
        }
    }
    let mut results = blast_hsp_results_new(query_info.num_queries);
    s_fill_results_from_compo_heaps(&mut results, &mut heaps, options.hitlist_size);
    // BLAST_ComputeTraceback_MT's outer database tail follows Kappa's heap
    // drain: apply the source Evalue/score/OID comparator (including its
    // underflow equivalence), then the original final subject-list limit.
    checked(
        blast_hsp_results_sort_by_evalue(Some(&mut results)) as i32,
        "final result sorting",
    )?;
    crate::core::blast_result_pipeline::s_blast_prune_extra_hits(
        &mut results,
        options.hitlist_size,
    );
    Ok(ProteinSearchResults {
        results,
        query_ids: batch.identifiers().to_vec(),
        subject_ids: database
            .sequences()
            .iter()
            .map(|s| s.id().to_owned())
            .collect(),
        query_info,
        ungapped_stats,
        gapped_stats,
        preliminary_subject_lists,
        traceback_subject_lists,
        early_terminated_subject_lists,
    })
}
