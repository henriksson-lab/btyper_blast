//! Native database-wide default blastn/megablast pipeline with in-memory results.
//! The original preliminary collector selects candidates before traceback;
//! final result heaps and Evalue ordering apply after selected HSPs are traced.
use crate::core::{
    blast_diagnostics::BlastUngappedStats,
    blast_hits::*,
    blast_parameters::blast_scoring_parameters_new,
    blast_program::EBlastProgramType,
    blast_query_info::BlastQueryInfo,
    blast_result_pipeline::{s_blast_prune_extra_hits, s_sort_hsp_list_by_oid},
    hspfilter_collector::{blast_hsp_collector_params_new, s_blast_hsp_collector_run},
};
use crate::nucleotide_setup::{
    expand_nucleotide_subject, NucleotideInitialSearch, NucleotideSetupError,
};
use crate::query_split::{blast_hsp_stream_merge, QuerySplitError, QuerySplitPlan};
use crate::{NucleotideDatabase, NucleotideTask, PreparedNucleotideBatch};
use std::fmt;
#[derive(Debug)]
pub enum NucleotideSearchError {
    Setup(NucleotideSetupError),
    QuerySplit(QuerySplitError),
    EmptyDatabase,
    DatabaseTooLarge,
    InvalidSubjectOrdinal(i32),
    MissingState(&'static str),
    Upstream { stage: &'static str, status: i32 },
}
impl fmt::Display for NucleotideSearchError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Setup(e) => e.fmt(f),
            Self::QuerySplit(e) => e.fmt(f),
            Self::EmptyDatabase => write!(f, "nucleotide search requires a nonempty database"),
            Self::DatabaseTooLarge => {
                write!(f, "database dimensions exceed native BLAST integer limits")
            }
            Self::InvalidSubjectOrdinal(oid) => {
                write!(f, "BLAST result references missing subject ordinal {oid}")
            }
            Self::MissingState(stage) => {
                write!(f, "BLAST {stage} did not initialize required state")
            }
            Self::Upstream { stage, status } => {
                write!(f, "BLAST {stage} failed with status {status}")
            }
        }
    }
}
impl std::error::Error for NucleotideSearchError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Setup(e) => Some(e),
            Self::QuerySplit(e) => Some(e),
            _ => None,
        }
    }
}
impl From<NucleotideSetupError> for NucleotideSearchError {
    fn from(error: NucleotideSetupError) -> Self {
        Self::Setup(error)
    }
}
impl From<QuerySplitError> for NucleotideSearchError {
    fn from(error: QuerySplitError) -> Self {
        Self::QuerySplit(error)
    }
}
fn checked(status: i32, stage: &'static str) -> Result<(), NucleotideSearchError> {
    if status == 0 {
        Ok(())
    } else {
        Err(NucleotideSearchError::Upstream { stage, status })
    }
}
/// HSP context indexes refer to query_info; each subject OID indexes subject_ids.
/// Query offsets retain native per-context coordinates (reverse-strand contexts
/// are not reformatted into the textual BLAST coordinate convention).
#[derive(Debug)]
pub struct NucleotideSearchResults {
    pub results: Box<BlastHSPResults>,
    pub query_ids: Vec<String>,
    pub subject_ids: Vec<String>,
    pub query_info: BlastQueryInfo,
    pub task: NucleotideTask,
    pub ungapped_stats: BlastUngappedStats,
    pub preliminary_subject_lists: usize,
    pub traceback_subject_lists: usize,
}
/// Source preliminary engine and collector for one complete query batch.
/// Query splitting invokes this for each local batch before stream merging;
/// final traceback is deliberately absent from this stage.
fn collect_preliminary_database(
    search: &mut NucleotideInitialSearch<'_>,
    database: &NucleotideDatabase,
) -> Result<Box<BlastHSPResults>, NucleotideSearchError> {
    let query_info = search.query_info().clone();
    let program = EBlastProgramType::Blastn;
    let options = search
        .hit_parameters
        .options
        .as_ref()
        .ok_or(NucleotideSearchError::MissingState("hit-saving options"))?
        .clone();
    let collector = blast_hsp_collector_params_new(Some(&options), 0, true)
        .ok_or(NucleotideSearchError::MissingState("collector parameters"))?;
    let mut scoring = None;
    checked(
        blast_scoring_parameters_new(
            Some(&search.scoring_options),
            &search.score_block,
            &mut scoring,
        ) as i32,
        "scoring parameters",
    )?;
    let scoring = scoring.ok_or(NucleotideSearchError::MissingState("scoring parameters"))?;
    let mut preliminary = blast_hsp_results_new(query_info.num_queries);
    for (ordinal, subject) in database.sequences().iter().enumerate() {
        let raw = if subject.len() <= crate::core::blast_gapalign::MAX_DBSEQ_LEN as usize {
            crate::core::blast_prelim_evaluate::PreliminaryContextResult::GappedScore(Some(
                search.preliminary_subject(subject)?,
            ))
        } else {
            let list = preliminary_chunked_subject(search, subject, ordinal as i32)?;
            crate::core::blast_prelim_evaluate::PreliminaryContextResult::Finished(list)
        };
        let mut raw = Some(raw);
        let mut list = None;
        checked(
            crate::core::blast_prelim_evaluate::s_blast_search_engine_core(
                program,
                &[subject.len() as i32],
                subject.len() as i32,
                ordinal as i32,
                &query_info,
                &scoring,
                &search.hit_parameters,
                &mut search.score_block,
                None,
                |_| Ok(raw.take().expect("one nucleotide subject context")),
                &mut list,
            ) as i32,
            "preliminary evaluation",
        )?;
        if let Some(list) = list {
            checked(
                s_blast_hsp_collector_run(Some(&mut preliminary), Some(&collector), Some(list)),
                "preliminary collector",
            )?;
            // Original engine feeds the preliminary heap's threshold back to
            // GetGappedScore only when a per-query low-score array is enabled.
            if let Some(low_scores) = search.hit_parameters.low_score.as_mut() {
                for query_index in 0..preliminary.num_queries as usize {
                    if let Some(hitlist) = preliminary.hitlist_array[query_index].as_ref() {
                        if hitlist.heapified {
                            low_scores[query_index] = low_scores[query_index]
                                .max((options.low_score_perc * hitlist.low_score as f64) as i32);
                        }
                    }
                }
            }
        }
    }
    Ok(preliminary)
}

/// Search a prepared query batch using its standard blastn or megablast task.
/// This assembles native kernels; full supported-result parity and performance
/// remain separate verification requirements.
pub fn search_nucleotide_database(
    batch: &PreparedNucleotideBatch,
    database: &NucleotideDatabase,
) -> Result<NucleotideSearchResults, NucleotideSearchError> {
    if database.is_empty() {
        return Err(NucleotideSearchError::EmptyDatabase);
    }
    if database.len() > i32::MAX as usize || database.total_length() > i64::MAX as u64 {
        return Err(NucleotideSearchError::DatabaseTooLarge);
    }
    let lengths: Vec<_> = (0..batch.query_info().num_queries as usize)
        .map(|query| batch.query_info().contexts[2 * query].query_length as usize)
        .collect();
    let plan = QuerySplitPlan::nucleotide(batch.task(), &lengths)?;
    let mut search = if plan.is_split() {
        NucleotideInitialSearch::new_for_traceback(batch, database)?
    } else {
        NucleotideInitialSearch::new(batch, database)?
    };
    let query_info = search.query_info().clone();
    let program = EBlastProgramType::Blastn;
    let options = search
        .hit_parameters
        .options
        .as_ref()
        .ok_or(NucleotideSearchError::MissingState("hit-saving options"))?
        .clone();
    let (mut preliminary, ungapped_stats) = if !plan.is_split() {
        let preliminary = collect_preliminary_database(&mut search, database)?;
        (preliminary, *search.stats())
    } else {
        // CQuerySplitter re-runs setup on each local sequence piece, while
        // SetEffectiveSearchSpace passes the full vector unchanged. Source
        // CalcEffLengths indexes that prefix by LOCAL context, not global map.
        let effective_spaces: Vec<_> = query_info
            .contexts
            .iter()
            .map(|context| context.eff_searchsp)
            .collect();
        let mut combined = blast_hsp_results_new(query_info.num_queries);
        let mut stats = BlastUngappedStats::default();
        for chunk in 0..plan.num_chunks() {
            let local_batch = plan.prepare_nucleotide_batch_chunk(chunk, batch)?;
            let mut local_search = match NucleotideInitialSearch::new_with_effective_searchspaces(
                &local_batch,
                database,
                &effective_spaces,
            ) {
                Ok(search) => search,
                // Original prelim_stage.cpp ignores only a split chunk whose
                // every context lacks calculable ungapped KA parameters.
                Err(NucleotideSetupError::Upstream {
                    stage: "ungapped statistics",
                    status: 1,
                }) => continue,
                Err(error) => return Err(error.into()),
            };
            let mut local = collect_preliminary_database(&mut local_search, database)?;
            let local_stats = local_search.stats();
            stats.lookup_hits = stats.lookup_hits.wrapping_add(local_stats.lookup_hits);
            stats.num_seqs_lookup_hits = stats
                .num_seqs_lookup_hits
                .wrapping_add(local_stats.num_seqs_lookup_hits);
            stats.init_extends = stats.init_extends.wrapping_add(local_stats.init_extends);
            stats.good_init_extends = stats
                .good_init_extends
                .wrapping_add(local_stats.good_init_extends);
            stats.num_seqs_passed = stats
                .num_seqs_passed
                .wrapping_add(local_stats.num_seqs_passed);
            checked(
                blast_hsp_stream_merge(
                    plan.split_block(),
                    chunk,
                    &mut local,
                    &mut combined,
                    program,
                ),
                "split query preliminary stream merge",
            )?;
        }
        (combined, stats)
    };
    // BlastHSPStreamClose transfers the collector's selected lists, assigns
    // query_index, sorts decreasing OID, then BatchRead consumes from the end.
    let mut selected = Vec::new();
    for query_index in 0..preliminary.num_queries as usize {
        if let Some(hitlist) = preliminary.hitlist_array[query_index].as_mut() {
            for entry in &mut hitlist.hsplist_array[..hitlist.hsplist_count as usize] {
                if let Some(mut list) = entry.take() {
                    list.query_index = query_index as i32;
                    selected.push(list);
                }
            }
            hitlist.hsplist_count = 0;
        }
    }
    selected.sort_by(|a, b| s_sort_hsp_list_by_oid(a, b));
    let preliminary_subject_lists = selected.len();
    let mut traceback_subject_lists = 0;
    let mut results = blast_hsp_results_new(query_info.num_queries);
    // Source BlastHSPStreamBatchRead groups all query lists for one subject.
    // Selected lists already have source OID order; only the owner lifetime
    // changes here, with one buffer bounded by the largest visited subject.
    let mut expanded = Vec::new();
    let mut expanded_oid = None;
    while let Some(list) = selected.pop() {
        let oid = list.oid;
        let subject = usize::try_from(oid)
            .ok()
            .and_then(|oid| database.sequence(oid))
            .ok_or(NucleotideSearchError::InvalidSubjectOrdinal(oid))?;
        let query_index = list.query_index;
        if expanded_oid != Some(oid) {
            expand_nucleotide_subject(subject, &mut expanded);
            expanded_oid = Some(oid);
        }
        let mut traced = search.traceback_preliminary_subject_expanded(&expanded, list)?;
        traceback_subject_lists += 1;
        traced.oid = oid;
        traced.query_index = query_index;
        if traced.hspcnt > 0 {
            checked(
                blast_hsp_results_insert_hsp_list(&mut results, Some(traced), options.hitlist_size)
                    as i32,
                "final result insertion",
            )?;
        }
    }
    checked(
        blast_hsp_results_sort_by_evalue(Some(&mut results)) as i32,
        "final result sorting",
    )?;
    s_blast_prune_extra_hits(&mut results, options.hitlist_size);
    Ok(NucleotideSearchResults {
        results,
        query_ids: batch.identifiers().to_vec(),
        subject_ids: database
            .sequences()
            .iter()
            .map(|s| s.id().to_owned())
            .collect(),
        query_info,
        task: batch.task(),
        ungapped_stats,
        preliminary_subject_lists,
        traceback_subject_lists,
    })
}

/// Native binding of the original OneContext chunk loop. Original packed data
/// remain borrowed through each chunk; ambiguous placeholders are never regenerated.
fn preliminary_chunked_subject(
    search: &mut NucleotideInitialSearch<'_>,
    subject: &crate::NucleotideSequence,
    oid: i32,
) -> Result<Option<Box<BlastHSPList>>, NucleotideSearchError> {
    use crate::core::{
        blast_prelim_evaluate::s_blast_search_engine_one_context, blast_subject_chunk::*,
    };
    let mut backup = SubjectSplitStruct::default();
    let mut chunk = SubjectChunk::default();
    s_backup_subject(subject.len() as i32, &mut backup, &mut chunk);
    let mut combined = None;
    while s_get_next_subject_chunk(&mut chunk, &mut backup, true, DBSEQ_CHUNK_OVERLAP)
        == SUBJECT_SPLIT_OK
    {
        let raw = search.preliminary_subject_chunk(subject, backup.offset, chunk.length)?;
        let mut list = None;
        checked(
            s_blast_search_engine_one_context(
                EBlastProgramType::Blastn,
                oid,
                || Ok(Some(raw)),
                &mut list,
            ) as i32,
            "chunk preliminary evaluation",
        )?;
        if let Some(list) = list.as_mut() {
            blast_hsp_list_adjust_offsets(list, backup.offset);
        }
        let overlap = if backup.offset == 0 {
            0
        } else {
            DBSEQ_CHUNK_OVERLAP
        };
        checked(
            blast_hsp_lists_merge(
                &mut list,
                &mut combined,
                i32::MAX,
                &[backup.offset],
                i32::MIN,
                overlap,
                true,
                false,
            ) as i32,
            "subject chunk merging",
        )?;
    }
    Ok(combined)
}
