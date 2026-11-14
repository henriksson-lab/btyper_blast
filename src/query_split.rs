//! Native ownership adapter for the ordinary-query CQuerySplitter path.
//! Chunk coordinates count biological query letters, not context sentinels.
//! Preliminary results must be merged before full-query traceback/CBS.
use crate::core::{
    blast_program::EBlastProgramType,
    blast_query_info::{blast_query_info_new, BlastQueryInfo},
    blast_query_split::{
        split_query_blk_new, split_query_calculate_num_chunks, split_query_get_chunk_size,
        x_compute_chunk_ranges, QuerySplitProgram, SQueryChunkBoundary, SSplitQueryBlk,
    },
    blast_setup::s_query_info_set_context,
};
use crate::NucleotideTask;
use crate::{
    PreparedNucleotideBatch, PreparedNucleotideQuery, PreparedProteinBatch, PreparedProteinQuery,
};
use std::{fmt, ops::Range};

#[derive(Debug, Clone, Eq, PartialEq)]
pub enum QuerySplitError {
    EmptyBatch,
    EmptyQuery(usize),
    SequenceTooLong,
    InvalidChunk(usize),
    QueryDimensions,
    Nucleotide(crate::query::QueryError),
    NucleotideBatch(crate::query_batch::QueryBatchError),
    Protein(crate::protein_query::ProteinQueryError),
}
impl fmt::Display for QuerySplitError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyBatch => write!(f, "query splitting requires a nonempty batch"),
            Self::EmptyQuery(i) => write!(f, "query {i} has no residues"),
            Self::SequenceTooLong => write!(f, "query split dimensions exceed BLAST offsets"),
            Self::InvalidChunk(i) => write!(f, "query split chunk {i} is out of range"),
            Self::QueryDimensions => write!(f, "query dimensions differ from the query split plan"),
            Self::Nucleotide(error) => error.fmt(f),
            Self::NucleotideBatch(error) => error.fmt(f),
            Self::Protein(error) => error.fmt(f),
        }
    }
}
impl std::error::Error for QuerySplitError {}

/// Native replacement for an ObjectMgr split query's sequence interval.
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct QueryPiece {
    pub query_index: usize,
    /// Half-open coordinates in the original forward query.
    pub range: Range<usize>,
}
#[derive(Debug, Clone)]
pub struct QuerySplitPlan {
    program: EBlastProgramType,
    lengths: Vec<usize>,
    chunk_size: usize,
    split_block: Box<SSplitQueryBlk>,
    pieces: Vec<Vec<QueryPiece>>,
}

/// s_SetSplitQuerySeqInterval, lowering CSeq_interval to a half-open range.
pub fn s_set_split_query_seq_interval(
    chunk: SQueryChunkBoundary,
    query_range: SQueryChunkBoundary,
    query_offset: usize,
) -> Range<usize> {
    assert!(chunk.left < query_range.right && query_range.left < chunk.right);
    let qstart = chunk.left as i64 - query_range.left as i64;
    let qend = chunk.right as i64 - query_range.right as i64;
    let start = qstart.max(0) as usize + query_offset;
    let end = if qend >= 0 {
        (query_range.right - query_range.left) as usize + query_offset
    } else {
        (chunk.right - query_range.left) as usize + query_offset
    };
    start..end
}

/// CQuerySplitter::x_ComputeQueryIndicesForChunks, ordinary full queries.
/// ObjectMgr scopes/IDs are represented by the original query index.
pub fn x_compute_query_indices_for_chunks(
    split: &mut SSplitQueryBlk,
    lengths: &[usize],
) -> Vec<Vec<QueryPiece>> {
    let mut query_ranges = Vec::with_capacity(lengths.len());
    let mut start = 0i32;
    for &length in lengths {
        let end = start + length as i32;
        query_ranges.push(SQueryChunkBoundary {
            left: start,
            right: end,
        });
        start = end;
    }
    let mut pieces = vec![Vec::new(); split.num_chunks as usize];
    for chunk_num in 0..split.num_chunks as usize {
        let chunk = split.chunk_bounds[chunk_num];
        for (qindex, &query_range) in query_ranges.iter().enumerate() {
            if chunk.left >= query_range.right || query_range.left >= chunk.right {
                continue;
            }
            split.chunk_query_map[chunk_num].push(qindex as u32);
            pieces[chunk_num].push(QueryPiece {
                query_index: qindex,
                range: s_set_split_query_seq_interval(chunk, query_range, 0),
            });
        }
    }
    pieces
}

/// Ordinary protein-query and both-strand nucleotide branches of
/// CQuerySplitter::x_ComputeQueryContextsForChunks. Translated queries are out
/// of scope: tblastn translates the subject, not the query.
pub fn x_compute_query_contexts_for_chunks(split: &mut SSplitQueryBlk, program: EBlastProgramType) {
    let contexts = match program {
        EBlastProgramType::Blastn => 2,
        EBlastProgramType::Tblastn => 1,
        _ => panic!("unsupported query split program"),
    };
    for chunk_num in 0..split.num_chunks as usize {
        for &query in &split.chunk_query_map[chunk_num] {
            for ctx in 0..contexts {
                split.chunk_ctx_map[chunk_num].push((contexts * query + ctx) as i32);
            }
        }
    }
}

/// CContextTranslator::GetAbsoluteContext.
pub fn get_absolute_context(split: &SSplitQueryBlk, chunk_num: usize, context: i32) -> i32 {
    split.chunk_ctx_map[chunk_num][context as usize]
}
/// CContextTranslator::GetContextInChunk.
pub fn get_context_in_chunk(
    split: &SSplitQueryBlk,
    chunk_num: usize,
    absolute_context: i32,
) -> i32 {
    split.chunk_ctx_map[chunk_num]
        .iter()
        .position(|&ctx| ctx == absolute_context)
        .map_or(-1, |ctx| ctx as i32)
}
/// CContextTranslator::GetStartingChunk; signed reverse traversal replaces the
/// original size_t-underflow/signed-cast termination condition.
pub fn get_starting_chunk(split: &SSplitQueryBlk, curr_chunk: usize, context: i32) -> i32 {
    let absolute_context = get_absolute_context(split, curr_chunk, context);
    if absolute_context == -1 {
        return -1;
    }
    let mut retval = curr_chunk;
    for previous in (0..curr_chunk).rev() {
        if get_context_in_chunk(split, previous, absolute_context) == -1 {
            break;
        }
        retval = previous;
    }
    retval as i32
}
pub fn s_is_plus_strand(info: &BlastQueryInfo, context: i32) -> bool {
    info.contexts[context as usize].frame >= 0
}
pub fn s_get_absolute_context_length(
    chunk_qinfo: &[BlastQueryInfo],
    chunk_num: i32,
    split: &SSplitQueryBlk,
    absolute_context: i32,
) -> usize {
    if chunk_num < 0 {
        return 0;
    }
    let pos = get_context_in_chunk(split, chunk_num as usize, absolute_context);
    if pos != -1 {
        chunk_qinfo[chunk_num as usize].contexts[pos as usize].query_length as usize
    } else {
        0
    }
}
/// Exact accumulation order of the nontranslated CQuerySplitter offset path.
pub fn x_compute_context_offsets_non_translated_queries(
    split: &mut SSplitQueryBlk,
    global_qinfo: &BlastQueryInfo,
    chunk_qinfo: &[BlastQueryInfo],
) {
    let overlap_size = split.chunk_overlap_sz;
    for chunk_num in 0..split.num_chunks as usize {
        let info = &chunk_qinfo[chunk_num];
        for _ in 0..info.first_context {
            split.chunk_offset_map[chunk_num].push(i32::MAX as u32);
        }
        for ctx in info.first_context..=info.last_context {
            let starting_chunk = get_starting_chunk(split, chunk_num, ctx);
            let absolute_context = get_absolute_context(split, chunk_num, ctx);
            if absolute_context == -1 || starting_chunk == -1 {
                assert!(!info.contexts[ctx as usize].is_valid);
                split.chunk_offset_map[chunk_num].push(i32::MAX as u32);
                continue;
            }
            let correction;
            if s_is_plus_strand(info, ctx) {
                let mut sum = 0usize;
                let mut c = chunk_num as i32;
                while c != starting_chunk {
                    let previous =
                        s_get_absolute_context_length(chunk_qinfo, c - 1, split, absolute_context);
                    let current =
                        s_get_absolute_context_length(chunk_qinfo, c, split, absolute_context);
                    let overlap = overlap_size.min(current);
                    sum += previous - overlap.min(previous);
                    c -= 1;
                }
                correction = sum;
            } else {
                let mut subtrahend = 0usize;
                let mut c = chunk_num as i32;
                while c >= starting_chunk && c >= 0 {
                    let previous =
                        s_get_absolute_context_length(chunk_qinfo, c - 1, split, absolute_context);
                    let current =
                        s_get_absolute_context_length(chunk_qinfo, c, split, absolute_context);
                    let overlap = overlap_size.min(current);
                    subtrahend += current - overlap.min(previous);
                    c -= 1;
                }
                correction = global_qinfo.contexts[absolute_context as usize].query_length as usize
                    - subtrahend;
            }
            assert!(info.contexts[ctx as usize].is_valid);
            split.chunk_offset_map[chunk_num].push(correction as u32);
        }
    }
}

fn query_info_for_lengths(program: EBlastProgramType, lengths: &[usize]) -> BlastQueryInfo {
    let mut info = blast_query_info_new(program, lengths.len() as i32).unwrap();
    let ncontexts = if program == EBlastProgramType::Blastn {
        2
    } else {
        1
    };
    for (qindex, &length) in lengths.iter().enumerate() {
        for context in 0..ncontexts {
            s_query_info_set_context(
                &mut info,
                (qindex * ncontexts + context) as u32,
                length as u32,
            );
        }
    }
    info.min_length = *lengths.iter().min().unwrap() as u32;
    info.max_length = *lengths.iter().max().unwrap() as u32;
    info
}
impl QuerySplitPlan {
    pub fn nucleotide(task: NucleotideTask, lengths: &[usize]) -> Result<Self, QuerySplitError> {
        Self::new(
            match task {
                NucleotideTask::Megablast => QuerySplitProgram::Megablast,
                NucleotideTask::Blastn => QuerySplitProgram::Blastn,
            },
            lengths,
        )
    }
    pub fn protein(lengths: &[usize]) -> Result<Self, QuerySplitError> {
        Self::new(QuerySplitProgram::Tblastn, lengths)
    }
    pub fn new(kind: QuerySplitProgram, lengths: &[usize]) -> Result<Self, QuerySplitError> {
        if lengths.is_empty() {
            return Err(QuerySplitError::EmptyBatch);
        }
        let mut total = 0usize;
        for (i, &length) in lengths.iter().enumerate() {
            if length == 0 {
                return Err(QuerySplitError::EmptyQuery(i));
            }
            total = total
                .checked_add(length)
                .ok_or(QuerySplitError::SequenceTooLong)?;
        }
        if total > i32::MAX as usize || lengths.len() > i32::MAX as usize / 2 {
            return Err(QuerySplitError::SequenceTooLong);
        }
        let program = if kind == QuerySplitProgram::Tblastn {
            EBlastProgramType::Tblastn
        } else {
            EBlastProgramType::Blastn
        };
        let mut chunk_size = split_query_get_chunk_size(kind);
        let num_chunks =
            split_query_calculate_num_chunks(program, &mut chunk_size, total, lengths.len());
        let mut split_block = split_query_blk_new(num_chunks, true).unwrap();
        // Source Split skips geometry for an unsplit batch. Native identity
        // maps make the adapter usable uniformly, without changing searches.
        x_compute_chunk_ranges(&mut split_block, chunk_size, total, program);
        let pieces = x_compute_query_indices_for_chunks(&mut split_block, lengths);
        x_compute_query_contexts_for_chunks(&mut split_block, program);
        let global = query_info_for_lengths(program, lengths);
        let infos = pieces
            .iter()
            .map(|chunk| {
                let lens: Vec<_> = chunk.iter().map(|p| p.range.len()).collect();
                query_info_for_lengths(program, &lens)
            })
            .collect::<Vec<_>>();
        x_compute_context_offsets_non_translated_queries(&mut split_block, &global, &infos);
        Ok(Self {
            program,
            lengths: lengths.to_vec(),
            chunk_size,
            split_block,
            pieces,
        })
    }
    pub fn is_split(&self) -> bool {
        self.split_block.num_chunks > 1
    }
    pub fn num_chunks(&self) -> usize {
        self.pieces.len()
    }
    pub fn chunk_size(&self) -> usize {
        self.chunk_size
    }
    pub fn split_block(&self) -> &SSplitQueryBlk {
        &self.split_block
    }
    pub fn pieces(&self, chunk: usize) -> Result<&[QueryPiece], QuerySplitError> {
        self.pieces
            .get(chunk)
            .map(Vec::as_slice)
            .ok_or(QuerySplitError::InvalidChunk(chunk))
    }
    pub fn program(&self) -> EBlastProgramType {
        self.program
    }
    pub fn query_lengths(&self) -> &[usize] {
        &self.lengths
    }

    /// Bounded ownership path: materialize only this chunk, with no temporary
    /// copy or filtering pass over all full-query residue buffers.
    pub fn prepare_nucleotide_batch_chunk(
        &self,
        chunk: usize,
        batch: &PreparedNucleotideBatch,
    ) -> Result<PreparedNucleotideBatch, QuerySplitError> {
        if self.program != EBlastProgramType::Blastn
            || batch.identifiers().len() != self.lengths.len()
            || self.lengths.iter().enumerate().any(|(i, &length)| {
                batch.query_info().contexts[2 * i].query_length as usize != length
            })
        {
            return Err(QuerySplitError::QueryDimensions);
        }
        let queries = self
            .pieces(chunk)?
            .iter()
            .map(|piece| {
                batch
                    .prepare_query_range(piece.query_index, piece.range.clone())
                    .map_err(QuerySplitError::Nucleotide)
            })
            .collect::<Result<Vec<_>, _>>()?;
        PreparedNucleotideBatch::new(&queries).map_err(QuerySplitError::NucleotideBatch)
    }
    pub fn prepare_protein_batch_chunk(
        &self,
        chunk: usize,
        batch: &PreparedProteinBatch,
    ) -> Result<PreparedProteinBatch, QuerySplitError> {
        if self.program != EBlastProgramType::Tblastn
            || batch.identifiers().len() != self.lengths.len()
            || self
                .lengths
                .iter()
                .enumerate()
                .any(|(i, &length)| batch.query_info().contexts[i].query_length as usize != length)
        {
            return Err(QuerySplitError::QueryDimensions);
        }
        let queries = self
            .pieces(chunk)?
            .iter()
            .map(|piece| {
                batch
                    .prepare_query_range(piece.query_index, piece.range.clone())
                    .map_err(QuerySplitError::Protein)
            })
            .collect::<Result<Vec<_>, _>>()?;
        PreparedProteinBatch::new(&queries).map_err(QuerySplitError::Protein)
    }

    /// Prepare chunk DUST from local residues and inherited full-query masks.
    /// This deliberately does not treat an already-filtered chunk as raw input.
    pub fn prepare_nucleotide_chunk(
        &self,
        chunk: usize,
        queries: &[PreparedNucleotideQuery],
    ) -> Result<PreparedNucleotideBatch, QuerySplitError> {
        if self.program != EBlastProgramType::Blastn
            || queries.len() != self.lengths.len()
            || queries
                .iter()
                .zip(&self.lengths)
                .any(|(query, &length)| query.len() != length)
        {
            return Err(QuerySplitError::QueryDimensions);
        }
        let mut prepared = Vec::new();
        for piece in self.pieces(chunk)? {
            let query = &queries[piece.query_index];
            let start = piece.range.start as u32;
            let end = piece.range.end as u32;
            let inherited = query
                .dust_masks()
                .iter()
                .filter_map(|&(left, right)| {
                    if left < end && right >= start {
                        Some((left.max(start) - start, right.min(end - 1) - start))
                    } else {
                        None
                    }
                })
                .collect::<Vec<_>>();
            prepared.push(
                PreparedNucleotideQuery::from_ncbi4na_with_additional_masks(
                    query.id(),
                    &query.ncbi4na()[piece.range.clone()],
                    query.task(),
                    &inherited,
                )
                .map_err(QuerySplitError::Nucleotide)?,
            );
        }
        PreparedNucleotideBatch::new(&prepared).map_err(QuerySplitError::NucleotideBatch)
    }
    /// Protein masks originate in CORE, not ObjectMgr factory user masks, so
    /// source setup recomputes SEG from each local unmasked protein sequence.
    pub fn prepare_protein_chunk(
        &self,
        chunk: usize,
        queries: &[PreparedProteinQuery],
    ) -> Result<PreparedProteinBatch, QuerySplitError> {
        if self.program != EBlastProgramType::Tblastn
            || queries.len() != self.lengths.len()
            || queries
                .iter()
                .zip(&self.lengths)
                .any(|(query, &length)| query.len() != length)
        {
            return Err(QuerySplitError::QueryDimensions);
        }
        let mut prepared = Vec::new();
        for piece in self.pieces(chunk)? {
            let query = &queries[piece.query_index];
            prepared.push(
                PreparedProteinQuery::from_ncbistdaa(
                    query.id(),
                    &query.ncbistdaa()[piece.range.clone()],
                )
                .map_err(QuerySplitError::Protein)?,
            );
        }
        PreparedProteinBatch::new(&prepared).map_err(QuerySplitError::Protein)
    }
}

/// SplitQuery_SetEffectiveSearchSpace after the original full-query calculator
/// has populated full_query_info. Keep the ENTIRE vector: original chunk setup
/// indexes its local context into this unchanged prefix, not a projected map.
pub fn split_query_set_effective_search_space(
    options: &mut crate::core::blast_effective_lengths::BlastEffectiveLengthsOptions,
    full_query_info: &BlastQueryInfo,
) {
    if options.num_searchspaces != 0
        && options
            .searchsp_eff
            .as_ref()
            .expect("search-space array required")[0]
            != 0
    {
        return;
    }
    let spaces: Vec<_> = (0..=full_query_info.last_context)
        .map(|index| full_query_info.contexts[index as usize].eff_searchsp)
        .collect();
    options.num_searchspaces = spaces.len() as i32;
    options.searchsp_eff = Some(spaces);
}

/// BlastHSPStreamMerge's results-level operation. The caller finalizes its
/// collectors before entry; native ownership replaces the two stream objects.
/// All score ties, heap limits and overlap stitching remain in HitListMerge.
pub fn blast_hsp_stream_merge(
    split: &SSplitQueryBlk,
    chunk_num: usize,
    old_results: &mut crate::core::blast_hits::BlastHSPResults,
    combined: &mut crate::core::blast_hits::BlastHSPResults,
    program: EBlastProgramType,
) -> i32 {
    use crate::core::{
        blast_hits::blast_hsp_list_sort_by_score,
        blast_query_split::blast_hit_list_merge,
        blast_util::{blast_context_to_frame, blast_get_number_of_contexts},
    };
    let contexts_per_query = blast_get_number_of_contexts(program) as i32;
    let queries = &split.chunk_query_map[chunk_num];
    let contexts = &split.chunk_ctx_map[chunk_num];
    let offsets = &split.chunk_offset_map[chunk_num];
    for i in 0..old_results.num_queries as usize {
        let global_query = queries[i] as usize;
        let Some(hitlist) = old_results.hitlist_array[i].as_mut() else {
            continue;
        };
        let mut split_points = [-1i32; 6];
        for j in 0..contexts_per_query as usize {
            let local_context = i * contexts_per_query as usize + j;
            if contexts[local_context] >= 0 {
                split_points[contexts[local_context] as usize % contexts_per_query as usize] =
                    offsets[local_context] as i32;
            }
        }
        for entry in &mut hitlist.hsplist_array[..hitlist.hsplist_count as usize] {
            let list = entry.as_mut().unwrap();
            for hsp in &mut list.hsp_array[..list.hspcnt as usize] {
                let hsp = hsp.as_mut().unwrap();
                let local_context = hsp.context as usize;
                hsp.context = contexts[local_context];
                hsp.query.offset += offsets[local_context] as i32;
                hsp.query.end += offsets[local_context] as i32;
                hsp.query.gapped_start += offsets[local_context] as i32;
                hsp.query.frame = blast_context_to_frame(program, hsp.context as u32) as i16;
            }
            list.query_index = global_query as i32;
        }
        blast_hit_list_merge(
            &mut old_results.hitlist_array[i],
            &mut combined.hitlist_array[global_query],
            contexts_per_query,
            &split_points,
            split.chunk_overlap_sz as i32,
            split.gapped_merge,
        );
    }
    for hitlist in combined.hitlist_array[..combined.num_queries as usize]
        .iter_mut()
        .flatten()
    {
        for list in hitlist.hsplist_array[..hitlist.hsplist_count as usize]
            .iter_mut()
            .flatten()
        {
            blast_hsp_list_sort_by_score(Some(list));
        }
    }
    0
}
