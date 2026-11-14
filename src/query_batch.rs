//! Native query ownership around translated context, sequence and mask kernels.
//! A batch uses one standard nucleotide task and both strands of each query.
//!
//! ```
//! use btyper_blast::{NucleotideTask, PreparedNucleotideBatch, PreparedNucleotideQuery};
//! let queries = [PreparedNucleotideQuery::new(
//!     "query", b"AGTCCGATGCTACGATCGTGCATAGCGTACCT", NucleotideTask::Blastn,
//! ).unwrap()];
//! let batch = PreparedNucleotideBatch::new(&queries).unwrap();
//! assert_eq!(batch.query_info().contexts.len(), 2);
//! let seeds = batch.seed_index().unwrap();
//! assert_eq!(seeds.query_blastna(), batch.blastna_contexts());
//! ```

use crate::core::blast_encoding::EBlastEncoding;
use crate::core::blast_filter::{blast_complement_mask_locations, BlastMaskLoc};
use crate::core::blast_lookup::SSeqRange;
use crate::core::blast_program::EBlastProgramType;
use crate::core::blast_query_info::{
    blast_query_info_new, query_info_get_seq_buf_len, BlastQueryInfo,
};
use crate::core::blast_setup::{
    get_sequence_nucleotide_both_strands, s_query_info_set_context, ESentinelType,
};
use crate::core::blast_util::get_reverse_nucl_sequence;
use crate::seed::{NucleotideSeedIndex, SeedError};
use crate::{NucleotideTask, PreparedNucleotideQuery};
use std::fmt;

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum QueryBatchError {
    EmptyBatch,
    MixedTasks,
    CombinedSequenceTooLong,
}
impl fmt::Display for QueryBatchError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyBatch => write!(f, "query batch is empty"),
            Self::MixedTasks => write!(
                f,
                "all queries in a batch must use the same nucleotide task"
            ),
            Self::CombinedSequenceTooLong => write!(
                f,
                "combined query contexts exceed the BLAST i32 buffer limit"
            ),
        }
    }
}
impl std::error::Error for QueryBatchError {}

/// Assembled, soft-masked default nucleotide contexts. This is query/seed
/// preparation; extension and complete alignment search remain separate work.
#[derive(Debug)]
pub struct PreparedNucleotideBatch {
    identifiers: Vec<String>,
    task: NucleotideTask,
    query_info: BlastQueryInfo,
    sequence: Vec<u8>,
    lookup_regions: Vec<SSeqRange>,
}

impl PreparedNucleotideBatch {
    pub fn new(queries: &[PreparedNucleotideQuery]) -> Result<Self, QueryBatchError> {
        let Some(first) = queries.first() else {
            return Err(QueryBatchError::EmptyBatch);
        };
        let task = first.task();
        let mut required = 1_u64;
        for query in queries {
            if query.task() != task {
                return Err(QueryBatchError::MixedTasks);
            }
            required += 2 * (query.len() as u64 + 1);
            if required > i32::MAX as u64 {
                return Err(QueryBatchError::CombinedSequenceTooLong);
            }
        }
        let program = EBlastProgramType::Blastn;
        let mut query_info = blast_query_info_new(program, queries.len() as i32).unwrap();
        let mut max_length = 0;
        let mut min_length = i32::MAX as u32;
        for (index, query) in queries.iter().enumerate() {
            let length = query.len() as u32;
            max_length = max_length.max(length);
            min_length = min_length.min(length);
            s_query_info_set_context(&mut query_info, 2 * index as u32, length);
            s_query_info_set_context(&mut query_info, 2 * index as u32 + 1, length);
        }
        query_info.max_length = max_length;
        query_info.min_length = min_length;
        let buflen = query_info_get_seq_buf_len(&query_info) as usize;
        // SetupQueries allocates one extra zero byte beyond the defined buffer.
        let mut sequence = vec![0; buflen + 1];
        let mut masks = BlastMaskLoc {
            total_size: query_info.last_context + 1,
            seqloc_array: vec![Vec::new(); query_info.last_context as usize + 1],
        };
        for (index, query) in queries.iter().enumerate() {
            let reverse = get_reverse_nucl_sequence(query.ncbi4na(), query.len() as i32);
            let both = get_sequence_nucleotide_both_strands(
                query.ncbi4na(),
                &reverse[1..reverse.len() - 1],
                EBlastEncoding::Nucleotide,
                ESentinelType::Sentinels,
            );
            let offset = query_info.contexts[2 * index].query_offset as usize;
            sequence[offset..offset + both.len()].copy_from_slice(&both);
            let forward_masks: Vec<_> = query
                .dust_masks()
                .iter()
                .map(|(left, right)| SSeqRange {
                    left: *left as i32,
                    right: *right as i32,
                })
                .collect();
            // CORE expects masks in forward-query coordinates for both
            // contexts; its complement function reverses the minus list.
            masks.seqloc_array[2 * index] = forward_masks.clone();
            masks.seqloc_array[2 * index + 1] = forward_masks;
        }
        let mut lookup_regions = Vec::new();
        blast_complement_mask_locations(
            program,
            &query_info,
            Some(&mut masks),
            Some(&mut lookup_regions),
        );
        Ok(Self {
            identifiers: queries.iter().map(|query| query.id().to_owned()).collect(),
            task,
            query_info,
            sequence,
            lookup_regions,
        })
    }

    pub fn identifiers(&self) -> &[String] {
        &self.identifiers
    }
    pub fn task(&self) -> NucleotideTask {
        self.task
    }
    pub fn query_info(&self) -> &BlastQueryInfo {
        &self.query_info
    }
    /// Defined allocation prefix, including leading/inter-context/trailing
    /// sentinels. The extra zero allocation byte is excluded.
    pub fn sequence_with_sentinels(&self) -> &[u8] {
        &self.sequence[..self.sequence.len() - 1]
    }
    /// Sequence-block view: both endpoint sentinels are excluded, while
    /// separators between strands and queries retain their original offsets.
    pub fn blastna_contexts(&self) -> &[u8] {
        &self.sequence[1..self.sequence.len() - 2]
    }
    pub fn lookup_regions(&self) -> &[SSeqRange] {
        &self.lookup_regions
    }
    pub fn seed_index(&self) -> Result<NucleotideSeedIndex, SeedError> {
        NucleotideSeedIndex::from_segments(self.blastna_contexts(), &self.lookup_regions, self.task)
    }

    /// Build only the requested local query interval. Recover inherited masks
    /// from the full forward lookup complement before running local DUST.
    pub(crate) fn prepare_query_range(
        &self,
        index: usize,
        range: std::ops::Range<usize>,
    ) -> Result<PreparedNucleotideQuery, crate::QueryError> {
        use crate::core::blast_encoding::BLASTNA_TO_NCBI4NA;
        let info = &self.query_info.contexts[2 * index];
        assert!(range.start < range.end && range.end <= info.query_length as usize);
        let context_start = info.query_offset;
        let context_end = context_start + info.query_length;
        let start = context_start + range.start as i32;
        let end = context_start + range.end as i32;
        let residues = self.blastna_contexts()[start as usize..end as usize]
            .iter()
            .map(|&code| BLASTNA_TO_NCBI4NA[code as usize])
            .collect::<Vec<_>>();
        let mut masks = Vec::new();
        let mut next = context_start;
        let mut add_mask = |left: i32, right: i32| {
            if left < end && right >= start && left <= right {
                masks.push((
                    (left.max(start) - start) as u32,
                    (right.min(end - 1) - start) as u32,
                ));
            }
        };
        for region in &self.lookup_regions {
            if region.left >= context_end
                || region.right < context_start
                || region.left > region.right
            {
                continue;
            }
            let left = region.left.max(context_start);
            let right = region.right.min(context_end - 1);
            if next < left {
                add_mask(next, left - 1);
            }
            next = right + 1;
        }
        if next < context_end {
            add_mask(next, context_end - 1);
        }
        PreparedNucleotideQuery::from_ncbi4na_with_additional_masks(
            &self.identifiers[index],
            &residues,
            self.task,
            &masks,
        )
    }
}
