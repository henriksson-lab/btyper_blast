//! Native ownership/dispatch adapter for the standard nucleotide seed phase.
//! Inputs are already assembled BLASTna contexts and included lookup regions.
//! Returned pairs are candidate seeds; extension and final search results are
//! not implemented here. This module does not promise a full search API.

use crate::core::blast_lookup::SSeqRange;
use crate::core::blast_nalookup::*;
use crate::core::blast_nascan::*;
use crate::core::blast_options::{QuerySetUpOptions, SBlastFilterOptions};
use crate::core::blast_program::EBlastProgramType;
use crate::core::lookup_util::estimate_num_table_entries;
use crate::{NucleotideSequence, NucleotideTask};
use std::fmt;

const OFFSET_ARRAY_SIZE: usize = 4096;

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum SeedError {
    InvalidQuery,
    InvalidRegions,
    ConstructionFailed(i16),
    InvalidOffset,
    InsufficientCapacity { minimum: usize },
}

impl fmt::Display for SeedError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidQuery => write!(
                f,
                "seed query must be nonempty BLASTna within the i32 length limit"
            ),
            Self::InvalidRegions => write!(
                f,
                "lookup regions must be ordered, disjoint inclusive query intervals"
            ),
            Self::ConstructionFailed(code) => {
                write!(f, "lookup construction failed with status {code}")
            }
            Self::InvalidOffset => write!(
                f,
                "scan offset must be nonnegative and aligned to the lookup stride"
            ),
            Self::InsufficientCapacity { minimum } => {
                write!(f, "seed buffer requires at least {minimum} pairs")
            }
        }
    }
}
impl std::error::Error for SeedError {}

#[derive(Debug)]
enum Table {
    Small(BlastSmallNaLookupTable),
    Ordinary(BlastNaLookupTable),
    Megablast(BlastMBLookupTable),
}

/// Standard-preset lookup construction and bounded subject scanning. Regions
/// retain full-context offsets; query masking must already be reflected in them.
/// Empty markers `[right+1, right]` emitted by BLAST mask setup are retained.
#[derive(Debug)]
pub struct NucleotideSeedIndex {
    table: Table,
    query: Vec<u8>,
    compressed_query: Option<Vec<u8>>,
}

impl NucleotideSeedIndex {
    pub fn from_segments(
        query: &[u8],
        regions: &[SSeqRange],
        task: NucleotideTask,
    ) -> Result<Self, SeedError> {
        if query.is_empty()
            || query.len() > i32::MAX as usize
            || query.iter().any(|byte| *byte > 15)
        {
            return Err(SeedError::InvalidQuery);
        }
        let mut previous_right = -1;
        for region in regions {
            if region.left < 0
                || region.right < -1
                || region.left > query.len() as i32
                || region.right >= query.len() as i32
                || region.left > region.right + 1
                || region.left <= previous_right
            {
                return Err(SeedError::InvalidRegions);
            }
            previous_right = region.right;
        }
        let options = LookupTableOptions {
            word_size: match task {
                NucleotideTask::Megablast => 28,
                NucleotideTask::Blastn => 11,
            },
            program_number: EBlastProgramType::Blastn,
            mb_template_length: 0,
            db_filter: false,
            stride: 0,
            max_db_word_count: 0,
        };
        let query_options = QuerySetUpOptions {
            filtering_options: Some(SBlastFilterOptions { mask_at_hash: true }),
            filter_string: None,
        };
        let mut max_offset = 0;
        let entries = estimate_num_table_entries(regions, &mut max_offset);
        let mut width = 0;
        let kind = blast_choose_na_lookup_table(&options, entries, max_offset, &mut width);
        let mut compressed_query = None;
        let table = match kind {
            ELookupTableType::SmallNaLookupTable => {
                match blast_small_na_lookup_table_new(query, regions, &options, width) {
                    Ok((lookup, compressed)) => {
                        compressed_query = Some(compressed);
                        Table::Small(lookup)
                    }
                    Err(_) => Table::Ordinary(blast_na_lookup_table_new(
                        query,
                        regions,
                        &options,
                        Some(&query_options),
                        width,
                    )),
                }
            }
            ELookupTableType::NaLookupTable => Table::Ordinary(blast_na_lookup_table_new(
                query,
                regions,
                &options,
                Some(&query_options),
                width,
            )),
            ELookupTableType::MBLookupTable => Table::Megablast(
                blast_mb_lookup_table_new(
                    query,
                    regions,
                    &options,
                    Some(&query_options),
                    entries,
                    width,
                )
                .map_err(SeedError::ConstructionFailed)?,
            ),
            _ => unreachable!("standard nucleotide presets selected an unsupported table"),
        };
        Ok(Self {
            table,
            query: query.to_vec(),
            compressed_query,
        })
    }

    pub fn query_blastna(&self) -> &[u8] {
        &self.query
    }
    pub fn compressed_query(&self) -> Option<&[u8]> {
        self.compressed_query.as_deref()
    }
    pub fn table_type(&self) -> ELookupTableType {
        match &self.table {
            Table::Small(_) => ELookupTableType::SmallNaLookupTable,
            Table::Ordinary(_) => ELookupTableType::NaLookupTable,
            Table::Megablast(_) => ELookupTableType::MBLookupTable,
        }
    }
    pub fn lookup_width(&self) -> i32 {
        match &self.table {
            Table::Small(t) => t.lut_word_length,
            Table::Ordinary(t) => t.lut_word_length,
            Table::Megablast(t) => t.lut_word_length,
        }
    }
    pub fn scan_step(&self) -> i32 {
        match &self.table {
            Table::Small(t) => t.scan_step,
            Table::Ordinary(t) => t.scan_step,
            Table::Megablast(t) => t.scan_step,
        }
    }
    pub fn longest_chain(&self) -> i32 {
        match &self.table {
            Table::Small(t) => t.longest_chain,
            Table::Ordinary(t) => t.longest_chain,
            Table::Megablast(t) => t.longest_chain,
        }
    }
    /// Same standard buffer allowance as upstream `GetOffsetArraySize`.
    pub fn recommended_capacity(&self) -> usize {
        OFFSET_ARRAY_SIZE + self.longest_chain() as usize
    }
    pub fn minimum_capacity(&self) -> usize {
        self.longest_chain() as usize + usize::from(matches!(&self.table, Table::Megablast(_)))
    }

    /// Scan from offset zero, then resume using the updated offset. This
    /// returns at most the buffer length and preserves original hit ordering.
    /// A subject shorter than the lookup width returns zero without scanning.
    pub fn scan_into(
        &self,
        subject: &NucleotideSequence,
        next_offset: &mut i32,
        pairs: &mut [BlastOffsetPair],
    ) -> Result<usize, SeedError> {
        if *next_offset < 0 || *next_offset % self.scan_step() != 0 {
            return Err(SeedError::InvalidOffset);
        }
        let last = subject.len() as i32 - self.lookup_width();
        if *next_offset > last {
            return Ok(0);
        }
        let minimum = self.minimum_capacity();
        if pairs.len() < minimum {
            return Err(SeedError::InsufficientCapacity { minimum });
        }
        let capacity = pairs.len().min(i32::MAX as usize) as i32;
        let mut range = [*next_offset, last];
        let packed = subject.packed_for_scanning();
        let view = PackedNucleotideSubject {
            sequence: packed,
            mask_type: ESubjectMaskingType::NoSubjMasking,
        };
        let hits = match &self.table {
            Table::Small(t) => match self.scan_step() {
                4 => s_blast_small_na_scan_subject_8_4(t, packed, pairs, capacity, &mut range),
                21 => {
                    s_blast_small_na_scan_subject_8_1_mod4(t, packed, pairs, capacity, &mut range)
                }
                _ => unreachable!(),
            },
            Table::Ordinary(t) => {
                if self.scan_step() == 4 {
                    s_blast_na_scan_subject_8_4(t, &view, pairs, capacity, &mut range)
                } else {
                    s_blast_na_scan_subject_any(t, &view, pairs, capacity, &mut range)
                }
            }
            Table::Megablast(t) => match (self.lookup_width(), self.scan_step()) {
                (10, 2) => s_mb_scan_subject_10_2(t, &view, pairs, capacity, &mut range),
                (11, 1) => s_mb_scan_subject_11_1_mod4(t, &view, pairs, capacity, &mut range),
                (11, 18) => s_mb_scan_subject_11_2_mod4(t, &view, pairs, capacity, &mut range),
                (12, 17) => s_mb_scan_subject_any(t, &view, pairs, capacity, &mut range),
                _ => unreachable!(),
            },
        };
        *next_offset = range[0];
        Ok(hits as usize)
    }
}

impl NucleotideSeedIndex {
    pub(crate) fn lookup_wrap(&self) -> crate::core::na_ungapped::NucleotideLookupWrap<'_> {
        use crate::core::na_ungapped::*;
        match &self.table {
            Table::Small(t) => NucleotideLookupWrap {
                lut: NucleotideLookup::Small(t),
                lookup_callback: s_small_na_lookup,
            },
            Table::Ordinary(t) => NucleotideLookupWrap {
                lut: NucleotideLookup::Ordinary(t),
                lookup_callback: s_na_lookup,
            },
            Table::Megablast(t) => NucleotideLookupWrap {
                lut: NucleotideLookup::Megablast(t),
                lookup_callback: s_mb_lookup,
            },
        }
    }
    pub(crate) fn scan_callback(&self) -> crate::core::na_ungapped::NucleotideScanCallback {
        use crate::core::na_ungapped::NucleotideScanCallback;
        match &self.table {
            Table::Small(t) => NucleotideScanCallback::Small(if t.scan_step == 4 {
                s_blast_small_na_scan_subject_8_4
            } else {
                s_blast_small_na_scan_subject_8_1_mod4
            }),
            Table::Ordinary(t) => NucleotideScanCallback::Ordinary(if t.scan_step == 4 {
                s_blast_na_scan_subject_8_4
            } else {
                s_blast_na_scan_subject_any
            }),
            Table::Megablast(t) => {
                NucleotideScanCallback::Megablast(match (t.lut_word_length, t.scan_step) {
                    (10, 2) => s_mb_scan_subject_10_2,
                    (11, 1) => s_mb_scan_subject_11_1_mod4,
                    (11, 18) => s_mb_scan_subject_11_2_mod4,
                    (12, 17) => s_mb_scan_subject_any,
                    _ => unreachable!("standard nucleotide scan selection"),
                })
            }
        }
    }
}
