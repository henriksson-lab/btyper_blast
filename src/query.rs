//! Query input adapters for the supported standard nucleotide tasks.
//!
//! DUST masks are kept separately from sequence residues: standard soft masking
//! suppresses lookup seeds while extension/scoring still sees the original bases.
//!
//! ```
//! use btyper_blast::{NucleotideTask, PreparedNucleotideQuery};
//!
//! let query = PreparedNucleotideQuery::new(
//!     "query_1", b"ACGTACGTACGTACGT", NucleotideTask::Megablast,
//! ).unwrap();
//! assert!(query.soft_masking());
//! assert_eq!(query.len(), 16);
//! ```

use crate::core::blast_encoding::{NCBI4NA_TO_BLASTNA, NCBI4NA_TO_IUPACNA};
use crate::core::blast_util::get_reverse_nucl_sequence;
use crate::database::nucleotide_value;
use crate::dustmask::symdust::{CSymDustMasker, DEFAULT_LEVEL, DEFAULT_LINKER, DEFAULT_WINDOW};
use std::fmt;

/// Seq_loc.cpp x_MergeRanges restricted to finite, unfuzzed forward ranges
/// and fStrand_Ignore|fMerge_All. No IDs or ObjectMgr objects are required.
pub(crate) fn x_merge_ranges(first: &mut (u32, u32), second: (u32, u32)) -> bool {
    if first.0 <= second.0 && first.1 >= second.1 {
        return true;
    }
    if first.0 >= second.0 && first.1 <= second.1 {
        *first = second;
        return true;
    }
    if first.0 <= second.1 && second.0 <= first.1 {
        first.0 = first.0.min(second.0);
        first.1 = first.1.max(second.1);
        return true;
    }
    if first.1 + 1 == second.0 {
        first.1 = second.1;
        return true;
    }
    false
}
/// Seq_loc.cpp x_RangesToSeq_loc, its single-ID, forward-strand sorted range
/// branch. Vector storage replaces Seq-loc nodes and omitted fuzz fields.
pub(crate) fn x_ranges_to_seq_loc(mut ranges: Vec<(u32, u32)>) -> Vec<(u32, u32)> {
    ranges.sort_by_key(|range| range.0);
    let mut output = Vec::new();
    let mut last = None;
    for range in ranges {
        if let Some(previous) = last.as_mut() {
            if x_merge_ranges(previous, range) {
                continue;
            }
            output.push(*previous);
        }
        last = Some(range);
    }
    if let Some(last) = last {
        output.push(last);
    }
    output
}
/// dust_filter.cpp s_CombineDustMasksWithUserProvidedMasks after DUST and
/// coordinate mapping. This scope has local finite ranges, one ID, no fuzz.
pub(crate) fn s_combine_dust_masks_with_user_provided_masks(
    dust_masks: &[(u32, u32)],
    inherited_masks: &[(u32, u32)],
) -> Vec<(u32, u32)> {
    if dust_masks.is_empty() {
        return inherited_masks.to_vec();
    }
    let mut combined = inherited_masks.to_vec();
    combined.extend_from_slice(dust_masks);
    x_ranges_to_seq_loc(combined)
}

/// Plain `blastn` defaults to Megablast. The traditional blastn task is distinct.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum NucleotideTask {
    #[default]
    Megablast,
    Blastn,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum QueryError {
    EmptyIdentifier,
    EmptySequence,
    SequenceTooLong,
    InvalidResidue { offset: usize, byte: u8 },
    InvalidNcbi4na { offset: usize, value: u8 },
    InvalidMask { left: u32, right: u32 },
}

impl fmt::Display for QueryError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyIdentifier => write!(f, "query identifier is empty"),
            Self::EmptySequence => write!(f, "query contains no residues"),
            Self::SequenceTooLong => write!(f, "query exceeds the BLAST i32 length domain"),
            Self::InvalidResidue { offset, byte } => write!(
                f,
                "invalid query nucleotide byte {byte:#04x} at offset {offset}"
            ),
            Self::InvalidNcbi4na { offset, value } => {
                write!(f, "invalid query NCBI4na value {value} at offset {offset}")
            }
            Self::InvalidMask { left, right } => {
                write!(f, "query mask interval {left}..={right} is out of range")
            }
        }
    }
}

impl std::error::Error for QueryError {}

/// Prepared two-strand query data, excluding leading/trailing sentinels. DUST
/// intervals use inclusive, zero-based coordinates on the corresponding strand.
/// This prepares inputs; it is not a completed search engine.
#[derive(Debug)]
pub struct PreparedNucleotideQuery {
    id: String,
    task: NucleotideTask,
    ncbi4na: Vec<u8>,
    forward_blastna: Vec<u8>,
    reverse_blastna: Vec<u8>,
    dust_masks: Vec<(u32, u32)>,
    reverse_dust_masks: Vec<(u32, u32)>,
    ignored_hyphens: Vec<usize>,
}

impl PreparedNucleotideQuery {
    /// Native split-query adapter: retain the clipped inherited full-query
    /// masks while also running DUST on the local sequence, as ObjectMgr does.
    /// Mask coordinates are inclusive and local to this sequence.
    pub fn from_ncbi4na_with_additional_masks(
        id: impl Into<String>,
        sequence: &[u8],
        task: NucleotideTask,
        additional_masks: &[(u32, u32)],
    ) -> Result<Self, QueryError> {
        for &(left, right) in additional_masks {
            if left > right || right as usize >= sequence.len() {
                return Err(QueryError::InvalidMask { left, right });
            }
        }
        let mut query = Self::from_ncbi4na(id, sequence, task)?;
        query.dust_masks =
            s_combine_dust_masks_with_user_provided_masks(&query.dust_masks, additional_masks);
        let last = query.len() as u32 - 1;
        query.reverse_dust_masks = query
            .dust_masks
            .iter()
            .rev()
            .map(|&(start, stop)| (last - stop, last - start))
            .collect();
        Ok(query)
    }
    /// Prepare ASCII IUPAC data. Like the default BLAST query reader, '-' is
    /// omitted and its input position is returned in ignored_hyphens(). This
    /// differs from database FASTA parsing, where a parsed gap becomes N.
    pub fn new(
        id: impl Into<String>,
        sequence: &[u8],
        task: NucleotideTask,
    ) -> Result<Self, QueryError> {
        if sequence.len() > i32::MAX as usize {
            return Err(QueryError::SequenceTooLong);
        }
        let mut encoded = Vec::with_capacity(sequence.len());
        let mut ignored_hyphens = Vec::new();
        for (offset, byte) in sequence.iter().copied().enumerate() {
            if byte == b'-' {
                ignored_hyphens.push(offset);
                continue;
            }
            encoded
                .push(nucleotide_value(byte).ok_or(QueryError::InvalidResidue { offset, byte })?);
        }
        let mut result = Self::from_ncbi4na(id, &encoded, task)?;
        result.ignored_hyphens = ignored_hyphens;
        Ok(result)
    }

    /// Prepare expanded NCBI4na (one value per byte), for example from a database
    /// read through the FASTA adapter. Raw gap code zero becomes N.
    pub fn from_ncbi4na(
        id: impl Into<String>,
        sequence: &[u8],
        task: NucleotideTask,
    ) -> Result<Self, QueryError> {
        let id = id.into();
        if id.trim().is_empty() {
            return Err(QueryError::EmptyIdentifier);
        }
        if sequence.is_empty() {
            return Err(QueryError::EmptySequence);
        }
        if sequence.len() > i32::MAX as usize {
            return Err(QueryError::SequenceTooLong);
        }
        let mut ncbi4na = Vec::with_capacity(sequence.len());
        for (offset, value) in sequence.iter().copied().enumerate() {
            if value > 15 {
                return Err(QueryError::InvalidNcbi4na { offset, value });
            }
            ncbi4na.push(if value == 0 { 15 } else { value });
        }
        let iupac: Vec<_> = ncbi4na
            .iter()
            .map(|base| NCBI4NA_TO_IUPACNA[*base as usize])
            .collect();
        // Upstream Blast_FindDustFilterLoc creates a new duster for each query.
        let mut masker = CSymDustMasker::new(DEFAULT_LEVEL, DEFAULT_WINDOW, DEFAULT_LINKER);
        let dust_masks = masker.mask(&iupac);
        let last = ncbi4na.len() as u32 - 1;
        let reverse_dust_masks = dust_masks
            .iter()
            .rev()
            .map(|(start, stop)| (last - stop, last - start))
            .collect();
        let forward_blastna = ncbi4na
            .iter()
            .map(|base| NCBI4NA_TO_BLASTNA[*base as usize])
            .collect();
        let reversed = get_reverse_nucl_sequence(&ncbi4na, ncbi4na.len() as i32);
        let reverse_blastna = reversed[1..reversed.len() - 1]
            .iter()
            .map(|base| NCBI4NA_TO_BLASTNA[*base as usize])
            .collect();
        Ok(Self {
            id,
            task,
            ncbi4na,
            forward_blastna,
            reverse_blastna,
            dust_masks,
            reverse_dust_masks,
            ignored_hyphens: Vec::new(),
        })
    }

    pub fn id(&self) -> &str {
        &self.id
    }
    pub fn task(&self) -> NucleotideTask {
        self.task
    }
    pub fn len(&self) -> usize {
        self.ncbi4na.len()
    }
    pub fn is_empty(&self) -> bool {
        self.ncbi4na.is_empty()
    }
    pub fn ncbi4na(&self) -> &[u8] {
        &self.ncbi4na
    }
    pub fn forward_blastna(&self) -> &[u8] {
        &self.forward_blastna
    }
    pub fn reverse_blastna(&self) -> &[u8] {
        &self.reverse_blastna
    }
    pub fn dust_masks(&self) -> &[(u32, u32)] {
        &self.dust_masks
    }
    pub fn reverse_dust_masks(&self) -> &[(u32, u32)] {
        &self.reverse_dust_masks
    }
    pub fn soft_masking(&self) -> bool {
        true
    }
    pub fn ignored_hyphens(&self) -> &[usize] {
        &self.ignored_hyphens
    }
}
