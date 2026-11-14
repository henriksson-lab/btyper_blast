//! Protein query ownership for standard tblastn: NCBIstdaa conversion, original
//! SEG masks and the original hard-mask/complement setup. Complete translated
//! search results are assembled by subsequent search stages.
use crate::core::{
    blast_encoding::*,
    blast_filter::*,
    blast_lookup::SSeqRange,
    blast_program::EBlastProgramType,
    blast_query_info::*,
    blast_seg::{seg_parameters_new_aa, seq_buffer_seg},
    blast_setup::s_query_info_set_context,
};
use std::{collections::VecDeque, fmt};

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ProteinQueryError {
    EmptyQuery,
    SequenceTooLong,
    InvalidResidue { offset: usize, byte: u8 },
    InvalidCode { offset: usize, code: u8 },
    Seg(i16),
    EmptyBatch,
    CombinedSequenceTooLong,
}
impl fmt::Display for ProteinQueryError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyQuery => write!(f, "protein query is empty"),
            Self::SequenceTooLong | Self::CombinedSequenceTooLong => {
                write!(f, "protein query exceeds the BLAST offset domain")
            }
            Self::InvalidResidue { offset, byte } => {
                write!(f, "invalid protein residue {byte:#04x} at offset {offset}")
            }
            Self::InvalidCode { offset, code } => {
                write!(f, "invalid NCBIstdaa code {code} at offset {offset}")
            }
            Self::Seg(status) => write!(f, "SEG failed with status {status}"),
            Self::EmptyBatch => write!(f, "protein query batch is empty"),
        }
    }
}
impl std::error::Error for ProteinQueryError {}
#[derive(Debug)]
pub struct PreparedProteinQuery {
    id: String,
    sequence: Vec<u8>,
    masks: Vec<SSeqRange>,
    ignored_hyphens: usize,
}
impl PreparedProteinQuery {
    pub fn new(id: impl Into<String>, sequence: &[u8]) -> Result<Self, ProteinQueryError> {
        let mut encoded = Vec::with_capacity(sequence.len());
        let mut ignored = 0;
        for (offset, &byte) in sequence.iter().enumerate() {
            if byte.is_ascii_whitespace() {
                continue;
            }
            if byte == b'-' {
                ignored += 1;
                continue;
            }
            let code = AMINOACID_TO_NCBISTDAA[byte.to_ascii_uppercase() as usize];
            if code == 0 {
                return Err(ProteinQueryError::InvalidResidue { offset, byte });
            }
            encoded.push(code);
        }
        let mut result = Self::from_ncbistdaa(id, encoded)?;
        result.ignored_hyphens = ignored;
        Ok(result)
    }
    pub fn from_ncbistdaa(
        id: impl Into<String>,
        sequence: impl AsRef<[u8]>,
    ) -> Result<Self, ProteinQueryError> {
        let sequence = sequence.as_ref();
        if sequence.is_empty() {
            return Err(ProteinQueryError::EmptyQuery);
        }
        if sequence.len() > i32::MAX as usize {
            return Err(ProteinQueryError::SequenceTooLong);
        }
        for (offset, &code) in sequence.iter().enumerate() {
            if code >= BLASTAA_SIZE as u8 {
                return Err(ProteinQueryError::InvalidCode { offset, code });
            }
        }
        let mut padded = sequence.to_vec();
        padded.push(K_PROT_SENTINEL);
        let mut masks = VecDeque::new();
        // BlastSetUp_Filter explicitly merges overlapping query SEG segments.
        // SeqBufferSeg's standalone defaults differ from this caller setting.
        let mut seg_parameters = seg_parameters_new_aa();
        seg_parameters.overlaps = true;
        let status = seq_buffer_seg(
            &padded,
            sequence.len() as i32,
            0,
            Some(&mut seg_parameters),
            &mut masks,
        );
        if status != 0 {
            return Err(ProteinQueryError::Seg(status));
        }
        Ok(Self {
            id: id.into(),
            sequence: sequence.to_vec(),
            masks: masks.into(),
            ignored_hyphens: 0,
        })
    }
    pub fn id(&self) -> &str {
        &self.id
    }
    pub fn len(&self) -> usize {
        self.sequence.len()
    }
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }
    pub fn ncbistdaa(&self) -> &[u8] {
        &self.sequence
    }
    pub fn seg_masks(&self) -> &[SSeqRange] {
        &self.masks
    }
    pub fn ignored_hyphens(&self) -> usize {
        self.ignored_hyphens
    }
}

#[derive(Debug)]
pub struct PreparedProteinBatch {
    identifiers: Vec<String>,
    query_info: BlastQueryInfo,
    sequence: MaskedQuerySequence,
    lookup_regions: Vec<SSeqRange>,
}
impl PreparedProteinBatch {
    pub fn new(queries: &[PreparedProteinQuery]) -> Result<Self, ProteinQueryError> {
        if queries.is_empty() {
            return Err(ProteinQueryError::EmptyBatch);
        }
        let total = queries
            .iter()
            .try_fold(0usize, |n, q| n.checked_add(q.len() + 1))
            .ok_or(ProteinQueryError::CombinedSequenceTooLong)?;
        if total > i32::MAX as usize - 2 || queries.len() > i32::MAX as usize {
            return Err(ProteinQueryError::CombinedSequenceTooLong);
        }
        let program = EBlastProgramType::Tblastn;
        let mut info = blast_query_info_new(program, queries.len() as i32).unwrap();
        for (i, q) in queries.iter().enumerate() {
            s_query_info_set_context(&mut info, i as u32, q.len() as u32);
        }
        info.max_length = queries.iter().map(|q| q.len() as u32).max().unwrap();
        info.min_length = queries.iter().map(|q| q.len() as u32).min().unwrap();
        let mut sequence = MaskedQuerySequence {
            sequence_start: vec![K_PROT_SENTINEL; query_info_get_seq_buf_len(&info) as usize + 1],
            sequence_start_nomask: None,
        };
        let mut mask = BlastMaskLoc {
            total_size: queries.len() as i32,
            seqloc_array: queries.iter().map(|q| q.masks.clone()).collect(),
        };
        for (i, q) in queries.iter().enumerate() {
            let offset = info.contexts[i].query_offset as usize + 1;
            sequence.sequence_start[offset..offset + q.len()].copy_from_slice(&q.sequence);
        }
        blast_set_up_mask_query(&mut sequence, &info, &mask, program);
        let mut lookup_regions = vec![];
        blast_complement_mask_locations(program, &info, Some(&mut mask), Some(&mut lookup_regions));
        Ok(Self {
            identifiers: queries.iter().map(|q| q.id.clone()).collect(),
            query_info: info,
            sequence,
            lookup_regions,
        })
    }
    pub fn identifiers(&self) -> &[String] {
        &self.identifiers
    }
    pub fn query_info(&self) -> &BlastQueryInfo {
        &self.query_info
    }
    pub fn sequence_with_sentinels(&self) -> &[u8] {
        &self.sequence.sequence_start[..self.sequence.sequence_start.len() - 1]
    }
    pub fn encoded_contexts(&self) -> &[u8] {
        &self.sequence.sequence_start[1..self.sequence.sequence_start.len() - 2]
    }
    pub fn original_with_sentinels(&self) -> &[u8] {
        self.sequence
            .sequence_start_nomask
            .as_deref()
            .unwrap_or_else(|| self.sequence_with_sentinels())
    }
    pub fn original_contexts(&self) -> &[u8] {
        let raw = self.original_with_sentinels();
        &raw[1..raw.len() - 1]
    }
    pub fn lookup_regions(&self) -> &[SSeqRange] {
        &self.lookup_regions
    }
    pub(crate) fn prepare_query_range(
        &self,
        index: usize,
        range: std::ops::Range<usize>,
    ) -> Result<PreparedProteinQuery, ProteinQueryError> {
        let info = &self.query_info.contexts[index];
        assert!(range.start < range.end && range.end <= info.query_length as usize);
        let start = info.query_offset as usize + range.start;
        let end = info.query_offset as usize + range.end;
        PreparedProteinQuery::from_ncbistdaa(
            &self.identifiers[index],
            &self.original_contexts()[start..end],
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn query_seg_merges_overlaps_like_blast_setup_filter() {
        // Fourth sequence in a fixed LCG stream. Original 2.17.0+
        // BlastSetUp_Filter returns these five merged inclusive intervals;
        // standalone SeqBufferSeg defaults instead return eight overlaps.
        let mut state = 21717u32;
        let mut sequence = vec![0; 320];
        for case in 0..4 {
            for (i, residue) in sequence.iter_mut().enumerate() {
                state = state.wrapping_mul(1664525).wrapping_add(1013904223);
                *residue = if i % 73 < 41 {
                    3 + ((state >> 16) % (3 + case)) as u8
                } else {
                    3 + ((state >> 16) % 18) as u8
                };
            }
        }
        let query = PreparedProteinQuery::from_ncbistdaa("overlapping_seg", sequence).unwrap();
        let masks: Vec<_> = query
            .seg_masks()
            .iter()
            .map(|range| (range.left, range.right))
            .collect();
        assert_eq!(
            masks,
            [(3, 41), (73, 113), (146, 186), (225, 259), (295, 319)]
        );
        let batch = PreparedProteinBatch::new(&[query]).unwrap();
        assert!(batch
            .lookup_regions()
            .iter()
            .all(|range| range.left <= range.right));
    }
}
