// NCBI BLAST+ 2.17.0+, selected blast_filter.c/h definitions.
// Authors: NCBI BLAST developers. See LICENSE-NCBI.
use super::blast_lookup::SSeqRange;
use super::blast_program::EBlastProgramType;
use super::blast_query_info::BlastQueryInfo;

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct BlastMaskLoc {
    pub total_size: i32,
    pub seqloc_array: Vec<Vec<SSeqRange>>,
}

/// Matches C `BlastIsReverseStrand`.
#[inline]
pub fn blast_is_reverse_strand(is_na: bool, context: i32) -> bool {
    is_na && context & 1 != 0
}

/// Matches C `BlastSeqLocListReverse`. Ordered slices replace linked nodes
/// and the temporary array of node pointers, preserving exact element order.
pub fn blast_seq_loc_list_reverse(head: Option<&mut [SSeqRange]>) {
    let Some(head) = head else {
        return;
    };
    head.reverse();
}

/// Matches C `BLAST_ComplementMaskLocations`. Reverse-context mask lists are
/// reversed in place, exactly as through the C function's const container.
/// Masks for nucleotide contexts enter in forward-query coordinates.
pub fn blast_complement_mask_locations(
    program_number: EBlastProgramType,
    query_info: &BlastQueryInfo,
    mut mask_loc: Option<&mut BlastMaskLoc>,
    complement_mask: Option<&mut Vec<SSeqRange>>,
) -> i16 {
    let k_is_nucl = matches!(
        program_number,
        EBlastProgramType::Blastn | EBlastProgramType::Mapping
    );
    let Some(complement_mask) = complement_mask else {
        return -1;
    };
    complement_mask.clear();
    for context in query_info.first_context..=query_info.last_context {
        let mut first = true;
        let mut last_interval_open = true;
        let mut left = 0;
        let info = &query_info.contexts[context as usize];
        if !info.is_valid {
            continue;
        }
        let start_offset = info.query_offset;
        let end_offset = info.query_length + start_offset - 1;
        assert!(start_offset <= end_offset);
        if mask_loc
            .as_ref()
            .is_none_or(|mask| mask.seqloc_array[context as usize].is_empty())
        {
            complement_mask.push(SSeqRange {
                left: start_offset,
                right: end_offset,
            });
            continue;
        }
        let locations = &mut mask_loc.as_mut().unwrap().seqloc_array[context as usize];
        if blast_is_reverse_strand(k_is_nucl, context) {
            blast_seq_loc_list_reverse(Some(locations));
        }
        for seq_range in locations {
            let (filter_start, filter_end) = if blast_is_reverse_strand(k_is_nucl, context) {
                (end_offset - seq_range.right, end_offset - seq_range.left)
            } else {
                (
                    start_offset + seq_range.left,
                    start_offset + seq_range.right,
                )
            };
            if first {
                last_interval_open = true;
                first = false;
                if filter_start > start_offset {
                    left = start_offset;
                } else {
                    left = filter_end + 1;
                    continue;
                }
            }
            let right = filter_start - 1;
            complement_mask.push(SSeqRange { left, right });
            if filter_end >= end_offset {
                last_interval_open = false;
                break;
            } else {
                left = filter_end + 1;
            }
        }
        if last_interval_open {
            complement_mask.push(SSeqRange {
                left,
                right: end_offset,
            });
        }
    }
    0
}

pub const K_NUCL_MASK: u8 = 14;
pub const K_PROT_MASK: u8 = 21;
/// Owned source sequence-block fields needed when hard masking query residues.
#[derive(Debug)]
pub struct MaskedQuerySequence {
    pub sequence_start: Vec<u8>,
    pub sequence_start_nomask: Option<Vec<u8>>,
}

pub fn blast_mask_the_residues(
    buffer: &mut [u8],
    length: i32,
    is_na: bool,
    mask_loc: &[SSeqRange],
    reverse: bool,
    offset: i32,
) {
    let letter = if is_na { K_NUCL_MASK } else { K_PROT_MASK };
    for mask in mask_loc {
        let (mut start, mut stop) = if reverse {
            (length - 1 - mask.right, length - 1 - mask.left)
        } else {
            (mask.left, mask.right)
        };
        start -= offset;
        stop -= offset;
        assert!(start < length);
        assert!(stop <= length);
        for index in start..=stop {
            buffer[index as usize] = letter;
        }
    }
}
pub fn blast_set_up_mask_query(
    query: &mut MaskedQuerySequence,
    query_info: &BlastQueryInfo,
    filter_maskloc: &BlastMaskLoc,
    program: EBlastProgramType,
) {
    let is_nucl = matches!(
        program,
        EBlastProgramType::Blastn | EBlastProgramType::Mapping
    );
    let mut has_mask = false;
    for index in 0..filter_maskloc.total_size as usize {
        if !filter_maskloc.seqloc_array[index].is_empty() {
            has_mask = true;
            break;
        }
    }
    if !has_mask {
        return;
    }
    let last = &query_info.contexts[query_info.last_context as usize];
    let total_length = (last.query_offset + last.query_length + 2) as usize;
    query.sequence_start_nomask = Some(query.sequence_start[..total_length].to_vec());
    for context in query_info.first_context..=query_info.last_context {
        let ctx = &query_info.contexts[context as usize];
        if !ctx.is_valid {
            continue;
        }
        blast_mask_the_residues(
            &mut query.sequence_start[ctx.query_offset as usize + 1..],
            ctx.query_length,
            is_nucl,
            &filter_maskloc.seqloc_array[context as usize],
            blast_is_reverse_strand(is_nucl, context),
            0,
        );
    }
}
