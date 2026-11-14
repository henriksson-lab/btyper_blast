// NCBI BLAST+ 2.17.0+, blast_query_info.c/h.
// Author: Christiam Camacho. See LICENSE-NCBI.

use super::blast_program::*;
use super::blast_util::{blast_context_to_frame, blast_get_number_of_contexts};

#[derive(Clone, Debug, Default, Eq, PartialEq)]
#[repr(C)]
pub struct BlastContextInfo {
    pub query_offset: i32,
    pub query_length: i32,
    pub eff_searchsp: i64,
    pub length_adjustment: i32,
    pub query_index: i32,
    pub frame: i8,
    pub is_valid: bool,
    pub segment_flags: i32,
}

/// Query contexts used by nucleotide and translated-subject searches. The
/// PHI-only pattern_info pointer is omitted; no PHI search is exposed.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct BlastQueryInfo {
    pub first_context: i32,
    pub last_context: i32,
    pub num_queries: i32,
    pub contexts: Vec<BlastContextInfo>,
    pub max_length: u32,
    pub min_length: u32,
}

/// Matches C `Blast_GetQueryIndexFromContext`.
pub fn blast_get_query_index_from_context(context: i32, program: EBlastProgramType) -> i32 {
    if program == EBlastProgramType::PsiTblastn || blast_query_is_protein(program) {
        context
    } else if blast_query_is_translated(program) {
        context / 6
    } else {
        context / 2
    }
}

/// Matches C `BlastQueryInfoNew`. Owned storage replaces calloc/free; pattern
/// bookkeeping is outside the supported programs. Nonpositive counts return None.
pub fn blast_query_info_new(
    program: EBlastProgramType,
    num_queries: i32,
) -> Option<BlastQueryInfo> {
    let k_num_contexts = blast_get_number_of_contexts(program);
    if num_queries <= 0 {
        return None;
    }
    assert_ne!(k_num_contexts, 0);
    let mut retval = BlastQueryInfo {
        num_queries,
        first_context: 0,
        last_context: num_queries * k_num_contexts as i32 - 1,
        ..Default::default()
    };
    retval.contexts = vec![BlastContextInfo::default(); retval.last_context as usize + 1];
    for i in 0..=retval.last_context {
        let context = &mut retval.contexts[i as usize];
        context.query_index = blast_get_query_index_from_context(i, program);
        assert_ne!(context.query_index, -1);
        context.frame = blast_context_to_frame(program, i as u32) as i8;
        assert_ne!(context.frame, i8::MAX);
        context.is_valid = true;
        if blast_program_is_mapping(program) {
            context.segment_flags = 0;
        }
    }
    Some(retval)
}

/// Matches C `BSearchContextInfo`, including its min/max-length shortcut and
/// behavior on sentinel offsets between contexts. Valid callers use n >= 0.
pub fn b_search_context_info(n: i32, a: &BlastQueryInfo) -> i32 {
    let size = a.last_context + 1;
    let (mut b, mut e) = if a.min_length > 0 && a.max_length > 0 && a.first_context == 0 {
        // C mixes signed n with Uint4 lengths here. Retain the conversions.
        (
            ((n as u32 / (a.max_length + 1)).min((size - 1) as u32)) as i32,
            ((n as u32 / (a.min_length + 1) + 1).min(size as u32)) as i32,
        )
    } else {
        (0, size)
    };
    assert!(e <= size);
    while b < e - 1 {
        let m = (b + e) / 2;
        if a.contexts[m as usize].query_offset > n {
            e = m;
        } else {
            b = m;
        }
    }
    b
}

/// Matches C `Blast_GetQueryIndexFromQueryOffset`.
pub fn blast_get_query_index_from_query_offset(
    query_offset: i32,
    program: EBlastProgramType,
    query_info: &BlastQueryInfo,
) -> i32 {
    let context = b_search_context_info(query_offset, query_info);
    blast_get_query_index_from_context(context, program)
}

/// Matches C `QueryInfo_GetSeqBufLen`. Includes both endpoint sentinels, with
/// the original special case for a trailing zero-length context.
pub fn query_info_get_seq_buf_len(qinfo: &BlastQueryInfo) -> u32 {
    let cinfo = &qinfo.contexts[qinfo.last_context as usize];
    (cinfo.query_offset + cinfo.query_length + if cinfo.query_length != 0 { 2 } else { 1 }) as u32
}
