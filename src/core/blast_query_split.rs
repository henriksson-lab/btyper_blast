// NCBI BLAST+ 2.17.0+, split_query.c and ordinary-query C++ geometry.
// Authors: Christiam Camacho and NCBI BLAST developers. See LICENSE-NCBI.
use super::blast_program::{blast_query_is_translated, EBlastProgramType};
pub const K_BAD_PARAMETER: i16 = -1;
pub const K_OUT_OF_MEMORY: i16 = -2;
pub const K_INVALID_CONTEXT: i32 = -1;
pub type SQueryChunkBoundary = super::blast_lookup::SSeqRange;
#[derive(Clone, Debug)]
pub struct SSplitQueryBlk {
    pub num_chunks: u32,
    pub chunk_query_map: Vec<Vec<u32>>,
    pub chunk_ctx_map: Vec<Vec<i32>>,
    pub chunk_offset_map: Vec<Vec<u32>>,
    pub chunk_bounds: Vec<SQueryChunkBoundary>,
    pub chunk_overlap_sz: usize,
    pub gapped_merge: bool,
}
pub fn split_query_blk_new(num_chunks: u32, gapped_merge: bool) -> Option<Box<SSplitQueryBlk>> {
    if num_chunks == 0 {
        return None;
    }
    Some(Box::new(SSplitQueryBlk {
        num_chunks,
        gapped_merge,
        chunk_query_map: vec![Vec::new(); num_chunks as usize],
        chunk_ctx_map: vec![Vec::new(); num_chunks as usize],
        chunk_offset_map: vec![Vec::new(); num_chunks as usize],
        chunk_bounds: vec![SQueryChunkBoundary { left: 0, right: 0 }; num_chunks as usize],
        chunk_overlap_sz: 0,
    }))
}
pub fn split_query_blk_free(blk: Option<Box<SSplitQueryBlk>>) -> Option<Box<SSplitQueryBlk>> {
    drop(blk);
    None
}
pub fn split_query_blk_allow_gap(blk: Option<&SSplitQueryBlk>) -> bool {
    blk.is_some_and(|b| b.gapped_merge)
}
pub fn split_query_blk_set_chunk_bounds(
    blk: Option<&mut SSplitQueryBlk>,
    chunk_num: u32,
    starting_offset: u32,
    ending_offset: u32,
) -> i16 {
    let Some(b) = blk.filter(|b| chunk_num < b.num_chunks) else {
        return K_BAD_PARAMETER;
    };
    b.chunk_bounds[chunk_num as usize] = SQueryChunkBoundary {
        left: starting_offset as i32,
        right: ending_offset as i32,
    };
    0
}
pub fn split_query_blk_get_chunk_bounds(
    blk: Option<&SSplitQueryBlk>,
    chunk_num: u32,
    starting_offset: Option<&mut usize>,
    ending_offset: Option<&mut usize>,
) -> i16 {
    let (Some(b), Some(start), Some(end)) = (
        blk.filter(|b| chunk_num < b.num_chunks),
        starting_offset,
        ending_offset,
    ) else {
        return K_BAD_PARAMETER;
    };
    *start = b.chunk_bounds[chunk_num as usize].left as usize;
    *end = b.chunk_bounds[chunk_num as usize].right as usize;
    0
}
pub fn split_query_blk_get_num_queries_for_chunk(
    blk: Option<&SSplitQueryBlk>,
    chunk_num: u32,
    num_queries: Option<&mut usize>,
) -> i16 {
    let (Some(b), Some(n)) = (blk.filter(|b| chunk_num < b.num_chunks), num_queries) else {
        return K_BAD_PARAMETER;
    };
    *n = b.chunk_query_map[chunk_num as usize].len();
    0
}
pub fn split_query_blk_add_query_to_chunk(
    blk: Option<&mut SSplitQueryBlk>,
    value: u32,
    chunk_num: u32,
) -> i16 {
    let Some(b) = blk.filter(|b| chunk_num < b.num_chunks) else {
        return K_BAD_PARAMETER;
    };
    let array = &mut b.chunk_query_map[chunk_num as usize];
    if array.try_reserve(1).is_err() {
        return 50; // DynamicInt4/Uint4Array_Append: BLASTERR_MEMORY.
    }
    array.push(value);
    0
}
pub fn split_query_blk_add_context_to_chunk(
    blk: Option<&mut SSplitQueryBlk>,
    value: i32,
    chunk_num: u32,
) -> i16 {
    let Some(b) = blk.filter(|b| chunk_num < b.num_chunks) else {
        return K_BAD_PARAMETER;
    };
    let array = &mut b.chunk_ctx_map[chunk_num as usize];
    if array.try_reserve(1).is_err() {
        return 50; // DynamicInt4/Uint4Array_Append: BLASTERR_MEMORY.
    }
    array.push(value);
    0
}
pub fn split_query_blk_add_context_offset_to_chunk(
    blk: Option<&mut SSplitQueryBlk>,
    value: u32,
    chunk_num: u32,
) -> i16 {
    let Some(b) = blk.filter(|b| chunk_num < b.num_chunks) else {
        return K_BAD_PARAMETER;
    };
    let array = &mut b.chunk_offset_map[chunk_num as usize];
    if array.try_reserve(1).is_err() {
        return 50; // DynamicInt4/Uint4Array_Append: BLASTERR_MEMORY.
    }
    array.push(value);
    0
}
pub fn split_query_blk_get_query_indices_for_chunk(
    blk: Option<&SSplitQueryBlk>,
    chunk_num: u32,
    output: Option<&mut Option<Vec<u32>>>,
) -> i16 {
    let (Some(b), Some(out)) = (blk.filter(|b| chunk_num < b.num_chunks), output) else {
        return K_BAD_PARAMETER;
    };
    *out = None;
    let values = &b.chunk_query_map[chunk_num as usize];
    let mut copy = Vec::new();
    if copy.try_reserve_exact(values.len() + 1).is_err() {
        return K_OUT_OF_MEMORY;
    }
    copy.extend_from_slice(values);
    copy.push(u32::MAX);
    *out = Some(copy);
    0
}
pub fn split_query_blk_get_context_offsets_for_chunk(
    blk: Option<&SSplitQueryBlk>,
    chunk_num: u32,
    output: Option<&mut Option<Vec<u32>>>,
) -> i16 {
    let (Some(b), Some(out)) = (blk.filter(|b| chunk_num < b.num_chunks), output) else {
        return K_BAD_PARAMETER;
    };
    *out = None;
    let values = &b.chunk_offset_map[chunk_num as usize];
    let mut copy = Vec::new();
    if copy.try_reserve_exact(values.len() + 1).is_err() {
        return K_OUT_OF_MEMORY;
    }
    copy.extend_from_slice(values);
    copy.push(u32::MAX);
    *out = Some(copy);
    0
}
pub fn split_query_blk_get_query_contexts_for_chunk(
    blk: Option<&SSplitQueryBlk>,
    chunk_num: u32,
    output: Option<&mut Option<Vec<i32>>>,
    num_query_contexts: Option<&mut u32>,
) -> i16 {
    let (Some(b), Some(out), Some(n)) = (
        blk.filter(|b| chunk_num < b.num_chunks),
        output,
        num_query_contexts,
    ) else {
        return K_BAD_PARAMETER;
    };
    *out = None;
    *n = 0;
    let values = &b.chunk_ctx_map[chunk_num as usize];
    let mut copy = Vec::new();
    if copy.try_reserve_exact(values.len()).is_err() {
        return K_OUT_OF_MEMORY;
    }
    copy.extend_from_slice(values);
    *n = values.len() as u32;
    *out = Some(copy);
    0
}
pub fn split_query_blk_set_chunk_overlap_size(
    blk: Option<&mut SSplitQueryBlk>,
    size: usize,
) -> i16 {
    let Some(b) = blk else { return K_BAD_PARAMETER };
    b.chunk_overlap_sz = size;
    0
}
pub fn split_query_blk_get_chunk_overlap_size(blk: Option<&SSplitQueryBlk>) -> usize {
    blk.map_or(K_BAD_PARAMETER as usize, |b| b.chunk_overlap_sz)
}
/// C++ EProgram subset; core Blastn does not distinguish its task presets.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum QuerySplitProgram {
    Blastn,
    Megablast,
    Tblastn,
}
pub fn split_query_get_chunk_size(program: QuerySplitProgram) -> usize {
    if let Ok(value) = std::env::var("CHUNK_SIZE") {
        if !value.trim().is_empty() {
            return value
                .trim()
                .parse::<i32>()
                .expect("CHUNK_SIZE must be an integer") as usize;
        }
    }
    match program {
        QuerySplitProgram::Blastn => 1_000_000,
        QuerySplitProgram::Megablast => 5_000_000,
        QuerySplitProgram::Tblastn => 20_000,
    }
}
pub fn split_query_get_overlap_chunk_size(program: EBlastProgramType) -> usize {
    if let Ok(value) = std::env::var("OVERLAP_CHUNK_SIZE") {
        if !value.trim().is_empty() {
            return value
                .trim()
                .parse::<i32>()
                .expect("OVERLAP_CHUNK_SIZE must be an integer") as usize;
        }
    }
    if blast_query_is_translated(program) {
        297
    } else {
        100
    }
}
pub fn split_query_should_split(
    program: EBlastProgramType,
    _chunk_size: usize,
    _concatenated_query_length: usize,
    num_queries: usize,
) -> bool {
    !matches!(
        program,
        EBlastProgramType::Mapping
            | EBlastProgramType::RpsBlast
            | EBlastProgramType::RpsTblastn
            | EBlastProgramType::PhiBlastp
            | EBlastProgramType::PhiBlastn
    ) && !(program == EBlastProgramType::Blastx && num_queries > 1)
}
pub fn split_query_calculate_num_chunks(
    program: EBlastProgramType,
    chunk_size: &mut usize,
    concatenated_query_length: usize,
    num_queries: usize,
) -> u32 {
    if !split_query_should_split(program, *chunk_size, concatenated_query_length, num_queries) {
        return 1;
    }
    let overlap = split_query_get_overlap_chunk_size(program);
    if blast_query_is_translated(program) {
        *chunk_size -= *chunk_size % 3
    }
    let num_chunks = if *chunk_size > overlap {
        (concatenated_query_length / (*chunk_size - overlap)) as u32
    } else {
        0
    };
    if num_chunks <= 1 {
        *chunk_size = concatenated_query_length;
        return 1;
    }
    if !blast_query_is_translated(program) {
        *chunk_size =
            (concatenated_query_length + (num_chunks as usize - 1) * overlap) / num_chunks as usize;
        if (num_chunks as usize) < chunk_size.wrapping_sub(overlap) {
            *chunk_size += 1
        }
    }
    num_chunks
}
/// Ordinary-query lowering of CQuerySplitter::x_ComputeChunkRanges.
pub fn x_compute_chunk_ranges(
    blk: &mut SSplitQueryBlk,
    chunk_size: usize,
    total_query_length: usize,
    program: EBlastProgramType,
) {
    let overlap = split_query_get_overlap_chunk_size(program);
    let mut chunk_start = 0usize;
    for chunk_num in 0..blk.num_chunks {
        let mut chunk_end = chunk_start + chunk_size;
        if chunk_end >= total_query_length
            || (chunk_end < total_query_length && chunk_num + 1 == blk.num_chunks)
        {
            chunk_end = total_query_length
        }
        assert!(chunk_end <= u32::MAX as usize && chunk_start <= u32::MAX as usize);
        split_query_blk_set_chunk_bounds(
            Some(blk),
            chunk_num,
            chunk_start as u32,
            chunk_end as u32,
        );
        // In an unsplit query shorter than overlap C's unsigned subtraction
        // wraps. The following break observes the terminal chunk immediately.
        chunk_start = chunk_start.wrapping_add(chunk_size.wrapping_sub(overlap));
        if chunk_start > total_query_length || chunk_end == total_query_length {
            break;
        }
    }
    blk.chunk_overlap_sz = if blast_query_is_translated(program) {
        overlap / 3
    } else {
        overlap
    };
}

/// Comparator from blast_hits.c (distinct from hspstream's reverse-OID comparator).
pub fn s_sort_hsp_list_by_oid(
    x: &super::blast_hits::BlastHSPList,
    y: &super::blast_hits::BlastHSPList,
) -> std::cmp::Ordering {
    x.oid.cmp(&y.oid)
}
pub fn blast_hit_list_merge(
    old_hit_list_ptr: &mut Option<Box<super::blast_hits::BlastHitList>>,
    combined_hit_list_ptr: &mut Option<Box<super::blast_hits::BlastHitList>>,
    contexts_per_query: i32,
    split_offsets: &[i32],
    chunk_overlap_size: i32,
    allow_gap: bool,
) -> i16 {
    use super::blast_hits::{blast_hit_list_new, blast_hit_list_update};
    let Some(mut old) = old_hit_list_ptr.take() else {
        return 0;
    };
    let Some(mut combined) = combined_hit_list_ptr.take() else {
        *combined_hit_list_ptr = Some(old);
        return 0;
    };
    let mut new = blast_hit_list_new(old.hsplist_max);
    let n1 = old.hsplist_count as usize;
    let n2 = combined.hsplist_count as usize;
    if n1 > 1 {
        old.hsplist_array[..n1]
            .sort_by(|x, y| s_sort_hsp_list_by_oid(x.as_deref().unwrap(), y.as_deref().unwrap()));
    }
    if n2 > 1 {
        combined.hsplist_array[..n2]
            .sort_by(|x, y| s_sort_hsp_list_by_oid(x.as_deref().unwrap(), y.as_deref().unwrap()));
    }
    let query_is_split = (0..contexts_per_query).any(|i| split_offsets[i as usize] > 0);
    assert!(chunk_overlap_size != 0);
    let (mut i, mut j) = (0, 0);
    while i < n1 && j < n2 {
        let oid1 = old.hsplist_array[i].as_ref().unwrap().oid;
        let oid2 = combined.hsplist_array[j].as_ref().unwrap().oid;
        if oid1 < oid2 {
            blast_hit_list_update(&mut new, old.hsplist_array[i].take().unwrap());
            i += 1;
        } else if oid1 > oid2 {
            blast_hit_list_update(&mut new, combined.hsplist_array[j].take().unwrap());
            j += 1;
        } else {
            let hsp_max = combined.hsplist_array[j].as_ref().unwrap().hsp_max;
            if query_is_split {
                super::blast_subject_chunk::blast_hsp_lists_merge(
                    &mut old.hsplist_array[i],
                    &mut combined.hsplist_array[j],
                    hsp_max,
                    split_offsets,
                    contexts_per_query,
                    chunk_overlap_size,
                    allow_gap,
                    false,
                );
            } else {
                super::blast_prelim_evaluate::blast_hsp_list_append(
                    &mut old.hsplist_array[i],
                    &mut combined.hsplist_array[j],
                    hsp_max,
                );
            }
            blast_hit_list_update(&mut new, combined.hsplist_array[j].take().unwrap());
            i += 1;
            j += 1;
        }
    }
    for index in i..n1 {
        blast_hit_list_update(&mut new, old.hsplist_array[index].take().unwrap());
    }
    for index in j..n2 {
        blast_hit_list_update(&mut new, combined.hsplist_array[index].take().unwrap());
    }
    old.hsplist_count = 0;
    combined.hsplist_count = 0;
    *combined_hit_list_ptr = Some(new);
    0
}
