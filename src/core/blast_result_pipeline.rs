// NCBI BLAST+ 2.17.0+, selected blast_traceback.c / blast_hspstream.c.
// Authors: NCBI BLAST developers. See LICENSE-NCBI.
use super::blast_hits::{BlastHSPList, BlastHSPResults};
use std::cmp::Ordering;
/// Source comparator orders decreasing OID; stream reads from the end.
pub fn s_sort_hsp_list_by_oid(x: &BlastHSPList, y: &BlastHSPList) -> Ordering {
    y.oid.cmp(&x.oid)
}
/// Source requires results already sorted by best Evalue.
pub fn s_blast_prune_extra_hits(results: &mut BlastHSPResults, hitlist_size: i32) {
    for query_index in 0..results.num_queries as usize {
        let Some(hit_list) = results.hitlist_array[query_index].as_deref_mut() else {
            continue;
        };
        for subject_index in hitlist_size as usize..hit_list.hsplist_count as usize {
            hit_list.hsplist_array[subject_index] = None;
        }
        hit_list.hsplist_count = hit_list.hsplist_count.min(hitlist_size);
    }
}
