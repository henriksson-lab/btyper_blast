// NCBI BLAST+ 2.17.0+, selected blast_diagnostics.c/h definitions.
// Authors: NCBI BLAST developers. See LICENSE-NCBI.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
#[repr(C)]
pub struct BlastUngappedStats {
    pub lookup_hits: i64,
    pub num_seqs_lookup_hits: i32,
    pub init_extends: i32,
    pub good_init_extends: i32,
    pub num_seqs_passed: i32,
}

pub fn blast_ungapped_stats_update(
    stats: Option<&mut BlastUngappedStats>,
    total_hits: i32,
    extended_hits: i32,
    saved_hits: i32,
) {
    let Some(stats) = stats else { return };
    if total_hits == 0 {
        return;
    }
    stats.lookup_hits = stats.lookup_hits.wrapping_add(total_hits as i64);
    stats.num_seqs_lookup_hits = stats.num_seqs_lookup_hits.wrapping_add(1);
    stats.init_extends = stats.init_extends.wrapping_add(extended_hits);
    stats.good_init_extends = stats.good_init_extends.wrapping_add(saved_hits);
    if saved_hits > 0 {
        stats.num_seqs_passed = stats.num_seqs_passed.wrapping_add(1);
    }
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
#[repr(C)]
pub struct BlastGappedStats {
    pub seqs_ungapped_passed: i32,
    pub extensions: i32,
    pub good_extensions: i32,
    pub num_seqs_passed: i32,
}
