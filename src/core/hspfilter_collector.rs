// NCBI BLAST+ 2.17.0+, default hspfilter_collector.c collector.
// Authors: NCBI BLAST developers. See LICENSE-NCBI.
use super::{
    blast_hits::*, blast_options::BlastHitSavingOptions, blast_program::EBlastProgramType,
    blast_query_info::blast_get_query_index_from_context,
};

#[derive(Clone, Debug)]
pub struct BlastHSPCollectorParams {
    pub program: EBlastProgramType,
    pub prelim_hitlist_size: i32,
    pub hsp_num_max: i32,
}

pub fn blast_hsp_collector_params_new(
    hit_options: Option<&BlastHitSavingOptions>,
    composition_based_stats: i32,
    gapped_calculation: bool,
) -> Option<BlastHSPCollectorParams> {
    let options = hit_options?;
    Some(BlastHSPCollectorParams {
        prelim_hitlist_size: get_prelim_hitlist_size(
            options.hitlist_size,
            composition_based_stats,
            gapped_calculation,
        ),
        hsp_num_max: blast_hsp_num_max(gapped_calculation, options.hsp_num_max),
        program: options.program_number,
    })
}

/// Collector state is passed as typed explicit borrows instead of the C
/// writer's void pointers. The incoming list transfers ownership on entry.
pub fn s_blast_hsp_collector_run(
    results: Option<&mut BlastHSPResults>,
    params: Option<&BlastHSPCollectorParams>,
    hsp_list: Option<Box<BlastHSPList>>,
) -> i32 {
    let Some(mut list) = hsp_list else { return 0 };
    let (Some(results), Some(params)) = (results, params) else {
        return -1;
    };
    assert!(matches!(
        params.program,
        EBlastProgramType::Blastn | EBlastProgramType::Tblastn
    ));
    if results.num_queries > 1 {
        let mut lists: Vec<Option<Box<BlastHSPList>>> =
            (0..results.num_queries).map(|_| None).collect();
        for index in 0..list.hspcnt as usize {
            let hsp = list.hsp_array[index].take().expect("nonnull live HSP");
            let query_index =
                blast_get_query_index_from_context(hsp.context, params.program) as usize;
            let target = lists[query_index].get_or_insert_with(|| {
                let mut target = blast_hsp_list_new(params.hsp_num_max);
                target.oid = list.oid;
                target
            });
            blast_hsp_list_save_hsp(target, Some(hsp));
        }
        list.hspcnt = 0;
        drop(list);
        for (index, target) in lists.into_iter().enumerate() {
            if let Some(target) = target {
                let hitlist = results.hitlist_array[index]
                    .get_or_insert_with(|| blast_hit_list_new(params.prelim_hitlist_size));
                blast_hit_list_update(hitlist, target);
            }
        }
    } else if list.hspcnt > 0 {
        let hitlist = results.hitlist_array[0]
            .get_or_insert_with(|| blast_hit_list_new(params.prelim_hitlist_size));
        blast_hit_list_update(hitlist, list);
    }
    0
}
