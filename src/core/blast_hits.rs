// NCBI BLAST+ 2.17.0+, selected blast_hits.c/h result-storage definitions.
// Authors: NCBI BLAST developers. See LICENSE-NCBI.
use super::blast_stat::BlastScoreBlk;
use super::gapinfo::{EGapAlignOpType, GapEditScript};
use std::cmp::Ordering;

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
#[repr(C)]
pub struct BlastSeg {
    pub frame: i16,
    pub offset: i32,
    pub end: i32,
    pub gapped_start: i32,
}
#[derive(Clone, Copy, Debug, Default)]
pub struct SPHIHspInfo {
    pub index: i32,
    pub length: i32,
}
/// Mapping-only payload is outside the blastn/tblastn storage scope.
#[derive(Debug, Default)]
pub struct BlastHSP {
    pub score: i32,
    pub num_ident: i32,
    pub bit_score: f64,
    pub evalue: f64,
    pub query: BlastSeg,
    pub subject: BlastSeg,
    pub context: i32,
    pub gap_info: Option<Box<GapEditScript>>,
    pub num: i32,
    pub comp_adjustment_method: i16,
    pub pat_info: Option<Box<SPHIHspInfo>>,
    pub num_positives: i32,
}
#[derive(Debug, Default)]
pub struct BlastHSPList {
    pub oid: i32,
    pub query_index: i32,
    pub hsp_array: Vec<Option<Box<BlastHSP>>>,
    pub hspcnt: i32,
    pub allocated: i32,
    pub hsp_max: i32,
    pub do_not_reallocate: bool,
    pub best_evalue: f64,
}
#[derive(Debug, Default)]
pub struct BlastHitList {
    pub hsplist_count: i32,
    pub hsplist_max: i32,
    pub worst_evalue: f64,
    pub low_score: i32,
    pub heapified: bool,
    pub hsplist_array: Vec<Option<Box<BlastHSPList>>>,
    pub hsplist_current: i32,
    pub num_hits: i32,
}
#[derive(Debug, Default)]
pub struct BlastHSPResults {
    pub num_queries: i32,
    pub hitlist_array: Vec<Option<Box<BlastHitList>>>,
}

pub fn blast_hsp_new() -> Box<BlastHSP> {
    Box::new(BlastHSP::default())
}
pub fn blast_hsp_init(
    query_start: i32,
    query_end: i32,
    subject_start: i32,
    subject_end: i32,
    query_gapped_start: i32,
    subject_gapped_start: i32,
    query_context: i32,
    query_frame: i16,
    subject_frame: i16,
    score: i32,
    gap_edit: Option<&mut Option<GapEditScript>>,
    ret_hsp: Option<&mut Option<Box<BlastHSP>>>,
) -> i16 {
    let Some(ret_hsp) = ret_hsp else {
        return -1;
    };
    let mut new_hsp = blast_hsp_new();
    *ret_hsp = None;
    new_hsp.query.offset = query_start;
    new_hsp.subject.offset = subject_start;
    new_hsp.query.end = query_end;
    new_hsp.subject.end = subject_end;
    new_hsp.query.gapped_start = query_gapped_start;
    new_hsp.subject.gapped_start = subject_gapped_start;
    new_hsp.context = query_context;
    new_hsp.query.frame = query_frame;
    new_hsp.subject.frame = subject_frame;
    new_hsp.score = score;
    if let Some(edit) = gap_edit {
        if let Some(edit) = edit.take() {
            new_hsp.gap_info = Some(Box::new(edit));
        }
    }
    *ret_hsp = Some(new_hsp);
    0
}
pub fn blast_hsp_calc_length_and_gaps(
    hsp: &BlastHSP,
    length_out: &mut i32,
    gaps_out: &mut i32,
    gap_opens_out: &mut i32,
) {
    let mut length = hsp.query.end - hsp.query.offset;
    let s_length = hsp.subject.end - hsp.subject.offset;
    let (mut gap_opens, mut gaps) = (0, 0);
    if let Some(esp) = &hsp.gap_info {
        for index in 0..esp.size as usize {
            if esp.op_type[index] == EGapAlignOpType::Del {
                length += esp.num[index];
                gaps += esp.num[index];
                gap_opens += 1;
            } else if esp.op_type[index] == EGapAlignOpType::Ins {
                gap_opens += 1;
                gaps += esp.num[index];
            }
        }
    } else if s_length > length {
        length = s_length;
    }
    *length_out = length;
    *gap_opens_out = gap_opens;
    *gaps_out = gaps;
}
pub fn score_compare_hs_ps(hsp1: Option<&BlastHSP>, hsp2: Option<&BlastHSP>) -> Ordering {
    let (hsp1, hsp2) = match (hsp1, hsp2) {
        (None, None) => return Ordering::Equal,
        (None, Some(_)) => return Ordering::Greater,
        (Some(_), None) => return Ordering::Less,
        (Some(a), Some(b)) => (a, b),
    };
    hsp2.score
        .cmp(&hsp1.score)
        .then(hsp1.subject.offset.cmp(&hsp2.subject.offset))
        .then(hsp2.subject.end.cmp(&hsp1.subject.end))
        .then(hsp1.query.offset.cmp(&hsp2.query.offset))
        .then(hsp2.query.end.cmp(&hsp1.query.end))
}
pub fn blast_hsp_list_is_sorted_by_score(hsp_list: Option<&BlastHSPList>) -> bool {
    let Some(list) = hsp_list else {
        return true;
    };
    if list.hspcnt <= 1 {
        return true;
    }
    for pair in list.hsp_array[..list.hspcnt as usize].windows(2) {
        if score_compare_hs_ps(pair[0].as_deref(), pair[1].as_deref()) == Ordering::Greater {
            return false;
        }
    }
    true
}
pub fn blast_hsp_list_sort_by_score(hsp_list: Option<&mut BlastHSPList>) {
    let Some(list) = hsp_list else {
        return;
    };
    if list.hspcnt <= 1 {
        return;
    }
    if !blast_hsp_list_is_sorted_by_score(Some(list)) {
        list.hsp_array[..list.hspcnt as usize]
            .sort_by(|a, b| score_compare_hs_ps(a.as_deref(), b.as_deref()));
    }
}
pub fn s_evalue_comp(evalue1: f64, evalue2: f64) -> Ordering {
    let epsilon = 1.0e-180;
    if evalue1 < epsilon && evalue2 < epsilon {
        return Ordering::Equal;
    }
    if evalue1 < evalue2 {
        Ordering::Less
    } else if evalue1 > evalue2 {
        Ordering::Greater
    } else {
        Ordering::Equal
    }
}
pub fn s_evalue_compare_hs_ps(h1: Option<&BlastHSP>, h2: Option<&BlastHSP>) -> Ordering {
    let (h1, h2) = match (h1, h2) {
        (None, None) => return Ordering::Equal,
        (None, Some(_)) => return Ordering::Greater,
        (Some(_), None) => return Ordering::Less,
        (Some(a), Some(b)) => (a, b),
    };
    let result = s_evalue_comp(h1.evalue, h2.evalue);
    if result != Ordering::Equal {
        return result;
    }
    score_compare_hs_ps(Some(h1), Some(h2))
}
pub fn blast_hsp_list_sort_by_evalue(hsp_list: Option<&mut BlastHSPList>) {
    let Some(list) = hsp_list else {
        return;
    };
    if list.hspcnt > 1 {
        let mut index = 0;
        while index < list.hspcnt - 1 {
            if s_evalue_compare_hs_ps(
                list.hsp_array[index as usize].as_deref(),
                list.hsp_array[index as usize + 1].as_deref(),
            ) == Ordering::Greater
            {
                break;
            }
            index += 1;
        }
        if index < list.hspcnt - 1 {
            list.hsp_array[..list.hspcnt as usize]
                .sort_by(|a, b| s_evalue_compare_hs_ps(a.as_deref(), b.as_deref()));
        }
    }
}
pub fn blast_hsp_list_new(hsp_max: i32) -> Box<BlastHSPList> {
    let mut list = Box::new(BlastHSPList::default());
    list.hsp_max = i32::MAX;
    if hsp_max > 0 {
        list.hsp_max = hsp_max;
    }
    list.allocated = 100.min(list.hsp_max);
    list.hsp_array = (0..list.allocated).map(|_| None).collect();
    list
}
pub fn blast_hsp_list_is_empty(hsp_list: Option<&BlastHSPList>) -> bool {
    hsp_list.is_some_and(|list| list.hspcnt == 0)
}
pub fn blast_hsp_list_swap(list1: &mut BlastHSPList, list2: &mut BlastHSPList) {
    std::mem::swap(list1, list2);
}

/// Typed bindings for the pointer arrays used by the original byte heap.
/// Offsets are measured in pointer elements instead of bytes.
pub enum BlastHeapView<'a> {
    Hsps(&'a mut [Option<Box<BlastHSP>>]),
    HspLists(&'a mut [Option<Box<BlastHSPList>>]),
}
#[derive(Clone, Copy)]
pub enum BlastHeapComparator {
    HspScore,
    HspListEvalue,
}
impl BlastHeapView<'_> {
    pub(crate) fn compare(&self, a: usize, b: usize, compar: BlastHeapComparator) -> Ordering {
        match (self, compar) {
            (Self::Hsps(items), BlastHeapComparator::HspScore) => {
                score_compare_hs_ps(items[a].as_deref(), items[b].as_deref())
            }
            (Self::HspLists(items), BlastHeapComparator::HspListEvalue) => {
                s_evalue_compare_hsp_lists(
                    items[a].as_deref().unwrap(),
                    items[b].as_deref().unwrap(),
                )
            }
            _ => panic!("heap array and comparator must agree"),
        }
    }
    pub(crate) fn swap(&mut self, a: usize, b: usize) {
        match self {
            Self::Hsps(items) => items.swap(a, b),
            Self::HspLists(items) => items.swap(a, b),
        }
    }
}
pub fn s_heapify(
    heap: &mut BlastHeapView<'_>,
    base0: usize,
    mut base: usize,
    lim: usize,
    last: usize,
    compar: BlastHeapComparator,
) {
    let mut left_son = base0 + 2 * (base - base0) + 1;
    while base <= lim {
        let large_son =
            if left_son == last || heap.compare(left_son, left_son + 1, compar) != Ordering::Less {
                left_son
            } else {
                left_son + 1
            };
        if heap.compare(base, large_son, compar) == Ordering::Less {
            heap.swap(base, large_son);
            base = large_son;
            left_son = base0 + 2 * (base - base0) + 1;
        } else {
            break;
        }
    }
}
pub fn s_create_heap(heap: &mut BlastHeapView<'_>, nel: usize, compar: BlastHeapComparator) {
    if nel < 2 {
        return;
    }
    let lim = (nel - 2) / 2;
    let last = nel - 1;
    let mut i = nel / 2;
    let mut base = i - 1;
    while i > 0 {
        s_heapify(heap, 0, base, lim, last, compar);
        i -= 1;
        if i > 0 {
            base -= 1;
        }
    }
}
pub fn s_blast_hsp_list_insert_hsp_in_heap(list: &mut BlastHSPList, hsp: Option<Box<BlastHSP>>) {
    if score_compare_hs_ps(hsp.as_deref(), list.hsp_array[0].as_deref()) == Ordering::Greater {
        return;
    }
    list.hsp_array[0] = hsp;
    if list.hspcnt >= 2 {
        s_heapify(
            &mut BlastHeapView::Hsps(&mut list.hsp_array),
            0,
            0,
            list.hspcnt as usize / 2 - 1,
            list.hspcnt as usize - 1,
            BlastHeapComparator::HspScore,
        );
    }
}
pub fn s_blast_check_best_evalue(list: &BlastHSPList) -> bool {
    if list.hspcnt == 0 {
        return true;
    }
    let mut best_evalue = i32::MAX as f64;
    for item in &list.hsp_array[..list.hspcnt as usize] {
        let value = item.as_ref().unwrap().evalue;
        if value < best_evalue {
            best_evalue = value;
        }
    }
    if (best_evalue - list.best_evalue).abs() / (best_evalue + 1.0e-200) > 0.01 {
        return false;
    }
    true
}
pub fn s_blast_get_best_evalue(list: &BlastHSPList) -> f64 {
    let mut best_evalue = i32::MAX as f64;
    for item in &list.hsp_array[..list.hspcnt as usize] {
        let value = item.as_ref().unwrap().evalue;
        if value < best_evalue {
            best_evalue = value;
        }
    }
    best_evalue
}
pub fn blast_hsp_list_save_hsp(list: &mut BlastHSPList, new_hsp: Option<Box<BlastHSP>>) -> i16 {
    let hspcnt = list.hspcnt;
    let mut hsp_allocated = list.allocated;
    let mut status = 0;
    if hspcnt >= hsp_allocated && !list.do_not_reallocate {
        let new_allocated = list.allocated.wrapping_mul(2).min(list.hsp_max);
        if new_allocated > list.allocated {
            if list
                .hsp_array
                .try_reserve_exact(new_allocated as usize - list.hsp_array.len())
                .is_err()
            {
                list.do_not_reallocate = true;
                status = -1;
            } else {
                list.hsp_array.resize_with(new_allocated as usize, || None);
                list.allocated = new_allocated;
                hsp_allocated = new_allocated;
            }
        } else {
            list.do_not_reallocate = true;
        }
        if list.do_not_reallocate {
            s_create_heap(
                &mut BlastHeapView::Hsps(&mut list.hsp_array),
                hspcnt as usize,
                BlastHeapComparator::HspScore,
            );
        }
    }
    if hspcnt < hsp_allocated {
        list.hsp_array[list.hspcnt as usize] = new_hsp;
        list.hspcnt += 1;
        return status;
    }
    s_blast_hsp_list_insert_hsp_in_heap(list, new_hsp);
    status
}
pub fn blast_hsp_list_purge_null_hs_ps(hsp_list: Option<&mut BlastHSPList>) -> i16 {
    let Some(list) = hsp_list else {
        return 0;
    };
    if list.hspcnt == 0 {
        return 0;
    }
    let mut index1 = 0;
    let hspcnt = list.hspcnt as usize;
    for index in 0..hspcnt {
        if list.hsp_array[index].is_some() {
            if index != index1 {
                list.hsp_array[index1] = list.hsp_array[index].take();
            }
            index1 += 1;
        }
    }
    for index in index1..hspcnt {
        list.hsp_array[index] = None;
    }
    list.hspcnt = index1 as i32;
    0
}
pub fn blast_hsp_list_reap_by_evalue(hsp_list: Option<&mut BlastHSPList>, cutoff: f64) -> i16 {
    let Some(list) = hsp_list else {
        return 0;
    };
    let mut hsp_cnt = 0;
    for index in 0..list.hspcnt as usize {
        if list.hsp_array[index].as_ref().unwrap().evalue > cutoff {
            list.hsp_array[index] = None;
        } else {
            if index > hsp_cnt {
                list.hsp_array[hsp_cnt] = list.hsp_array[index].take();
            }
            hsp_cnt += 1;
        }
    }
    list.hspcnt = hsp_cnt as i32;
    0
}
pub fn blast_hsp_list_reap_by_raw_score(
    hsp_list: Option<&mut BlastHSPList>,
    cutoff_score: i32,
) -> i16 {
    let Some(list) = hsp_list else {
        return 0;
    };
    let mut hsp_cnt = 0;
    for index in 0..list.hspcnt as usize {
        if list.hsp_array[index].as_ref().unwrap().score < cutoff_score {
            list.hsp_array[index] = None;
        } else {
            if index > hsp_cnt {
                list.hsp_array[hsp_cnt] = list.hsp_array[index].take();
            }
            hsp_cnt += 1;
        }
    }
    list.hspcnt = hsp_cnt as i32;
    0
}
pub fn blast_trim_hsp_list_by_max_hsps(
    hsp_list: Option<&mut BlastHSPList>,
    max_hsps_per_subject: i32,
) -> i16 {
    let Some(list) = hsp_list else {
        return 0;
    };
    if max_hsps_per_subject == 0 || list.hspcnt <= max_hsps_per_subject {
        return 0;
    }
    assert!(max_hsps_per_subject >= 0);
    for index in max_hsps_per_subject..list.hspcnt {
        list.hsp_array[index as usize] = None;
    }
    list.hspcnt = max_hsps_per_subject;
    0
}
pub fn blast_hsp_list_adjust_offsets(list: &mut BlastHSPList, offset: i32) {
    if offset == 0 {
        return;
    }
    for hsp in &mut list.hsp_array[..list.hspcnt as usize] {
        let hsp = hsp.as_mut().unwrap();
        hsp.subject.offset = hsp.subject.offset.wrapping_add(offset);
        hsp.subject.end = hsp.subject.end.wrapping_add(offset);
        hsp.subject.gapped_start = hsp.subject.gapped_start.wrapping_add(offset);
    }
}
pub fn blast_hsp_list_adjust_odd_blastn_scores(
    hsp_list: Option<&mut BlastHSPList>,
    gapped_calculation: bool,
    round_down: bool,
) {
    let Some(list) = hsp_list else {
        return;
    };
    if list.hspcnt == 0 || !gapped_calculation || !round_down {
        return;
    }
    for hsp in &mut list.hsp_array[..list.hspcnt as usize] {
        hsp.as_mut().unwrap().score &= !1;
    }
    blast_hsp_list_sort_by_score(Some(list));
}
pub fn blast_hsp_list_get_bit_scores(
    hsp_list: Option<&mut BlastHSPList>,
    gapped_calculation: bool,
    sbp: &BlastScoreBlk,
) -> i16 {
    let Some(list) = hsp_list else {
        return 1;
    };
    let kbp = if gapped_calculation {
        sbp.kbp_gap.as_deref().expect("gapped blocks required")
    } else {
        sbp.ungapped_kbp()
    };
    for item in &mut list.hsp_array[..list.hspcnt as usize] {
        let hsp = item.as_mut().unwrap();
        let k = kbp[hsp.context as usize].as_ref().unwrap();
        hsp.bit_score = (hsp.score as f64 * k.lambda - k.log_k) / std::f64::consts::LN_2;
    }
    0
}
pub fn s_evalue_compare_hsp_lists(h1: &BlastHSPList, h2: &BlastHSPList) -> Ordering {
    if h1.hspcnt == 0 && h2.hspcnt == 0 {
        return Ordering::Equal;
    } else if h1.hspcnt == 0 {
        return Ordering::Greater;
    } else if h2.hspcnt == 0 {
        return Ordering::Less;
    }
    let retval = s_evalue_comp(h1.best_evalue, h2.best_evalue);
    if retval != Ordering::Equal {
        return retval;
    }
    let score1 = h1.hsp_array[0].as_ref().unwrap().score;
    let score2 = h2.hsp_array[0].as_ref().unwrap().score;
    if score1 > score2 {
        return Ordering::Less;
    }
    if score1 < score2 {
        return Ordering::Greater;
    }
    h2.oid.cmp(&h1.oid)
}
pub fn s_evalue_compare_hsp_lists_rev(h1: &BlastHSPList, h2: &BlastHSPList) -> Ordering {
    s_evalue_compare_hsp_lists(h2, h1)
}
pub fn blast_hit_list_new(hitlist_size: i32) -> Box<BlastHitList> {
    Box::new(BlastHitList {
        hsplist_max: hitlist_size,
        low_score: i32::MAX,
        ..Default::default()
    })
}
pub fn s_blast_hit_list_purge(hit_list: Option<&mut BlastHitList>) {
    let Some(list) = hit_list else {
        return;
    };
    let hsplist_count = list.hsplist_count as usize;
    let mut index = 0;
    while index < hsplist_count && list.hsplist_array[index].as_ref().unwrap().hspcnt > 0 {
        index += 1;
    }
    list.hsplist_count = index as i32;
    while index < hsplist_count {
        list.hsplist_array[index] = None;
        index += 1;
    }
}
pub fn s_blast_hit_list_insert_hsp_list_in_heap(
    list: &mut BlastHitList,
    hsp_list: Box<BlastHSPList>,
) {
    list.hsplist_array[0] = Some(hsp_list);
    if list.hsplist_count >= 2 {
        s_heapify(
            &mut BlastHeapView::HspLists(&mut list.hsplist_array),
            0,
            0,
            list.hsplist_count as usize / 2 - 1,
            list.hsplist_count as usize - 1,
            BlastHeapComparator::HspListEvalue,
        );
    }
    let worst = list.hsplist_array[0].as_ref().unwrap();
    list.worst_evalue = worst.best_evalue;
    list.low_score = worst.hsp_array[0].as_ref().unwrap().score;
}
pub fn s_blast_hit_list_grow_hsp_list_array(list: &mut BlastHitList) -> i16 {
    if list.hsplist_current >= list.hsplist_max {
        return 1;
    }
    if list.hsplist_current <= 0 {
        list.hsplist_current = 100;
    } else {
        list.hsplist_current = list.hsplist_current.wrapping_mul(2).min(list.hsplist_max);
    }
    if list
        .hsplist_array
        .try_reserve_exact(list.hsplist_current as usize - list.hsplist_array.len())
        .is_err()
    {
        return 50;
    }
    list.hsplist_array
        .resize_with(list.hsplist_current as usize, || None);
    0
}
pub fn blast_hit_list_update(list: &mut BlastHitList, mut hsp_list: Box<BlastHSPList>) -> i16 {
    hsp_list.best_evalue = s_blast_get_best_evalue(&hsp_list);
    debug_assert!(s_blast_check_best_evalue(&hsp_list));
    if list.hsplist_count < list.hsplist_max {
        if list.hsplist_current == list.hsplist_count {
            let status = s_blast_hit_list_grow_hsp_list_array(list);
            if status != 0 {
                return status;
            }
        }
        let best = hsp_list.best_evalue;
        let score = hsp_list.hsp_array[0].as_ref().unwrap().score;
        list.hsplist_array[list.hsplist_count as usize] = Some(hsp_list);
        list.hsplist_count += 1;
        if best > list.worst_evalue {
            list.worst_evalue = best;
        }
        if score < list.low_score {
            list.low_score = score;
        }
    } else {
        if !list.heapified {
            for item in &mut list.hsplist_array[..list.hsplist_count as usize] {
                let item = item.as_mut().unwrap();
                blast_hsp_list_sort_by_evalue(Some(item));
                item.best_evalue = s_blast_get_best_evalue(item);
            }
            s_create_heap(
                &mut BlastHeapView::HspLists(&mut list.hsplist_array),
                list.hsplist_count as usize,
                BlastHeapComparator::HspListEvalue,
            );
            list.heapified = true;
        }
        blast_hsp_list_sort_by_evalue(Some(&mut hsp_list));
        hsp_list.best_evalue = s_blast_get_best_evalue(&hsp_list);
        let order = s_evalue_compare_hsp_lists(list.hsplist_array[0].as_ref().unwrap(), &hsp_list);
        if order != Ordering::Less {
            s_blast_hit_list_insert_hsp_list_in_heap(list, hsp_list);
        }
    }
    0
}
pub fn blast_hit_list_purge_null_hsp_lists(hit_list: Option<&mut BlastHitList>) -> i16 {
    let Some(list) = hit_list else {
        return 0;
    };
    if list.hsplist_count == 0 {
        return 0;
    }
    let mut index1 = 0;
    let count = list.hsplist_count as usize;
    for index in 0..count {
        if list.hsplist_array[index].is_some() {
            if index != index1 {
                list.hsplist_array[index1] = list.hsplist_array[index].take();
            }
            index1 += 1;
        }
    }
    for index in index1..count {
        list.hsplist_array[index] = None;
    }
    list.hsplist_count = index1 as i32;
    0
}
pub fn blast_hit_list_sort_by_evalue(hit_list: Option<&mut BlastHitList>) -> i16 {
    let Some(list) = hit_list else {
        return 0;
    };
    if list.hsplist_count > 1 {
        list.hsplist_array[..list.hsplist_count as usize]
            .sort_by(|a, b| s_evalue_compare_hsp_lists(a.as_ref().unwrap(), b.as_ref().unwrap()));
    }
    s_blast_hit_list_purge(Some(list));
    0
}
pub fn blast_hsp_results_new(num_queries: i32) -> Box<BlastHSPResults> {
    assert!(num_queries >= 0);
    Box::new(BlastHSPResults {
        num_queries,
        hitlist_array: (0..num_queries).map(|_| None).collect(),
    })
}
pub fn blast_hsp_results_sort_by_evalue(results: Option<&mut BlastHSPResults>) -> i16 {
    let Some(results) = results else {
        return 0;
    };
    for item in &mut results.hitlist_array[..results.num_queries as usize] {
        if let Some(list) = item {
            if list.hsplist_count > 1 && !list.hsplist_array.is_empty() {
                list.hsplist_array[..list.hsplist_count as usize].sort_by(|a, b| {
                    s_evalue_compare_hsp_lists(a.as_ref().unwrap(), b.as_ref().unwrap())
                });
            }
            s_blast_hit_list_purge(Some(list));
        }
    }
    0
}
pub fn blast_hsp_results_reverse_sort(results: &mut BlastHSPResults) -> i16 {
    for list in results.hitlist_array[..results.num_queries as usize]
        .iter_mut()
        .flatten()
    {
        if list.hsplist_count > 1 {
            list.hsplist_array[..list.hsplist_count as usize].sort_by(|a, b| {
                s_evalue_compare_hsp_lists_rev(a.as_ref().unwrap(), b.as_ref().unwrap())
            });
        }
        s_blast_hit_list_purge(Some(list));
    }
    0
}
pub fn blast_hsp_results_reverse_order(results: &mut BlastHSPResults) -> i16 {
    for list in results.hitlist_array[..results.num_queries as usize]
        .iter_mut()
        .flatten()
    {
        if list.hsplist_count > 1 {
            for i in 0..list.hsplist_count as usize / 2 {
                list.hsplist_array
                    .swap(i, list.hsplist_count as usize - i - 1);
            }
        }
    }
    0
}
pub fn blast_hsp_results_insert_hsp_list(
    results: &mut BlastHSPResults,
    hsp_list: Option<Box<BlastHSPList>>,
    hitlist_size: i32,
) -> i16 {
    let Some(hsp_list) = hsp_list else {
        return 0;
    };
    if hsp_list.hspcnt == 0 {
        return 0;
    }
    let query = hsp_list.query_index;
    assert!(query < results.num_queries);
    let target = &mut results.hitlist_array[query as usize];
    if target.is_none() {
        *target = Some(blast_hit_list_new(hitlist_size));
    }
    blast_hit_list_update(target.as_deref_mut().unwrap(), hsp_list);
    0
}

pub fn get_prelim_hitlist_size(
    hitlist_size: i32,
    composition_based_stats: i32,
    gapped_calculation: bool,
) -> i32 {
    let mut prelim_hitlist_size = hitlist_size;
    let adaptive = std::env::var_os("ADAPTIVE_CBS").is_some();
    if composition_based_stats != 0 {
        if adaptive {
            if hitlist_size < 1000 {
                prelim_hitlist_size = prelim_hitlist_size.wrapping_add(1000).max(1500);
            } else {
                prelim_hitlist_size = prelim_hitlist_size.wrapping_mul(2).wrapping_add(50);
            }
        } else if hitlist_size <= 500 {
            prelim_hitlist_size = 1050;
        } else {
            prelim_hitlist_size = prelim_hitlist_size.wrapping_mul(2).wrapping_add(50);
        }
    } else if gapped_calculation {
        prelim_hitlist_size = prelim_hitlist_size
            .wrapping_mul(2)
            .max(10)
            .min(prelim_hitlist_size.wrapping_add(50));
    }
    prelim_hitlist_size
}
pub fn blast_hsp_num_max(_gapped_calculation: bool, hsp_num_max: i32) -> i32 {
    if hsp_num_max <= 0 {
        i32::MAX
    } else {
        hsp_num_max
    }
}
