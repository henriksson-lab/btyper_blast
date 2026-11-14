// NCBI BLAST+ 2.17.0+, default unmasked subject splitting and preliminary merge.
// Authors: NCBI BLAST developers. See LICENSE-NCBI.
use super::blast_gapalign::MAX_DBSEQ_LEN;
use super::blast_hits::{blast_hsp_list_purge_null_hs_ps, BlastHSP, BlastHSPList};
use super::blast_lookup::SSeqRange;
use super::blast_prelim_evaluate::s_blast_hsp_lists_combine_by_score;
pub const DBSEQ_CHUNK_OVERLAP: i32 = 100;
pub const OVERLAP_DIAG_CLOSE: i32 = 10;
pub const SUBJECT_SPLIT_DONE: i16 = 0;
pub const SUBJECT_SPLIT_OK: i16 = 1;
#[derive(Clone, Copy, Debug)]
pub struct SubjectSplitStruct {
    pub initialized: bool,
    pub full_range: SSeqRange,
    pub offset: i32,
    pub next: i32,
}
#[derive(Clone, Copy, Debug)]
pub struct SubjectChunk {
    pub length: i32,
    pub sequence_offset: i32,
    pub range: SSeqRange,
    pub chunk: i32,
}
/// Default eNoSubjMasking backup; immutable original storage stays with caller.
pub fn s_backup_subject(length: i32, backup: &mut SubjectSplitStruct, subject: &mut SubjectChunk) {
    if backup.initialized {
        return;
    }
    backup.initialized = true;
    backup.full_range = SSeqRange {
        left: 0,
        right: length,
    };
    backup.offset = 0;
    backup.next = 0;
    subject.chunk = -1;
}
pub fn s_get_next_subject_chunk(
    subject: &mut SubjectChunk,
    backup: &mut SubjectSplitStruct,
    is_nucleotide: bool,
    chunk_overlap: i32,
) -> i16 {
    let dbseq_chunk_overlap = if chunk_overlap > 0 {
        chunk_overlap
    } else {
        DBSEQ_CHUNK_OVERLAP
    };
    if backup.next >= backup.full_range.right {
        return SUBJECT_SPLIT_DONE;
    }
    let residual = if is_nucleotide { backup.next % 4 } else { 0 };
    backup.offset = backup.next - residual;
    subject.sequence_offset = if is_nucleotide {
        backup.offset / 4
    } else {
        backup.offset
    };
    if backup.offset as i64 + (MAX_DBSEQ_LEN as i64) < backup.full_range.right as i64 {
        subject.length = MAX_DBSEQ_LEN as i32;
        backup.next = backup.offset + MAX_DBSEQ_LEN as i32 - dbseq_chunk_overlap;
    } else {
        subject.length = backup.full_range.right - backup.offset;
        backup.next = backup.full_range.right;
    }
    subject.chunk += 1;
    if backup.offset == 0 && residual == 0 && backup.next == backup.full_range.right {
        subject.range = backup.full_range;
        return SUBJECT_SPLIT_OK;
    }
    subject.range = SSeqRange {
        left: residual,
        right: subject.length,
    };
    SUBJECT_SPLIT_OK
}
pub fn s_hsp_start_diag(hsp: &BlastHSP) -> i32 {
    hsp.query.offset - hsp.subject.offset
}
pub fn s_hsp_end_diag(hsp: &BlastHSP) -> i32 {
    hsp.query.end - hsp.subject.end
}
fn contained(a: i32, b: i32, c: i32, d: i32, e: i32, f: i32) -> bool {
    (a <= c && b >= c) && (d <= f && e >= f)
}
pub fn s_blast_merge_two_hs_ps(hsp1: &mut BlastHSP, hsp2: &BlastHSP, allow_gap: bool) -> bool {
    assert!(hsp1.gap_info.is_none() || hsp2.gap_info.is_none());
    if !allow_gap
        && hsp1.subject.offset - hsp2.subject.offset - hsp1.query.offset + hsp2.query.offset != 0
    {
        return false;
    }
    if hsp1.subject.frame != hsp2.subject.frame {
        return false;
    }
    if contained(
        hsp1.query.offset,
        hsp1.query.end,
        hsp2.query.offset,
        hsp1.subject.offset,
        hsp1.subject.end,
        hsp2.subject.offset,
    ) || contained(
        hsp1.query.offset,
        hsp1.query.end,
        hsp2.query.end,
        hsp1.subject.offset,
        hsp1.subject.end,
        hsp2.subject.end,
    ) {
        let score_density = (hsp1.score + hsp2.score) as f64 * 1.0
            / ((hsp1.query.end - hsp1.query.offset) + (hsp2.query.end - hsp2.query.offset)) as f64;
        hsp1.query.offset = hsp1.query.offset.min(hsp2.query.offset);
        hsp1.subject.offset = hsp1.subject.offset.min(hsp2.subject.offset);
        hsp1.query.end = hsp1.query.end.max(hsp2.query.end);
        hsp1.subject.end = hsp1.subject.end.max(hsp2.subject.end);
        if hsp2.score > hsp1.score {
            hsp1.query.gapped_start = hsp2.query.gapped_start;
            hsp1.subject.gapped_start = hsp2.subject.gapped_start;
            hsp1.score = hsp2.score;
        }
        hsp1.score =
            ((score_density * (hsp1.query.end - hsp1.query.offset) as f64) as i32).max(hsp1.score);
        return true;
    }
    false
}
pub fn blast_hsp_lists_merge(
    hsp_list_ptr: &mut Option<Box<BlastHSPList>>,
    combined_hsp_list_ptr: &mut Option<Box<BlastHSPList>>,
    hsp_num_max: i32,
    split_offsets: &[i32],
    contexts_per_query: i32,
    chunk_overlap_size: i32,
    allow_gap: bool,
    short_reads: bool,
) -> i16 {
    assert!(contexts_per_query != 0);
    assert!(
        !short_reads,
        "mapping/short-read merge is outside standard blastn/tblastn"
    );
    if hsp_list_ptr.as_ref().is_none_or(|list| list.hspcnt == 0) {
        return 0;
    }
    if combined_hsp_list_ptr.is_none() {
        *combined_hsp_list_ptr = hsp_list_ptr.take();
        return 0;
    }
    let combined = combined_hsp_list_ptr.as_mut().unwrap();
    let list = hsp_list_ptr.as_mut().unwrap();
    let (mut hspcnt1, mut hspcnt2) = (0usize, 0usize);
    for index1 in 0..combined.hspcnt as usize {
        let hsp = combined.hsp_array[index1].as_ref().unwrap();
        let overlaps = if contexts_per_query < 0 {
            hsp.subject.end > split_offsets[0]
        } else {
            let offset = split_offsets[(hsp.context % contexts_per_query) as usize];
            offset >= 0
                && ((hsp.query.frame >= 0 && hsp.query.end > offset)
                    || (hsp.query.frame < 0 && hsp.query.offset < offset + chunk_overlap_size))
        };
        if overlaps {
            combined.hsp_array.swap(hspcnt1, index1);
            hspcnt1 += 1;
        }
    }
    for index2 in 0..list.hspcnt as usize {
        let hsp = list.hsp_array[index2].as_ref().unwrap();
        let overlaps = if contexts_per_query < 0 {
            hsp.subject.offset < split_offsets[0] + chunk_overlap_size
        } else {
            let offset = split_offsets[(hsp.context % contexts_per_query) as usize];
            offset >= 0
                && ((hsp.query.frame < 0 && hsp.query.end > offset)
                    || (hsp.query.frame >= 0 && hsp.query.offset < offset + chunk_overlap_size))
        };
        if overlaps {
            list.hsp_array.swap(hspcnt2, index2);
            hspcnt2 += 1;
        }
    }
    if hspcnt1 > 0 && hspcnt2 > 0 {
        for index1 in 0..hspcnt1 {
            let hsp1 = combined.hsp_array[index1].as_mut().unwrap();
            for index2 in 0..hspcnt2 {
                let Some(hsp2) = list.hsp_array[index2].as_ref() else {
                    continue;
                };
                if hsp1.context != hsp2.context {
                    continue;
                }
                let (end_diag, start_diag) = if contexts_per_query < 0 || hsp1.query.frame >= 0 {
                    (s_hsp_end_diag(hsp1), s_hsp_start_diag(hsp2))
                } else {
                    (s_hsp_end_diag(hsp2), s_hsp_start_diag(hsp1))
                };
                if (end_diag - start_diag).abs() < OVERLAP_DIAG_CLOSE
                    && s_blast_merge_two_hs_ps(hsp1, hsp2, allow_gap)
                {
                    list.hsp_array[index2] = None;
                }
            }
        }
        blast_hsp_list_purge_null_hs_ps(Some(list));
    }
    let mut new_hspcnt = (list.hspcnt + combined.hspcnt).min(hsp_num_max);
    if new_hspcnt >= combined.allocated - 1 && !combined.do_not_reallocate {
        let new_allocated = (2 * new_hspcnt).min(hsp_num_max);
        if new_allocated > combined.allocated {
            let additional = new_allocated as usize - combined.hsp_array.len();
            if combined.hsp_array.try_reserve_exact(additional).is_err() {
                combined.do_not_reallocate = true;
            } else {
                combined
                    .hsp_array
                    .resize_with(new_allocated as usize, || None);
                combined.allocated = new_allocated;
            }
        } else {
            combined.do_not_reallocate = true;
        }
        new_hspcnt = new_hspcnt.min(combined.allocated);
    }
    s_blast_hsp_lists_combine_by_score(list, combined, new_hspcnt);
    *hsp_list_ptr = None;
    0
}

impl Default for SubjectSplitStruct {
    fn default() -> Self {
        Self {
            initialized: false,
            full_range: SSeqRange { left: 0, right: 0 },
            offset: 0,
            next: 0,
        }
    }
}
impl Default for SubjectChunk {
    fn default() -> Self {
        Self {
            length: 0,
            sequence_offset: 0,
            range: SSeqRange { left: 0, right: 0 },
            chunk: 0,
        }
    }
}
