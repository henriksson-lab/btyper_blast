// NCBI BLAST+ 2.17.0+, selected aa_ungapped.c definitions.
// Authors: NCBI BLAST developers. See LICENSE-NCBI.
use super::blast_aalookup::BlastAaLookupTable;
use super::blast_aascan::AminoAcidSubject;
use super::blast_diagnostics::{blast_ungapped_stats_update, BlastUngappedStats};
use super::blast_extend::{
    blast_extend_word_exit, blast_init_hit_list_sort_by_score, blast_save_init_hsp,
    BlastExtendWord, BlastInitHitList,
};
use super::blast_nascan::BlastOffsetPair;
use super::blast_parameters::BlastInitialWordParameters;
use super::blast_query_info::{b_search_context_info, BlastQueryInfo};

/// The sequence-block fields read by ordinary protein extension. Storage can
/// include a readable tail while length remains the original logical length.
#[derive(Clone, Copy, Debug)]
pub struct ProteinSequence<'a> {
    pub sequence: &'a [u8],
    pub length: i32,
}

pub fn s_blast_aa_extend_right(
    matrix: &[Vec<i32>],
    subject: &ProteinSequence<'_>,
    query: &ProteinSequence<'_>,
    s_off: i32,
    q_off: i32,
    dropoff: i32,
    length: &mut i32,
    mut maxscore: i32,
    s_last_off: &mut i32,
) -> i32 {
    let mut best_i = -1;
    let mut score = maxscore;
    let n = (subject.length - s_off).min(query.length - q_off);
    let mut i = 0;
    while i < n {
        score += matrix[query.sequence[(q_off + i) as usize] as usize]
            [subject.sequence[(s_off + i) as usize] as usize];
        if score > maxscore {
            maxscore = score;
            best_i = i;
        }
        if score <= 0 || maxscore - score >= dropoff {
            break;
        }
        i += 1;
    }
    *length = best_i + 1;
    *s_last_off = s_off + i;
    maxscore
}
pub fn s_blast_aa_extend_left(
    matrix: &[Vec<i32>],
    subject: &ProteinSequence<'_>,
    query: &ProteinSequence<'_>,
    s_off: i32,
    q_off: i32,
    dropoff: i32,
    length: &mut i32,
    mut maxscore: i32,
) -> i32 {
    let mut score = maxscore;
    let n = s_off.min(q_off);
    let mut best_i = n + 1;
    let s = s_off - n;
    let q = q_off - n;
    let mut i = n;
    while i >= 0 {
        score += matrix[query.sequence[(q + i) as usize] as usize]
            [subject.sequence[(s + i) as usize] as usize];
        if score > maxscore {
            maxscore = score;
            best_i = i;
        }
        if maxscore - score >= dropoff {
            break;
        }
        i -= 1;
    }
    *length = n - best_i + 1;
    maxscore
}
pub fn s_blast_aa_extend_one_hit(
    matrix: &[Vec<i32>],
    subject: &ProteinSequence<'_>,
    query: &ProteinSequence<'_>,
    s_off: i32,
    q_off: i32,
    dropoff: i32,
    hsp_q: &mut i32,
    hsp_s: &mut i32,
    hsp_len: &mut i32,
    word_size: i32,
    use_pssm: bool,
    s_last_off: &mut i32,
) -> i32 {
    assert!(!use_pssm, "PSSM/RPS extension is outside ordinary tblastn");
    let mut score = 0;
    let mut sum = 0;
    let mut left_disp = 0;
    let mut right_disp = 0;
    let mut q_left_off = q_off;
    let mut q_right_off = q_off + word_size;
    let mut q_best_left_off = q_off;
    for i in 0..word_size {
        sum += matrix[query.sequence[(q_off + i) as usize] as usize]
            [subject.sequence[(s_off + i) as usize] as usize];
        if sum > score {
            score = sum;
            q_best_left_off = q_left_off;
            q_right_off = q_off + i;
        } else if sum <= 0 {
            sum = 0;
            q_left_off = q_off + i + 1;
        }
    }
    let init_hit_width = q_right_off - q_left_off + 1;
    q_left_off = q_best_left_off;
    let s_left_off = q_left_off + (s_off - q_off);
    let s_right_off = q_right_off + (s_off - q_off);
    let left_score = s_blast_aa_extend_left(
        matrix,
        subject,
        query,
        s_left_off - 1,
        q_left_off - 1,
        dropoff,
        &mut left_disp,
        score,
    );
    let total_score = s_blast_aa_extend_right(
        matrix,
        subject,
        query,
        s_right_off + 1,
        q_right_off + 1,
        dropoff,
        &mut right_disp,
        left_score,
        s_last_off,
    );
    *hsp_q = q_left_off - left_disp;
    *hsp_s = s_left_off - left_disp;
    *hsp_len = left_disp + right_disp + init_hit_width;
    total_score
}
pub fn s_blast_aa_extend_two_hit(
    matrix: &[Vec<i32>],
    subject: &ProteinSequence<'_>,
    query: &ProteinSequence<'_>,
    s_left_off: i32,
    mut s_right_off: i32,
    mut q_right_off: i32,
    dropoff: i32,
    hsp_q: &mut i32,
    hsp_s: &mut i32,
    hsp_len: &mut i32,
    use_pssm: bool,
    word_size: i32,
    right_extend: &mut bool,
    s_last_off: &mut i32,
) -> i32 {
    assert!(!use_pssm, "PSSM/RPS extension is outside ordinary tblastn");
    let mut left_d = 0;
    let mut right_d = 0;
    let mut left_score = 0;
    let mut right_score = 0;
    let mut score = 0;
    for i in 0..word_size {
        score += matrix[query.sequence[(q_right_off + i) as usize] as usize]
            [subject.sequence[(s_right_off + i) as usize] as usize];
        if score > left_score {
            left_score = score;
            right_d = i + 1;
        }
    }
    q_right_off += right_d;
    s_right_off += right_d;
    right_d = 0;
    *right_extend = false;
    *s_last_off = s_right_off;
    left_score = s_blast_aa_extend_left(
        matrix,
        subject,
        query,
        s_right_off - 1,
        q_right_off - 1,
        dropoff,
        &mut left_d,
        0,
    );
    if left_d >= s_right_off - s_left_off {
        *right_extend = true;
        right_score = s_blast_aa_extend_right(
            matrix,
            subject,
            query,
            s_right_off,
            q_right_off,
            dropoff,
            &mut right_d,
            left_score,
            s_last_off,
        );
    }
    *hsp_q = q_right_off - left_d;
    *hsp_s = s_right_off - left_d;
    *hsp_len = left_d + right_d;
    left_score.max(right_score)
}

pub fn s_blast_aa_word_finder_two_hit(
    subject: &AminoAcidSubject<'_>,
    subject_length: i32,
    query: &ProteinSequence<'_>,
    lookup: &BlastAaLookupTable,
    ewp: &mut BlastExtendWord,
    matrix: &[Vec<i32>],
    word_params: &BlastInitialWordParameters,
    query_info: &BlastQueryInfo,
    offset_pairs: &mut [BlastOffsetPair],
    array_size: i32,
    ungapped_hsps: &mut BlastInitHitList,
    ungapped_stats: Option<&mut BlastUngappedStats>,
) -> i16 {
    assert!(
        !lookup.use_pssm,
        "PSSM/RPS lookup is outside ordinary tblastn"
    );
    assert!(
        array_size >= lookup.longest_chain && array_size as usize <= offset_pairs.len(),
        "word-finder buffer must hold its longest hit chain"
    );
    let diag = ewp
        .diag_table
        .as_mut()
        .expect("protein word finder requires original diagonal array");
    let diag_offset = diag.offset;
    let diag_mask = diag.diag_mask;
    let window = diag.window;
    let wordsize = lookup.word_length;
    let scansub = lookup
        .scansub_callback
        .expect("select the original protein scanner first");
    let mut totalhits = 0i32;
    let mut hits_extended = 0i32;
    let mut hsp_q = 0;
    let mut hsp_s = 0;
    let mut hsp_len = 0;
    let subject_block = ProteinSequence {
        sequence: subject.sequence,
        length: subject_length,
    };
    let mut scan_range = [
        0,
        subject.seq_ranges[0].left,
        subject.seq_ranges[0].right - wordsize,
    ];
    if scan_range[2] < scan_range[1] {
        scan_range[2] = scan_range[1];
    }
    while scan_range[1] <= scan_range[2] {
        let hits = scansub(lookup, subject, offset_pairs, array_size, &mut scan_range);
        totalhits = totalhits.wrapping_add(hits);
        for pair in &offset_pairs[..hits as usize] {
            let query_offset = pair.q_off as u32;
            let subject_offset = pair.s_off as u32;
            let diag_coord =
                (query_offset.wrapping_sub(subject_offset) & diag_mask as u32) as usize;
            if diag.hit_level_array[diag_coord].flag() != 0 {
                if (subject_offset.wrapping_add(diag_offset as u32) as i32)
                    < diag.hit_level_array[diag_coord].last_hit()
                {
                    continue;
                } else {
                    diag.hit_level_array[diag_coord]
                        .set_last_hit(subject_offset.wrapping_add(diag_offset as u32) as i32);
                    diag.hit_level_array[diag_coord].set_flag(0);
                }
            } else {
                let last_hit = diag.hit_level_array[diag_coord].last_hit() - diag_offset;
                let diff = subject_offset.wrapping_sub(last_hit as u32) as i32;
                if diff >= window {
                    diag.hit_level_array[diag_coord]
                        .set_last_hit(subject_offset.wrapping_add(diag_offset as u32) as i32);
                    continue;
                }
                if diff < wordsize {
                    continue;
                }
                let curr_context = b_search_context_info(query_offset as i32, query_info) as usize;
                if query_offset.wrapping_sub(diff as u32)
                    < query_info.contexts[curr_context].query_offset as u32
                {
                    diag.hit_level_array[diag_coord]
                        .set_last_hit(subject_offset.wrapping_add(diag_offset as u32) as i32);
                    continue;
                }
                let cutoffs = &word_params.cutoffs[curr_context];
                let mut right_extend = false;
                let mut s_last_off = 0;
                let score = s_blast_aa_extend_two_hit(
                    matrix,
                    &subject_block,
                    query,
                    last_hit + wordsize,
                    subject_offset as i32,
                    query_offset as i32,
                    cutoffs.x_dropoff,
                    &mut hsp_q,
                    &mut hsp_s,
                    &mut hsp_len,
                    false,
                    wordsize,
                    &mut right_extend,
                    &mut s_last_off,
                );
                hits_extended = hits_extended.wrapping_add(1);
                if score >= cutoffs.cutoff_score {
                    blast_save_init_hsp(
                        ungapped_hsps,
                        hsp_q,
                        hsp_s,
                        query_offset as i32,
                        subject_offset as i32,
                        hsp_len,
                        score,
                    );
                }
                if right_extend {
                    diag.hit_level_array[diag_coord].set_flag(1);
                    diag.hit_level_array[diag_coord]
                        .set_last_hit(s_last_off - (wordsize - 1) + diag_offset);
                } else {
                    diag.hit_level_array[diag_coord]
                        .set_last_hit(subject_offset.wrapping_add(diag_offset as u32) as i32);
                }
            }
        }
    }
    blast_extend_word_exit(Some(ewp), subject_length);
    blast_ungapped_stats_update(
        ungapped_stats,
        totalhits,
        hits_extended,
        ungapped_hsps.total,
    );
    0
}
pub fn s_blast_aa_word_finder_one_hit(
    subject: &AminoAcidSubject<'_>,
    subject_length: i32,
    query: &ProteinSequence<'_>,
    lookup: &BlastAaLookupTable,
    ewp: &mut BlastExtendWord,
    matrix: &[Vec<i32>],
    word_params: &BlastInitialWordParameters,
    query_info: &BlastQueryInfo,
    offset_pairs: &mut [BlastOffsetPair],
    array_size: i32,
    ungapped_hsps: &mut BlastInitHitList,
    ungapped_stats: Option<&mut BlastUngappedStats>,
) -> i16 {
    assert!(
        !lookup.use_pssm,
        "PSSM/RPS lookup is outside ordinary tblastn"
    );
    assert!(
        array_size >= lookup.longest_chain && array_size as usize <= offset_pairs.len(),
        "word-finder buffer must hold its longest hit chain"
    );
    let diag = ewp
        .diag_table
        .as_mut()
        .expect("protein word finder requires original diagonal array");
    let diag_offset = diag.offset;
    let diag_mask = diag.diag_mask;
    let scansub = lookup
        .scansub_callback
        .expect("select the original protein scanner first");
    let wordsize = lookup.word_length;
    let mut totalhits = 0i32;
    let mut hits_extended = 0i32;
    let mut hsp_q = 0;
    let mut hsp_s = 0;
    let mut hsp_len = 0;
    let subject_block = ProteinSequence {
        sequence: subject.sequence,
        length: subject_length,
    };
    let mut scan_range = [
        0,
        subject.seq_ranges[0].left,
        subject.seq_ranges[0].right - wordsize,
    ];
    while scan_range[1] <= scan_range[2] {
        let hits = scansub(lookup, subject, offset_pairs, array_size, &mut scan_range);
        totalhits = totalhits.wrapping_add(hits);
        for pair in &offset_pairs[..hits as usize] {
            let query_offset = pair.q_off as u32;
            let subject_offset = pair.s_off as u32;
            let diag_coord =
                (subject_offset.wrapping_sub(query_offset) & diag_mask as u32) as usize;
            let diff = subject_offset
                .wrapping_sub((diag.hit_level_array[diag_coord].last_hit() - diag_offset) as u32)
                as i32;
            if diff >= 0 {
                let curr_context = b_search_context_info(query_offset as i32, query_info) as usize;
                let cutoffs = &word_params.cutoffs[curr_context];
                let mut s_last_off = 0;
                let score = s_blast_aa_extend_one_hit(
                    matrix,
                    &subject_block,
                    query,
                    subject_offset as i32,
                    query_offset as i32,
                    cutoffs.x_dropoff,
                    &mut hsp_q,
                    &mut hsp_s,
                    &mut hsp_len,
                    wordsize,
                    false,
                    &mut s_last_off,
                );
                if score >= cutoffs.cutoff_score {
                    blast_save_init_hsp(
                        ungapped_hsps,
                        hsp_q,
                        hsp_s,
                        query_offset as i32,
                        subject_offset as i32,
                        hsp_len,
                        score,
                    );
                }
                diag.hit_level_array[diag_coord]
                    .set_last_hit(s_last_off - (wordsize - 1) + diag_offset);
                hits_extended = hits_extended.wrapping_add(1);
            }
        }
    }
    blast_extend_word_exit(Some(ewp), subject_length);
    blast_ungapped_stats_update(
        ungapped_stats,
        totalhits,
        hits_extended,
        ungapped_hsps.total,
    );
    0
}
/// Ordinary protein word finder. Subject codes may be a translated nucleotide
/// frame; length/ranges stay in that protein frame's coordinate system.
pub fn blast_aa_word_finder(
    subject: &AminoAcidSubject<'_>,
    subject_length: i32,
    query: &ProteinSequence<'_>,
    query_info: &BlastQueryInfo,
    lookup: &BlastAaLookupTable,
    matrix: &[Vec<i32>],
    word_params: &BlastInitialWordParameters,
    ewp: &mut BlastExtendWord,
    offset_pairs: &mut [BlastOffsetPair],
    offset_array_size: i32,
    init_hitlist: &mut BlastInitHitList,
    ungapped_stats: Option<&mut BlastUngappedStats>,
) -> i16 {
    let status = if ewp
        .diag_table
        .as_ref()
        .expect("protein diagonal array required")
        .multiple_hits
    {
        s_blast_aa_word_finder_two_hit(
            subject,
            subject_length,
            query,
            lookup,
            ewp,
            matrix,
            word_params,
            query_info,
            offset_pairs,
            offset_array_size,
            init_hitlist,
            ungapped_stats,
        )
    } else {
        s_blast_aa_word_finder_one_hit(
            subject,
            subject_length,
            query,
            lookup,
            ewp,
            matrix,
            word_params,
            query_info,
            offset_pairs,
            offset_array_size,
            init_hitlist,
            ungapped_stats,
        )
    };
    blast_init_hit_list_sort_by_score(init_hitlist);
    status
}
