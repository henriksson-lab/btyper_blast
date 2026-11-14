// NCBI BLAST+ 2.17.0+, blast_nalookup.c and blast_options.h.
// Authors: NCBI BLAST developers. See LICENSE-NCBI.

use super::blast_lookup::{blast_lookup_index_query_exact_matches, SSeqRange};
use super::blast_nascan::NaScanSubject;
use super::blast_options::{s_blast_filter_options_mask_at_hash, QuerySetUpOptions};
use super::blast_program::{blast_program_is_mapping, EBlastProgramType};
use super::blast_util::blast_compress_blastna_sequence;
use super::lookup_util::ilog2;

pub const NA_HITS_PER_CELL: usize = 3;

/// Contiguous megablast table fields used by the translated scan kernels.
/// Query offsets are one-based in the heads/links, with zero as terminator.
#[derive(Debug, Default)]
pub struct BlastMBLookupTable {
    pub word_length: i32,
    pub lut_word_length: i32,
    pub scan_step: i32,
    pub hashsize: i64,
    pub longest_chain: i32,
    pub hashtable: Vec<i32>,
    pub next_pos: Vec<i32>,
    pub pv_array: Vec<u32>,
    pub pv_array_bts: i32,
    pub stride: bool,
    pub masked_locations: Vec<SSeqRange>,
}

/// Matches C `s_RemovePolyAWords`. Mapping is outside the search presets;
/// this helper is retained as a dependency of the unchanged filling function.
/// The inherited width-16 branch is outside the verified domain.
pub fn s_remove_poly_a_words(mb_lt: &mut BlastMBLookupTable) -> i16 {
    let word_size = mb_lt.lut_word_length;
    mb_lt.hashtable[0] = 0;
    mb_lt.hashtable[((1_i64 << (2 * word_size)) - 1) as usize] = 0;
    if word_size < 16 {
        return 0;
    }
    assert_eq!(word_size, 16);
    for i in 1_i64..4 {
        for k in 0..word_size {
            mb_lt.hashtable[(i << (k * 2)) as usize] = 0;
        }
    }
    for i in 0_u32..3 {
        for k in 0..word_size {
            let word = (0xffffffff_u32 ^ (3_u32 << (k * 2))) | (i << (k * 2));
            mb_lt.hashtable[word as usize] = 0;
        }
    }
    0
}

/// Matches C `s_FillContigMBTable`. Query bytes exclude the leading sentinel;
/// signed offsets below retain the original sequence_start-relative indexing.
/// Included intervals should be ordered/disjoint, as produced by BLAST setup.
pub fn s_fill_contig_mb_table(
    query: &[u8],
    location: &[SSeqRange],
    mb_lt: &mut BlastMBLookupTable,
    lookup_options: &LookupTableOptions,
    counts: Option<&[u8]>,
) -> i16 {
    let k_lut_word_length = mb_lt.lut_word_length;
    let k_lut_mask = mb_lt.hashsize - 1;
    let full_word_size = mb_lt.word_length;
    let k_compression_factor = 2048;
    let k_db_filter = lookup_options.db_filter;
    mb_lt.next_pos = vec![0; query.len() + 1];
    let mut helper_array = vec![0_u32; (mb_lt.hashsize / k_compression_factor) as usize];
    if k_db_filter {
        mb_lt.pv_array[..(mb_lt.hashsize >> mb_lt.pv_array_bts) as usize].fill(0);
    }
    for loc in location {
        let mut from = loc.left;
        let to = loc.right - k_lut_word_length;
        let mut ecode = 0_i64;
        let max_word_count = lookup_options.max_db_word_count;
        let mut shift = 0;
        let mut pos_shift = 0;
        if lookup_options.stride > 0 {
            shift = lookup_options.stride - 1;
            pos_shift = k_lut_word_length + 1;
        }
        if full_word_size > loc.right - loc.left + 1 {
            continue;
        }
        let mut seq = from;
        let mut pos = seq + k_lut_word_length;
        from -= k_lut_word_length - 2;
        let last_offset = to + 2;
        let mut index = from;
        while index <= last_offset {
            seq += 1;
            let val = query[(seq - 1) as usize];
            if val & 0xfc != 0 {
                ecode = 0;
                pos = seq + k_lut_word_length;
                index += 1;
                continue;
            }
            ecode = ((ecode << 2) & k_lut_mask) + i64::from(val);
            if seq < pos {
                index += 1;
                continue;
            }
            if k_db_filter {
                let counts = counts.expect("database word counts required by db_filter");
                if ecode & 1 == 0 {
                    if counts[(ecode / 2) as usize] >> 4 >= max_word_count {
                        index += 1;
                        continue;
                    }
                } else if counts[(ecode / 2) as usize] & 0xf >= max_word_count {
                    index += 1;
                    continue;
                }
            }
            if mb_lt.hashtable[ecode as usize] == 0 {
                mb_lt.pv_array[(ecode >> mb_lt.pv_array_bts) as usize] |= 1_u32 << (ecode & 31);
            } else {
                let bucket = (ecode / k_compression_factor) as usize;
                helper_array[bucket] = helper_array[bucket].wrapping_add(1);
            }
            mb_lt.next_pos[index as usize] = mb_lt.hashtable[ecode as usize];
            mb_lt.hashtable[ecode as usize] = index;
            index += shift;
            seq += shift;
            pos = seq + pos_shift;
            index += 1;
        }
    }
    if blast_program_is_mapping(lookup_options.program_number) {
        s_remove_poly_a_words(mb_lt);
    }
    let mut longest_chain = 2;
    for count in helper_array {
        longest_chain = longest_chain.max(count);
    }
    mb_lt.longest_chain = longest_chain as i32;
    0
}

/// C `BlastMBLookupTableNew`, restricted to its supported contiguous,
/// non-database-filtered path at widths 9..12. Discontiguous templates and
/// mapper database filtering are rejected rather than translating their
/// unused dependencies.
pub fn blast_mb_lookup_table_new(
    query: &[u8],
    location: &[SSeqRange],
    lookup_options: &LookupTableOptions,
    query_options: Option<&QuerySetUpOptions>,
    approx_table_entries: i32,
    lut_width: i32,
) -> Result<BlastMBLookupTable, i16> {
    assert_eq!(
        lookup_options.mb_template_length, 0,
        "discontiguous lookup is outside scope"
    );
    assert!(
        !lookup_options.db_filter,
        "database word-count filtering is outside scope"
    );
    assert!(
        !blast_program_is_mapping(lookup_options.program_number),
        "mapping is outside scope"
    );
    if location.is_empty() {
        return Err(-1);
    }
    assert!(
        (9..=12).contains(&lut_width),
        "lookup widths outside 9..12 are outside scope"
    );
    let mut mb_lt = BlastMBLookupTable {
        word_length: lookup_options.word_size,
        stride: lookup_options.stride > 0,
        lut_word_length: lut_width,
        hashsize: 1_i64 << (2 * lut_width),
        ..Default::default()
    };
    mb_lt.hashtable = vec![0; mb_lt.hashsize as usize];
    if mb_lt.word_length > mb_lt.lut_word_length && s_has_mask_at_hash_enabled(query_options) {
        mb_lt.masked_locations = s_seq_loc_list_invert(location, query.len() as i32);
    }
    let k_target_pv_size = 131072;
    let mut pv_size = if mb_lt.lut_word_length <= 12 {
        if mb_lt.hashsize <= 8 * k_target_pv_size {
            (mb_lt.hashsize >> 5) as i32
        } else {
            (k_target_pv_size / 4) as i32
        }
    } else {
        (k_target_pv_size * 64 / 4) as i32
    };
    if approx_table_entries <= 15000 || approx_table_entries >= 800000 {
        pv_size /= 2;
    }
    mb_lt.pv_array_bts = ilog2(mb_lt.hashsize / i64::from(pv_size));
    mb_lt.pv_array = vec![0; pv_size as usize];
    mb_lt.scan_step = mb_lt.word_length - mb_lt.lut_word_length + 1;
    let status = s_fill_contig_mb_table(query, location, &mut mb_lt, lookup_options, None);
    if status != 0 {
        return Err(-1);
    }
    Ok(mb_lt)
}

/// The C payload union occupies three words. For overflow cells, word zero
/// stores the cursor and the other words retain their calloc initialization.
#[derive(Clone, Copy, Debug, Default)]
#[repr(C)]
pub struct NaLookupBackboneCell {
    pub num_used: i32,
    pub payload: [i32; NA_HITS_PER_CELL],
}

#[derive(Debug, Default)]
pub struct BlastNaLookupTable {
    pub mask: i32,
    pub word_length: i32,
    pub lut_word_length: i32,
    pub scan_step: i32,
    pub backbone_size: i32,
    pub longest_chain: i32,
    pub thick_backbone: Vec<NaLookupBackboneCell>,
    pub overflow: Vec<i32>,
    pub overflow_size: i32,
    pub pv: Vec<u32>,
    pub masked_locations: Vec<SSeqRange>,
    pub scansub_callback: Option<NaScanSubject>,
}

/// Matches C `s_HasMaskAtHashEnabled`, including the deprecated string path.
pub fn s_has_mask_at_hash_enabled(query_options: Option<&QuerySetUpOptions>) -> bool {
    let Some(query_options) = query_options else {
        return false;
    };
    if s_blast_filter_options_mask_at_hash(query_options.filtering_options.as_ref()) {
        return true;
    }
    if let Some(filter_string) = &query_options.filter_string {
        if filter_string
            .split(|byte| *byte == 0)
            .next()
            .unwrap()
            .contains(&b'm')
        {
            return true;
        }
    }
    false
}

/// Matches C `s_BlastNaLookupFinalize`: at most three offsets live directly
/// in a cell; longer lists retain their cursor in a contiguous overflow array.
pub fn s_blast_na_lookup_finalize(
    thin_backbone: &mut [Option<Vec<i32>>],
    lookup: &mut BlastNaLookupTable,
) {
    let mut overflow_cells_needed = 0;
    let mut overflow_cursor = 0;
    let mut longest_chain = 0;
    lookup.thick_backbone = vec![NaLookupBackboneCell::default(); lookup.backbone_size as usize];
    lookup.pv = vec![0; (lookup.backbone_size as usize >> 5) + 1];
    for chain in thin_backbone
        .iter()
        .take(lookup.backbone_size as usize)
        .flatten()
    {
        let num_hits = chain[1];
        if num_hits > NA_HITS_PER_CELL as i32 {
            overflow_cells_needed += num_hits;
        }
        longest_chain = longest_chain.max(num_hits);
    }
    lookup.longest_chain = longest_chain;
    if overflow_cells_needed > 0 {
        lookup.overflow = vec![0; overflow_cells_needed as usize];
    }
    for (i, cell) in thin_backbone
        .iter_mut()
        .enumerate()
        .take(lookup.backbone_size as usize)
    {
        let Some(chain) = cell.take() else {
            continue;
        };
        let num_hits = chain[1];
        lookup.thick_backbone[i].num_used = num_hits;
        lookup.pv[i >> 5] |= 1_u32 << (i & 31);
        if num_hits <= NA_HITS_PER_CELL as i32 {
            for j in 0..num_hits {
                lookup.thick_backbone[i].payload[j as usize] = chain[j as usize + 2];
            }
        } else {
            lookup.thick_backbone[i].payload[0] = overflow_cursor;
            for j in 0..num_hits {
                lookup.overflow[overflow_cursor as usize] = chain[j as usize + 2];
                overflow_cursor += 1;
            }
        }
    }
    lookup.overflow_size = overflow_cursor;
}

/// Matches C `BlastNaLookupTableNew`. Allocation ownership replaces its
/// pointer output/status (the C function always returns zero on success).
pub fn blast_na_lookup_table_new(
    query: &[u8],
    locations: &[SSeqRange],
    opt: &LookupTableOptions,
    query_options: Option<&QuerySetUpOptions>,
    lut_width: i32,
) -> BlastNaLookupTable {
    let mut lookup = BlastNaLookupTable {
        word_length: opt.word_size,
        lut_word_length: lut_width,
        backbone_size: 1 << (2 * lut_width),
        mask: (1 << (2 * lut_width)) - 1,
        scan_step: opt.word_size - lut_width + 1,
        ..Default::default()
    };
    let mut thin_backbone = vec![None; lookup.backbone_size as usize];
    blast_lookup_index_query_exact_matches(
        &mut thin_backbone,
        lookup.word_length,
        2,
        lookup.lut_word_length,
        query,
        locations,
    );
    if !locations.is_empty()
        && lookup.word_length > lookup.lut_word_length
        && s_has_mask_at_hash_enabled(query_options)
    {
        lookup.masked_locations = s_seq_loc_list_invert(locations, query.len() as i32);
    }
    s_blast_na_lookup_finalize(&mut thin_backbone, &mut lookup);
    lookup
}

/// Original small-table storage. Callback selection is supplied by the caller
/// until the remaining selected table/scanner implementations are translated.
#[derive(Debug, Default)]
pub struct BlastSmallNaLookupTable {
    pub mask: i32,
    pub word_length: i32,
    pub lut_word_length: i32,
    pub scan_step: i32,
    pub backbone_size: i32,
    pub longest_chain: i32,
    pub final_backbone: Vec<i16>,
    pub overflow: Vec<i16>,
    pub overflow_size: i32,
    pub masked_locations: Vec<SSeqRange>,
}

/// Matches C `s_BlastSmallNaLookupFinalize`. Consumes the thin chains, retaining
/// the negative overflow cursors and trailing -1 terminators. The compressed
/// query is returned separately from the lookup table, as in the C query block.
pub fn s_blast_small_na_lookup_finalize(
    thin_backbone: &mut [Option<Vec<i32>>],
    lookup: &mut BlastSmallNaLookupTable,
    query: &[u8],
    compressed_query: &mut Option<Vec<u8>>,
) -> i32 {
    let mut overflow_cells_needed = 2;
    let mut overflow_cursor: i32 = 2;
    let mut longest_chain = 0;
    for chain in thin_backbone
        .iter()
        .take(lookup.backbone_size as usize)
        .flatten()
    {
        let num_hits = chain[1];
        if num_hits > 1 {
            overflow_cells_needed += num_hits + 1;
        }
        longest_chain = longest_chain.max(num_hits);
    }
    if overflow_cells_needed >= 32768 {
        for chain in thin_backbone.iter_mut().take(lookup.backbone_size as usize) {
            *chain = None;
        }
        return -1;
    }
    *compressed_query = Some(blast_compress_blastna_sequence(query));
    lookup.final_backbone = vec![0; lookup.backbone_size as usize];
    lookup.longest_chain = longest_chain;
    lookup.overflow = vec![0; overflow_cells_needed as usize];
    for (i, cell) in thin_backbone
        .iter_mut()
        .enumerate()
        .take(lookup.backbone_size as usize)
    {
        let Some(chain) = cell.take() else {
            lookup.final_backbone[i] = -1;
            continue;
        };
        let num_hits = chain[1];
        if num_hits == 1 {
            lookup.final_backbone[i] = chain[2] as i16;
        } else {
            lookup.final_backbone[i] = -overflow_cursor as i16;
            for j in 0..num_hits {
                lookup.overflow[overflow_cursor as usize] = chain[j as usize + 2] as i16;
                overflow_cursor += 1;
            }
            lookup.overflow[overflow_cursor as usize] = -1;
            overflow_cursor += 1;
        }
    }
    lookup.overflow_size = overflow_cursor;
    0
}

/// Matches C `s_SeqLocListInvert`, including its short leading-interval rule
/// and unconditional subsequent intervals, which can have left > right.
pub fn s_seq_loc_list_invert(locations: &[SSeqRange], length: i32) -> Vec<SSeqRange> {
    assert!(!locations.is_empty());
    let mut result = Vec::new();
    let stop = 0.max(locations[0].left - 1);
    if stop > 2 {
        result.push(SSeqRange {
            left: 0,
            right: stop,
        });
    }
    for (i, location) in locations.iter().enumerate() {
        let start = location.right + 1;
        let stop = if i + 1 < locations.len() {
            locations[i + 1].left - 1
        } else {
            length - 1
        };
        result.push(SSeqRange {
            left: start,
            right: stop,
        });
    }
    result
}

/// Matches C `BlastSmallNaLookupTableNew`. The unused query-options parameter
/// is omitted. Allocation ownership replaces the C destructor and free calls;
/// its -1 overflow-limit status is retained as the error value.
pub fn blast_small_na_lookup_table_new(
    query: &[u8],
    locations: &[SSeqRange],
    opt: &LookupTableOptions,
    lut_width: i32,
) -> Result<(BlastSmallNaLookupTable, Vec<u8>), i32> {
    let mut lookup = BlastSmallNaLookupTable {
        word_length: opt.word_size,
        lut_word_length: lut_width,
        backbone_size: 1 << (2 * lut_width),
        mask: (1 << (2 * lut_width)) - 1,
        scan_step: opt.word_size - lut_width + 1,
        ..Default::default()
    };
    let mut thin_backbone = vec![None; lookup.backbone_size as usize];
    blast_lookup_index_query_exact_matches(
        &mut thin_backbone,
        lookup.word_length,
        2,
        lookup.lut_word_length,
        query,
        locations,
    );
    if !locations.is_empty() && lookup.word_length > lookup.lut_word_length {
        lookup.masked_locations = s_seq_loc_list_invert(locations, query.len() as i32);
    }
    let mut compressed = None;
    let status =
        s_blast_small_na_lookup_finalize(&mut thin_backbone, &mut lookup, query, &mut compressed);
    if status != 0 {
        return Err(status);
    }
    Ok((lookup, compressed.unwrap()))
}

/// Original discriminants, including unsupported modes, for table selection.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[repr(i32)]
pub enum ELookupTableType {
    MBLookupTable = 0,
    SmallNaLookupTable = 1,
    NaLookupTable = 2,
    AaLookupTable = 3,
    CompressedAaLookupTable = 4,
    PhiLookupTable = 5,
    PhiNaLookupTable = 6,
    RPSLookupTable = 7,
    IndexedMBLookupTable = 8,
    MixedMBLookupTable = 9,
    NaHashLookupTable = 10,
}

/// Only the fields read by `BlastChooseNaLookupTable`. This is not yet the
/// full options structure or a public search preset.
#[derive(Clone, Copy, Debug)]
pub struct LookupTableOptions {
    pub word_size: i32,
    pub mb_template_length: i32,
    pub program_number: EBlastProgramType,
    pub db_filter: bool,
    pub stride: i32,
    pub max_db_word_count: u8,
}

/// Matches C `BlastChooseNaLookupTable`, retaining every width threshold and
/// small-table offset limit. Selection does not imply all table builders exist.
pub fn blast_choose_na_lookup_table(
    lookup_options: &LookupTableOptions,
    approx_table_entries: i32,
    max_q_off: i32,
    lut_width: &mut i32,
) -> ELookupTableType {
    use ELookupTableType::*;
    assert!(lookup_options.word_size >= 4);
    if lookup_options.mb_template_length > 0 {
        *lut_width = lookup_options.word_size;
        return MBLookupTable;
    }
    if blast_program_is_mapping(lookup_options.program_number)
        && lookup_options.word_size >= 16
        && lookup_options.db_filter
    {
        *lut_width = 16;
        return NaHashLookupTable;
    }
    let (mut lut_type, width) = match lookup_options.word_size {
        4..=6 => (SmallNaLookupTable, lookup_options.word_size),
        7 => (
            SmallNaLookupTable,
            if approx_table_entries < 250 { 6 } else { 7 },
        ),
        8 => (
            SmallNaLookupTable,
            if approx_table_entries < 8500 { 7 } else { 8 },
        ),
        9 => {
            if approx_table_entries < 1250 {
                (SmallNaLookupTable, 7)
            } else if approx_table_entries < 21000 {
                (SmallNaLookupTable, 8)
            } else {
                (MBLookupTable, 9)
            }
        }
        10 => {
            if approx_table_entries < 1250 {
                (SmallNaLookupTable, 7)
            } else if approx_table_entries < 8500 {
                (SmallNaLookupTable, 8)
            } else if approx_table_entries < 18000 {
                (MBLookupTable, 9)
            } else {
                (MBLookupTable, 10)
            }
        }
        11 => {
            if approx_table_entries < 12000 {
                (SmallNaLookupTable, 8)
            } else if approx_table_entries < 180000 {
                (MBLookupTable, 10)
            } else {
                (MBLookupTable, 11)
            }
        }
        12 => {
            if approx_table_entries < 8500 {
                (SmallNaLookupTable, 8)
            } else if approx_table_entries < 18000 {
                (MBLookupTable, 9)
            } else if approx_table_entries < 60000 {
                (MBLookupTable, 10)
            } else if approx_table_entries < 900000 {
                (MBLookupTable, 11)
            } else {
                (MBLookupTable, 12)
            }
        }
        _ => {
            if approx_table_entries < 8500 {
                (SmallNaLookupTable, 8)
            } else if approx_table_entries < 300000 {
                (MBLookupTable, 11)
            } else {
                (MBLookupTable, 12)
            }
        }
    };
    *lut_width = width;
    if lut_type == SmallNaLookupTable && (approx_table_entries >= 32767 || max_q_off >= 32768) {
        lut_type = NaLookupTable;
    }
    lut_type
}
