// NCBI BLAST+ 2.17.0+, selected blast_aalookup.c/h definitions.
// Authors: NCBI BLAST developers. See LICENSE-NCBI.
use super::blast_aascan::AaScanSubject;
use super::blast_encoding::BLASTAA_SIZE;
use super::blast_lookup::{
    blast_lookup_add_word_hit, blast_lookup_index_query_exact_matches, SSeqRange, ThinBackbone,
};
use super::lookup_util::ilog2;

pub const AA_HITS_PER_CELL: usize = 3;
pub const BLAST_WORDSIZE_PROT: i32 = 3;
pub const BLAST_WORD_THRESHOLD_TBLASTN: i32 = 13;
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
#[repr(i32)]
pub enum EBoneType {
    #[default]
    Backbone = 0,
    Smallbone = 1,
}
#[derive(Clone, Copy, Debug, Default)]
#[repr(C)]
pub struct AaLookupBackboneCell {
    pub num_used: i32,
    pub payload: [i32; AA_HITS_PER_CELL],
}
/// Original small cell is 12 bytes/alignment 4: a u16 count, two padding bytes,
/// then an eight-byte union. Native-byte accessors preserve both union views.
#[derive(Clone, Copy, Debug, Default)]
#[repr(C, align(4))]
pub struct AaLookupSmallboneCell {
    pub num_used: u16,
    pub padding: u16,
    pub payload: [u8; 8],
}
impl AaLookupSmallboneCell {
    pub(crate) fn set_entry(&mut self, index: usize, value: u16) {
        self.payload[index * 2..index * 2 + 2].copy_from_slice(&value.to_ne_bytes());
    }
    pub(crate) fn overflow_cursor(&self) -> i32 {
        i32::from_ne_bytes(self.payload[..4].try_into().unwrap())
    }
    pub(crate) fn set_overflow_cursor(&mut self, value: i32) {
        self.payload[..4].copy_from_slice(&value.to_ne_bytes());
    }
}
#[derive(Debug, Default)]
pub struct BlastAaLookupTable {
    pub threshold: i32,
    pub mask: i32,
    pub charsize: i32,
    pub word_length: i32,
    pub lut_word_length: i32,
    pub alphabet_size: i32,
    pub backbone_size: i32,
    pub longest_chain: i32,
    pub thin_backbone: ThinBackbone,
    pub bone_type: EBoneType,
    pub thick_backbone: Vec<AaLookupBackboneCell>,
    pub small_backbone: Vec<AaLookupSmallboneCell>,
    pub overflow: Vec<i32>,
    pub small_overflow: Vec<u16>,
    pub overflow_size: i32,
    pub pv: Vec<u32>,
    pub use_pssm: bool,
    pub scansub_callback: Option<AaScanSubject>,
    pub neighbor_matches: i32,
    pub exact_matches: i32,
}
pub struct NeighborInfo<'a, 'b> {
    pub lookup: &'a mut BlastAaLookupTable,
    pub query_word: &'b [u8],
    pub subject_word: &'a mut [u8],
    pub alphabet_size: i32,
    pub wordsize: i32,
    pub charsize: i32,
    pub matrix: &'b [Vec<i32>],
    pub row_max: &'b [i32],
    pub offset_list: &'b [i32],
    pub threshold: i32,
    pub query_bias: i32,
}
/// Only the two fields read from LookupTableOptions are explicit inputs.
pub fn blast_aa_lookup_table_new(word_size: i32, threshold: f64) -> BlastAaLookupTable {
    assert!(
        word_size > 0 && word_size <= 6,
        "ordinary AA lookup index must fit a signed 32-bit word"
    );
    let mut lookup = BlastAaLookupTable::default();
    lookup.charsize = ilog2(BLASTAA_SIZE as i64) + 1;
    lookup.word_length = word_size;
    for i in 0..lookup.word_length {
        lookup.backbone_size |= ((BLASTAA_SIZE - 1) as i32) << (i * lookup.charsize);
    }
    lookup.backbone_size += 1;
    lookup.mask = (1 << (word_size * lookup.charsize)) - 1;
    lookup.alphabet_size = BLASTAA_SIZE as i32;
    lookup.threshold = threshold as i32;
    lookup.thin_backbone = vec![None; lookup.backbone_size as usize];
    lookup
}
pub fn blast_aa_lookup_finalize(lookup: &mut BlastAaLookupTable, bone_type: EBoneType) -> i32 {
    let mut overflow_cells_needed = 0;
    let mut overflow_cursor = 0;
    let mut longest_chain = 0;
    for i in 0..lookup.backbone_size as usize {
        if let Some(chain) = &lookup.thin_backbone[i] {
            if chain[1] > AA_HITS_PER_CELL as i32 {
                overflow_cells_needed += chain[1];
            }
            if chain[1] > longest_chain {
                longest_chain = chain[1];
            }
        }
    }
    lookup.overflow_size = overflow_cells_needed;
    lookup.longest_chain = longest_chain;
    lookup.bone_type = bone_type;
    if bone_type == EBoneType::Backbone {
        lookup.thick_backbone =
            vec![AaLookupBackboneCell::default(); lookup.backbone_size as usize];
        lookup.pv = vec![0; (lookup.backbone_size as usize >> 5) + 1];
        if overflow_cells_needed > 0 {
            lookup.overflow = vec![0; overflow_cells_needed as usize];
        }
        for i in 0..lookup.backbone_size as usize {
            if let Some(chain) = lookup.thin_backbone[i].take() {
                lookup.pv[i >> 5] |= 1u32 << (i & 31);
                lookup.thick_backbone[i].num_used = chain[1];
                if chain[1] <= AA_HITS_PER_CELL as i32 {
                    for j in 0..chain[1] as usize {
                        lookup.thick_backbone[i].payload[j] = chain[j + 2];
                    }
                } else {
                    lookup.thick_backbone[i].payload[0] = overflow_cursor;
                    let dest = overflow_cursor as usize;
                    overflow_cursor += chain[1];
                    for j in 0..chain[1] as usize {
                        lookup.overflow[dest + j] = chain[j + 2];
                    }
                }
            } else {
                lookup.thick_backbone[i].num_used = 0;
            }
        }
    } else {
        lookup.small_backbone =
            vec![AaLookupSmallboneCell::default(); lookup.backbone_size as usize];
        lookup.pv = vec![0; (lookup.backbone_size as usize >> 5) + 1];
        if overflow_cells_needed > 0 {
            lookup.small_overflow = vec![0; overflow_cells_needed as usize];
        }
        for i in 0..lookup.backbone_size as usize {
            if let Some(chain) = lookup.thin_backbone[i].take() {
                lookup.pv[i >> 5] |= 1u32 << (i & 31);
                lookup.small_backbone[i].num_used = chain[1] as u16;
                if chain[1] <= AA_HITS_PER_CELL as i32 {
                    for j in 0..chain[1] as usize {
                        lookup.small_backbone[i].set_entry(j, chain[j + 2] as u16);
                    }
                } else {
                    lookup.small_backbone[i].set_overflow_cursor(overflow_cursor);
                    let dest = overflow_cursor as usize;
                    overflow_cursor += chain[1];
                    for j in 0..chain[1] as usize {
                        lookup.small_overflow[dest + j] = chain[j + 2] as u16;
                    }
                }
            } else {
                lookup.small_backbone[i].num_used = 0;
            }
        }
    }
    lookup.thin_backbone = Vec::new();
    0
}
pub fn blast_aa_lookup_index_query(
    lookup: &mut BlastAaLookupTable,
    matrix: &[Vec<i32>],
    query: &[u8],
    location: &[SSeqRange],
    query_bias: i32,
) {
    assert!(
        !lookup.use_pssm,
        "PSSM/RPS lookup is outside the standard tblastn path"
    );
    s_add_neighboring_words(lookup, matrix, query, query_bias, location);
}
pub fn s_add_neighboring_words(
    lookup: &mut BlastAaLookupTable,
    matrix: &[Vec<i32>],
    query: &[u8],
    query_bias: i32,
    location: &[SSeqRange],
) {
    assert!(lookup.alphabet_size <= BLASTAA_SIZE as i32);
    let mut row_max = [0i32; BLASTAA_SIZE];
    for i in 0..lookup.alphabet_size as usize {
        row_max[i] = matrix[i][0];
        for j in 1..lookup.alphabet_size as usize {
            row_max[i] = row_max[i].max(matrix[i][j]);
        }
    }
    let mut exact_backbone = vec![None; lookup.backbone_size as usize];
    blast_lookup_index_query_exact_matches(
        &mut exact_backbone,
        lookup.word_length,
        lookup.charsize,
        lookup.word_length,
        query,
        location,
    );
    for i in 0..lookup.backbone_size as usize {
        if let Some(offsets) = exact_backbone[i].take() {
            s_add_word_hits(lookup, matrix, query, &offsets, query_bias, &row_max);
        }
    }
}
pub fn s_add_word_hits(
    lookup: &mut BlastAaLookupTable,
    matrix: &[Vec<i32>],
    query: &[u8],
    offset_list: &[i32],
    query_bias: i32,
    row_max: &[i32],
) {
    let w = &query[offset_list[2] as usize..];
    let mut s = [0u8; 32];
    let mut score = matrix[w[0] as usize][w[0] as usize];
    for i in 1..lookup.word_length as usize {
        score += matrix[w[i] as usize][w[i] as usize];
    }
    if lookup.threshold == 0 || score < lookup.threshold {
        for i in 0..offset_list[1] as usize {
            blast_lookup_add_word_hit(
                &mut lookup.thin_backbone,
                lookup.word_length,
                lookup.charsize,
                w,
                query_bias + offset_list[i + 2],
            );
        }
    }
    if lookup.threshold == 0 {
        return;
    }
    let alphabet_size = lookup.alphabet_size;
    let wordsize = lookup.word_length;
    let charsize = lookup.charsize;
    let threshold = lookup.threshold;
    score = row_max[w[0] as usize];
    for i in 1..wordsize as usize {
        score += row_max[w[i] as usize];
    }
    let mut info = NeighborInfo {
        lookup,
        query_word: w,
        subject_word: &mut s,
        alphabet_size,
        wordsize,
        charsize,
        matrix,
        row_max,
        offset_list,
        threshold,
        query_bias,
    };
    s_add_word_hits_core(&mut info, score, 0);
}
pub fn s_add_word_hits_core(info: &mut NeighborInfo<'_, '_>, mut score: i32, current_pos: i32) {
    let alphabet_size = info.alphabet_size;
    let threshold = info.threshold;
    let query_word = info.query_word;
    let matrix = info.matrix;
    score -= info.row_max[query_word[current_pos as usize] as usize];
    let row = &matrix[query_word[current_pos as usize] as usize];
    if current_pos == info.wordsize - 1 {
        let offset_list = info.offset_list;
        let query_bias = info.query_bias;
        let wordsize = info.wordsize;
        let charsize = info.charsize;
        for i in 0..alphabet_size as usize {
            if score + row[i] >= threshold {
                info.subject_word[current_pos as usize] = i as u8;
                for j in 0..offset_list[1] as usize {
                    blast_lookup_add_word_hit(
                        &mut info.lookup.thin_backbone,
                        wordsize,
                        charsize,
                        info.subject_word,
                        query_bias + offset_list[j + 2],
                    );
                }
            }
        }
        return;
    }
    for i in 0..alphabet_size as usize {
        if score + row[i] >= threshold {
            info.subject_word[current_pos as usize] = i as u8;
            s_add_word_hits_core(info, score + row[i], current_pos + 1);
        }
    }
}
