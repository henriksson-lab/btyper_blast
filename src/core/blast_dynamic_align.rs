// NCBI BLAST+ 2.17.0+, selected blast_gapalign.c/h and gapinfo.h definitions.
// Authors: NCBI BLAST developers. See LICENSE-NCBI.
use super::blast_extend::BlastInitHSP;
use super::blast_parameters::BlastScoringParameters;
use super::blast_util::FENCE_SENTRY;
use super::gapinfo::{gap_prelim_edit_block_add, EGapAlignOpType, GapPrelimEditBlock};

pub const MININT: i32 = i32::MIN / 2;
pub const CHUNKSIZE: i32 = 2_097_152;
pub const SCRIPT_SUB: u8 = EGapAlignOpType::Sub as u8;
pub const SCRIPT_GAP_IN_A: u8 = EGapAlignOpType::Del as u8;
pub const SCRIPT_GAP_IN_B: u8 = EGapAlignOpType::Ins as u8;
pub const SCRIPT_OP_MASK: u8 = 0x07;
pub const SCRIPT_EXTEND_GAP_A: u8 = 0x10;
pub const SCRIPT_EXTEND_GAP_B: u8 = 0x40;

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
#[repr(C)]
pub struct BlastGapDP {
    pub best: i32,
    pub best_gap: i32,
}
#[derive(Debug)]
pub struct GapStateArrayStruct {
    pub length: i32,
    pub used: i32,
    pub state_array: Vec<u8>,
}
/// Sequence pointer binding with a signed origin. This preserves original
/// sentinel/anchor reads at negative relative offsets without copying data or
/// constructing an out-of-allocation address.
#[derive(Clone, Copy, Debug)]
pub struct DynamicSequence<'a> {
    pub data: &'a [u8],
    pub offset: i32,
}
impl DynamicSequence<'_> {
    pub(crate) fn at(&self, index: i32) -> u8 {
        self.data[(self.offset + index) as usize]
    }
    pub(crate) fn shifted(&self, amount: i32) -> Self {
        Self {
            data: self.data,
            offset: self.offset + amount,
        }
    }
}
#[derive(Debug, Default)]
pub struct BlastDynamicAlign<'a> {
    pub matrix: &'a [Vec<i32>],
    pub pssm: Option<&'a [Vec<i32>]>,
    pub position_based: bool,
    pub gap_x_dropoff: i32,
    pub dp_mem: Vec<BlastGapDP>,
    pub dp_mem_alloc: i32,
    pub state_struct: Vec<GapStateArrayStruct>,
    pub query_start: i32,
    pub query_stop: i32,
    pub subject_start: i32,
    pub subject_stop: i32,
    pub score: i32,
}
#[derive(Clone, Copy)]
struct EditScriptRow {
    block: usize,
    start: usize,
}

/// The returned stable arena index replaces the original linked-node pointer.
pub fn s_gap_get_state(head: &mut Vec<GapStateArrayStruct>, mut length: i32) -> usize {
    assert!(length >= 0);
    let chunksize = CHUNKSIZE.max(length + length / 3);
    length += length / 3;
    if head.is_empty() {
        head.push(GapStateArrayStruct {
            length: chunksize,
            used: 0,
            state_array: vec![0; chunksize as usize],
        });
        return 0;
    }
    for index in 0..head.len() {
        if length < head[index].length - head[index].used {
            return index;
        }
        if head[index].used == 0 {
            head[index].state_array = vec![0; chunksize as usize];
            head[index].length = chunksize;
            return index;
        }
    }
    head.push(GapStateArrayStruct {
        length: chunksize,
        used: 0,
        state_array: vec![0; chunksize as usize],
    });
    head.len() - 1
}
pub fn s_gap_purge_state(state_struct: &mut [GapStateArrayStruct]) -> bool {
    for state in state_struct {
        state.used = 0;
    }
    true
}
/// Original traceback DP kernel. Pointer origins retain the C convention:
/// forward A starts at relative1, forward B may read the end sentinel, and
/// reverse B may read relative-1. Supply the original surrounding buffers.
pub fn align_ex(
    a: &DynamicSequence<'_>,
    b: &DynamicSequence<'_>,
    m: i32,
    n: i32,
    a_offset: &mut i32,
    b_offset: &mut i32,
    edit_block: &mut GapPrelimEditBlock,
    gap_align: &mut BlastDynamicAlign<'_>,
    score_params: &BlastScoringParameters,
    query_offset: i32,
    reversed: bool,
    reverse_sequence: bool,
    mut fence_hit: Option<&mut bool>,
) -> i32 {
    let matrix = gap_align.matrix;
    let pssm = if gap_align.position_based {
        Some(gap_align.pssm.expect("position-based matrix required"))
    } else {
        None
    };
    *a_offset = 0;
    *b_offset = 0;
    let gap_open = score_params.gap_open;
    let gap_extend = score_params.gap_extend;
    let gap_open_extend = gap_open + gap_extend;
    let x_dropoff = gap_align.gap_x_dropoff.max(gap_open_extend);
    if n <= 0 || m <= 0 {
        return 0;
    }
    s_gap_purge_state(&mut gap_align.state_struct);
    let mut edit_script_num_rows = 100;
    let mut edit_script = vec![EditScriptRow { block: 0, start: 0 }; edit_script_num_rows];
    let mut edit_start_offset = vec![0; edit_script_num_rows];
    let num_extra_cells = if gap_extend > 0 {
        x_dropoff / gap_extend + 3
    } else {
        n + 3
    };
    if num_extra_cells > gap_align.dp_mem_alloc {
        gap_align.dp_mem_alloc = (num_extra_cells + 100).max(2 * gap_align.dp_mem_alloc);
        gap_align.dp_mem = vec![BlastGapDP::default(); gap_align.dp_mem_alloc as usize];
    }
    let first_state = s_gap_get_state(&mut gap_align.state_struct, num_extra_cells);
    edit_script[0] = EditScriptRow {
        block: first_state,
        start: 0,
    };
    edit_start_offset[0] = 0;
    let mut score = -gap_open_extend;
    gap_align.dp_mem[0] = BlastGapDP {
        best: 0,
        best_gap: -gap_open_extend,
    };
    let mut i = 1;
    while i <= n {
        if score < -x_dropoff {
            break;
        }
        gap_align.dp_mem[i as usize] = BlastGapDP {
            best: score,
            best_gap: score - gap_open_extend,
        };
        score -= gap_extend;
        gap_align.state_struct[first_state].state_array[i as usize] = SCRIPT_GAP_IN_A;
        i += 1;
    }
    gap_align.state_struct[first_state].used = i + 1;
    let mut b_size = i;
    let mut best_score = 0;
    let mut first_b_index = 0;
    let b_increment = if reverse_sequence { -1 } else { 1 };
    for a_index in 1..=m {
        let state = s_gap_get_state(
            &mut gap_align.state_struct,
            if gap_extend > 0 {
                b_size - first_b_index + num_extra_cells
            } else {
                n + 3 - first_b_index
            },
        );
        if a_index as usize == edit_script_num_rows {
            edit_script_num_rows *= 2;
            edit_script.resize(edit_script_num_rows, EditScriptRow { block: 0, start: 0 });
            edit_start_offset.resize(edit_script_num_rows, 0);
        }
        edit_script[a_index as usize] = EditScriptRow {
            block: state,
            start: (gap_align.state_struct[state].used + 1) as usize,
        };
        edit_start_offset[a_index as usize] = first_b_index;
        let row = edit_script[a_index as usize];
        let orig_b_index = first_b_index;
        let matrix_row = if let Some(pssm) = pssm {
            if reversed || reverse_sequence {
                &pssm[(m - a_index) as usize]
            } else {
                &pssm[(a_index + query_offset) as usize]
            }
        } else if reverse_sequence {
            &matrix[a.at(m - a_index) as usize]
        } else {
            &matrix[a.at(a_index) as usize]
        };
        let mut b_ptr = if reverse_sequence {
            n - first_b_index
        } else {
            first_b_index
        };
        score = MININT;
        let mut score_gap_row = MININT;
        let mut last_b_index = first_b_index;
        let mut b_index = first_b_index;
        while b_index < b_size {
            b_ptr += b_increment;
            let mut score_gap_col = gap_align.dp_mem[b_index as usize].best_gap;
            if b.at(b_ptr) == FENCE_SENTRY {
                if let Some(fence) = fence_hit.as_deref_mut() {
                    *fence = true;
                }
                break;
            }
            let next_score =
                gap_align.dp_mem[b_index as usize].best + matrix_row[b.at(b_ptr) as usize];
            let mut script = SCRIPT_SUB;
            let script_col = SCRIPT_EXTEND_GAP_B;
            let script_row = SCRIPT_EXTEND_GAP_A;
            if score < score_gap_col {
                script = SCRIPT_GAP_IN_B;
                score = score_gap_col;
            }
            if score < score_gap_row {
                script = SCRIPT_GAP_IN_A;
                score = score_gap_row;
            }
            if best_score - score > x_dropoff {
                if first_b_index == b_index {
                    first_b_index += 1;
                } else {
                    gap_align.dp_mem[b_index as usize].best = MININT;
                }
            } else {
                last_b_index = b_index;
                if score > best_score {
                    best_score = score;
                    *a_offset = a_index;
                    *b_offset = b_index;
                }
                score_gap_row -= gap_extend;
                score_gap_col -= gap_extend;
                if score_gap_col < score - gap_open_extend {
                    gap_align.dp_mem[b_index as usize].best_gap = score - gap_open_extend;
                } else {
                    gap_align.dp_mem[b_index as usize].best_gap = score_gap_col;
                    script += script_col;
                }
                if score_gap_row < score - gap_open_extend {
                    score_gap_row = score - gap_open_extend;
                } else {
                    script += script_row;
                }
                gap_align.dp_mem[b_index as usize].best = score;
            }
            score = next_score;
            gap_align.state_struct[row.block].state_array
                [row.start + (b_index - orig_b_index) as usize] = script;
            b_index += 1;
        }
        if first_b_index == b_size || fence_hit.as_deref().copied() == Some(true) {
            break;
        }
        if last_b_index + num_extra_cells + 3 >= gap_align.dp_mem_alloc {
            gap_align.dp_mem_alloc =
                (last_b_index + num_extra_cells + 100).max(2 * gap_align.dp_mem_alloc);
            gap_align
                .dp_mem
                .resize(gap_align.dp_mem_alloc as usize, BlastGapDP::default());
        }
        if last_b_index < b_size - 1 {
            b_size = last_b_index + 1;
        } else {
            while score_gap_row >= best_score - x_dropoff && b_size <= n {
                gap_align.dp_mem[b_size as usize] = BlastGapDP {
                    best: score_gap_row,
                    best_gap: score_gap_row - gap_open_extend,
                };
                score_gap_row -= gap_extend;
                gap_align.state_struct[row.block].state_array
                    [row.start + (b_size - orig_b_index) as usize] = SCRIPT_GAP_IN_A;
                b_size += 1;
            }
        }
        gap_align.state_struct[state].used += b_index.max(b_size) - orig_b_index + 1;
        if b_size <= n {
            gap_align.dp_mem[b_size as usize] = BlastGapDP {
                best: MININT,
                best_gap: MININT,
            };
            b_size += 1;
        }
    }
    let mut a_index = *a_offset;
    let mut b_index = *b_offset;
    let mut script = SCRIPT_SUB;
    if fence_hit.as_deref().copied() != Some(true) {
        while a_index > 0 || b_index > 0 {
            let row = edit_script[a_index as usize];
            let next_script = gap_align.state_struct[row.block].state_array
                [row.start + (b_index - edit_start_offset[a_index as usize]) as usize];
            script = match script {
                SCRIPT_GAP_IN_A => {
                    if next_script & SCRIPT_EXTEND_GAP_A != 0 {
                        SCRIPT_GAP_IN_A
                    } else {
                        next_script & SCRIPT_OP_MASK
                    }
                }
                SCRIPT_GAP_IN_B => {
                    if next_script & SCRIPT_EXTEND_GAP_B != 0 {
                        SCRIPT_GAP_IN_B
                    } else {
                        next_script & SCRIPT_OP_MASK
                    }
                }
                _ => next_script & SCRIPT_OP_MASK,
            };
            if script == SCRIPT_GAP_IN_A {
                b_index -= 1;
            } else if script == SCRIPT_GAP_IN_B {
                a_index -= 1;
            } else {
                a_index -= 1;
                b_index -= 1;
            }
            let op = match script {
                SCRIPT_GAP_IN_A => EGapAlignOpType::Del,
                SCRIPT_GAP_IN_B => EGapAlignOpType::Ins,
                SCRIPT_SUB => EGapAlignOpType::Sub,
                _ => panic!("invalid edit operation"),
            };
            gap_prelim_edit_block_add(edit_block, op, 1);
        }
    }
    best_score
}
/// Original generic matrix/PSSM DP entry point. With score_only=false this
/// delegates to ALIGN_EX; score-only inputs must be initialized matrix codes.
pub fn blast_semi_gapped_align(
    a: &DynamicSequence<'_>,
    b: &DynamicSequence<'_>,
    m: i32,
    n: i32,
    a_offset: &mut i32,
    b_offset: &mut i32,
    score_only: bool,
    edit_block: Option<&mut GapPrelimEditBlock>,
    gap_align: &mut BlastDynamicAlign<'_>,
    score_params: &BlastScoringParameters,
    query_offset: i32,
    reversed: bool,
    reverse_sequence: bool,
    fence_hit: Option<&mut bool>,
) -> i32 {
    if !score_only {
        return align_ex(
            a,
            b,
            m,
            n,
            a_offset,
            b_offset,
            edit_block.expect("traceback output required"),
            gap_align,
            score_params,
            query_offset,
            reversed,
            reverse_sequence,
            fence_hit,
        );
    }
    let matrix = gap_align.matrix;
    let pssm = if gap_align.position_based {
        Some(gap_align.pssm.expect("position-based matrix required"))
    } else {
        None
    };
    *a_offset = 0;
    *b_offset = 0;
    let gap_open = score_params.gap_open;
    let gap_extend = score_params.gap_extend;
    let gap_open_extend = gap_open + gap_extend;
    let x_dropoff = gap_align.gap_x_dropoff.max(gap_open_extend);
    if n <= 0 || m <= 0 {
        return 0;
    }
    let num_extra_cells = if gap_extend > 0 {
        x_dropoff / gap_extend + 3
    } else {
        n + 3
    };
    if num_extra_cells > gap_align.dp_mem_alloc {
        gap_align.dp_mem_alloc = (num_extra_cells + 100).max(2 * gap_align.dp_mem_alloc);
        gap_align.dp_mem = vec![BlastGapDP::default(); gap_align.dp_mem_alloc as usize];
    }
    let mut score = -gap_open_extend;
    gap_align.dp_mem[0] = BlastGapDP {
        best: 0,
        best_gap: -gap_open_extend,
    };
    let mut i = 1;
    while i <= n {
        if score < -x_dropoff {
            break;
        }
        gap_align.dp_mem[i as usize] = BlastGapDP {
            best: score,
            best_gap: score - gap_open_extend,
        };
        score -= gap_extend;

        i += 1;
    }

    let mut b_size = i;
    let mut best_score = 0;
    let mut first_b_index = 0;
    let b_increment = if reverse_sequence { -1 } else { 1 };
    for a_index in 1..=m {
        let matrix_row = if let Some(pssm) = pssm {
            if reversed || reverse_sequence {
                &pssm[(m - a_index) as usize]
            } else {
                &pssm[(a_index + query_offset) as usize]
            }
        } else if reverse_sequence {
            &matrix[a.at(m - a_index) as usize]
        } else {
            &matrix[a.at(a_index) as usize]
        };
        let mut b_ptr = if reverse_sequence {
            n - first_b_index
        } else {
            first_b_index
        };
        score = MININT;
        let mut score_gap_row = MININT;
        let mut last_b_index = first_b_index;
        let mut b_index = first_b_index;
        while b_index < b_size {
            b_ptr += b_increment;
            let mut score_gap_col = gap_align.dp_mem[b_index as usize].best_gap;

            let next_score =
                gap_align.dp_mem[b_index as usize].best + matrix_row[b.at(b_ptr) as usize];

            if score < score_gap_col {
                score = score_gap_col;
            }
            if score < score_gap_row {
                score = score_gap_row;
            }
            if best_score - score > x_dropoff {
                if first_b_index == b_index {
                    first_b_index += 1;
                } else {
                    gap_align.dp_mem[b_index as usize].best = MININT;
                }
            } else {
                last_b_index = b_index;
                if score > best_score {
                    best_score = score;
                    *a_offset = a_index;
                    *b_offset = b_index;
                }
                score_gap_row -= gap_extend;
                score_gap_col -= gap_extend;
                gap_align.dp_mem[b_index as usize].best_gap =
                    (score - gap_open_extend).max(score_gap_col);
                score_gap_row = (score - gap_open_extend).max(score_gap_row);
                gap_align.dp_mem[b_index as usize].best = score;
            }
            score = next_score;

            b_index += 1;
        }
        if first_b_index == b_size {
            break;
        }
        if last_b_index + num_extra_cells + 3 >= gap_align.dp_mem_alloc {
            gap_align.dp_mem_alloc =
                (last_b_index + num_extra_cells + 100).max(2 * gap_align.dp_mem_alloc);
            gap_align
                .dp_mem
                .resize(gap_align.dp_mem_alloc as usize, BlastGapDP::default());
        }
        if last_b_index < b_size - 1 {
            b_size = last_b_index + 1;
        } else {
            while score_gap_row >= best_score - x_dropoff && b_size <= n {
                gap_align.dp_mem[b_size as usize] = BlastGapDP {
                    best: score_gap_row,
                    best_gap: score_gap_row - gap_open_extend,
                };
                score_gap_row -= gap_extend;

                b_size += 1;
            }
        }

        if b_size <= n {
            gap_align.dp_mem[b_size as usize] = BlastGapDP {
                best: MININT,
                best_gap: MININT,
            };
            b_size += 1;
        }
    }
    best_score
}

/// Packed preliminary DP. Reverse A is an aligned packed prefix; forward A
/// points to the byte before the first packed byte. B retains query sentinels.
pub fn s_blast_align_packed_nucl(
    b: &DynamicSequence<'_>,
    a: &DynamicSequence<'_>,
    n: i32,
    m: i32,
    b_offset: &mut i32,
    a_offset: &mut i32,
    gap_align: &mut BlastDynamicAlign<'_>,
    score_params: &BlastScoringParameters,
    reverse_sequence: bool,
    x_dropoff: i32,
) -> i32 {
    let matrix = gap_align.matrix;
    *a_offset = 0;
    *b_offset = 0;
    let gap_open = score_params.gap_open;
    let gap_extend = score_params.gap_extend;
    let gap_open_extend = gap_open + gap_extend;
    let x_dropoff = x_dropoff.max(gap_open_extend);
    if n <= 0 || m <= 0 {
        return 0;
    }
    let num_extra_cells = if gap_extend > 0 {
        x_dropoff / gap_extend + 3
    } else {
        n + 3
    };
    if num_extra_cells > gap_align.dp_mem_alloc {
        gap_align.dp_mem_alloc = (num_extra_cells + 100).max(2 * gap_align.dp_mem_alloc);
        gap_align.dp_mem = vec![BlastGapDP::default(); gap_align.dp_mem_alloc as usize];
    }
    let mut score = -gap_open_extend;
    gap_align.dp_mem[0] = BlastGapDP {
        best: 0,
        best_gap: -gap_open_extend,
    };
    let mut i = 1;
    while i <= n {
        if score < -x_dropoff {
            break;
        }
        gap_align.dp_mem[i as usize] = BlastGapDP {
            best: score,
            best_gap: score - gap_open_extend,
        };
        score -= gap_extend;

        i += 1;
    }

    let mut b_size = i;
    let mut best_score = 0;
    let mut first_b_index = 0;
    let b_increment = if reverse_sequence { -1 } else { 1 };
    for a_index in 1..=m {
        let a_base_pair = if reverse_sequence {
            (a.at((m - a_index) / 4) >> (2 * ((a_index - 1) % 4))) & 3
        } else {
            (a.at(1 + (a_index - 1) / 4) >> (2 * (3 - (a_index - 1) % 4))) & 3
        };
        let matrix_row = &matrix[a_base_pair as usize];
        let mut b_ptr = if reverse_sequence {
            n - first_b_index
        } else {
            first_b_index
        };
        score = MININT;
        let mut score_gap_row = MININT;
        let mut last_b_index = first_b_index;
        let mut b_index = first_b_index;
        while b_index < b_size {
            b_ptr += b_increment;
            let mut score_gap_col = gap_align.dp_mem[b_index as usize].best_gap;

            let next_score =
                gap_align.dp_mem[b_index as usize].best + matrix_row[b.at(b_ptr) as usize];

            if score < score_gap_col {
                score = score_gap_col;
            }
            if score < score_gap_row {
                score = score_gap_row;
            }
            if best_score - score > x_dropoff {
                if first_b_index == b_index {
                    first_b_index += 1;
                } else {
                    gap_align.dp_mem[b_index as usize].best = MININT;
                }
            } else {
                last_b_index = b_index;
                if score > best_score {
                    best_score = score;
                    *a_offset = a_index;
                    *b_offset = b_index;
                }
                score_gap_row -= gap_extend;
                score_gap_col -= gap_extend;
                gap_align.dp_mem[b_index as usize].best_gap =
                    (score - gap_open_extend).max(score_gap_col);
                score_gap_row = (score - gap_open_extend).max(score_gap_row);
                gap_align.dp_mem[b_index as usize].best = score;
            }
            score = next_score;

            b_index += 1;
        }
        if first_b_index == b_size {
            break;
        }
        if last_b_index + num_extra_cells + 3 >= gap_align.dp_mem_alloc {
            gap_align.dp_mem_alloc =
                (last_b_index + num_extra_cells + 100).max(2 * gap_align.dp_mem_alloc);
            gap_align
                .dp_mem
                .resize(gap_align.dp_mem_alloc as usize, BlastGapDP::default());
        }
        if last_b_index < b_size - 1 {
            b_size = last_b_index + 1;
        } else {
            while score_gap_row >= best_score - x_dropoff && b_size <= n {
                gap_align.dp_mem[b_size as usize] = BlastGapDP {
                    best: score_gap_row,
                    best_gap: score_gap_row - gap_open_extend,
                };
                score_gap_row -= gap_extend;

                b_size += 1;
            }
        }

        if b_size <= n {
            gap_align.dp_mem[b_size as usize] = BlastGapDP {
                best: MININT,
                best_gap: MININT,
            };
            b_size += 1;
        }
    }
    best_score
}

/// Original packed-nucleotide preliminary extension around one initial HSP.
/// Query and subject origins correspond to their sequence-block pointers.
pub fn s_blast_dyn_prog_nt_gapped_alignment(
    query: &DynamicSequence<'_>,
    query_length: i32,
    subject: &DynamicSequence<'_>,
    subject_length: i32,
    gap_align: &mut BlastDynamicAlign<'_>,
    score_params: &BlastScoringParameters,
    init_hsp: &BlastInitHSP,
) -> i16 {
    let mut x_dropoff = gap_align.gap_x_dropoff;
    if init_hsp
        .ungapped_data
        .as_ref()
        .expect("ungapped initial HSP required")
        .score
        < x_dropoff
    {
        x_dropoff = init_hsp.ungapped_data.as_ref().unwrap().score;
    }
    let offset_adjustment = 4 - init_hsp.offsets.s_off % 4;
    let mut q_length = init_hsp.offsets.q_off + offset_adjustment;
    let mut s_length = init_hsp.offsets.s_off + offset_adjustment;
    if q_length > query_length || s_length > subject_length {
        q_length -= 4;
        s_length -= 4;
    }
    let mut private_q_start = 0;
    let mut private_s_start = 0;
    let score_left = s_blast_align_packed_nucl(
        query,
        subject,
        q_length,
        s_length,
        &mut private_q_start,
        &mut private_s_start,
        gap_align,
        score_params,
        true,
        x_dropoff,
    );
    if score_left < 0 {
        return -1;
    }
    gap_align.query_start = q_length - private_q_start;
    gap_align.subject_start = s_length - private_s_start;
    let mut score_right = 0;
    if q_length < query_length && s_length < subject_length {
        let mut q_stop = gap_align.query_stop;
        let mut s_stop = gap_align.subject_stop;
        score_right = s_blast_align_packed_nucl(
            &query.shifted(q_length - 1),
            &subject.shifted((s_length + 3) / 4 - 1),
            query_length - q_length,
            subject_length - s_length,
            &mut q_stop,
            &mut s_stop,
            gap_align,
            score_params,
            false,
            x_dropoff,
        );
        gap_align.query_stop = q_stop;
        gap_align.subject_stop = s_stop;
        if score_right < 0 {
            return -1;
        }
        gap_align.query_stop += q_length;
        gap_align.subject_stop += s_length;
    } else {
        gap_align.query_stop = q_length;
        gap_align.subject_stop = s_length;
    }
    gap_align.score = score_right + score_left;
    0
}
