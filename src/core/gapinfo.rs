// NCBI BLAST+ 2.17.0+, selected functions from gapinfo.c and gapinfo.h.
// Author: Ilya Dondoshansky. See LICENSE-NCBI.

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[repr(i32)]
pub enum EGapAlignOpType {
    Del = 0,
    Del2 = 1,
    Del1 = 2,
    Sub = 3,
    Ins1 = 4,
    Ins2 = 5,
    Ins = 6,
    Decline = 7,
    Invalid = 8,
}

#[derive(Debug)]
pub struct GapEditScript {
    pub op_type: Vec<EGapAlignOpType>,
    pub num: Vec<i32>,
    pub size: i32,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct GapPrelimEditScript {
    pub op_type: EGapAlignOpType,
    pub num: i32,
}

#[derive(Debug)]
pub struct GapPrelimEditBlock {
    pub edit_ops: Vec<GapPrelimEditScript>,
    pub num_ops_allocated: i32,
    pub num_ops: i32,
    pub last_op: EGapAlignOpType,
}

/// Matches C `GapEditScriptNew`; zero initialization retains C's Del opcode.
pub fn gap_edit_script_new(size: i32) -> Option<GapEditScript> {
    if size <= 0 {
        return None;
    }
    Some(GapEditScript {
        op_type: vec![EGapAlignOpType::Del; size as usize],
        num: vec![0; size as usize],
        size,
    })
}

/// Matches C `GapEditScriptDelete`; consuming ownership replaces explicit frees.
pub fn gap_edit_script_delete(old: Option<GapEditScript>) -> Option<GapEditScript> {
    drop(old);
    None
}

/// Matches C `s_GapPrelimEditBlockRealloc` including its equality threshold.
pub fn s_gap_prelim_edit_block_realloc(edit_block: &mut GapPrelimEditBlock, total_ops: i32) -> i16 {
    if edit_block.num_ops_allocated <= total_ops {
        let new_size = total_ops * 2;
        if edit_block
            .edit_ops
            .try_reserve_exact(new_size as usize - edit_block.edit_ops.len())
            .is_err()
        {
            return -1;
        }
        edit_block.edit_ops.resize(
            new_size as usize,
            GapPrelimEditScript {
                op_type: EGapAlignOpType::Del,
                num: 0,
            },
        );
        edit_block.num_ops_allocated = new_size;
    }
    0
}

/// Matches C `s_GapPrelimEditBlockAddNew`.
pub fn s_gap_prelim_edit_block_add_new(
    edit_block: &mut GapPrelimEditBlock,
    op_type: EGapAlignOpType,
    num_ops: i32,
) -> i16 {
    if s_gap_prelim_edit_block_realloc(edit_block, edit_block.num_ops + 2) != 0 {
        return -1;
    }
    assert_ne!(op_type, EGapAlignOpType::Invalid);
    edit_block.last_op = op_type;
    edit_block.edit_ops[edit_block.num_ops as usize].op_type = op_type;
    edit_block.edit_ops[edit_block.num_ops as usize].num = num_ops;
    edit_block.num_ops += 1;
    0
}

/// Matches C `GapPrelimEditBlockAdd`.
pub fn gap_prelim_edit_block_add(
    edit_block: &mut GapPrelimEditBlock,
    op_type: EGapAlignOpType,
    num_ops: i32,
) {
    if num_ops == 0 {
        return;
    }
    if edit_block.last_op == op_type {
        edit_block.edit_ops[edit_block.num_ops as usize - 1].num += num_ops;
    } else {
        s_gap_prelim_edit_block_add_new(edit_block, op_type, num_ops);
    }
}

/// Matches C `GapPrelimEditBlockNew` including the initial 200-slot allocation.
pub fn gap_prelim_edit_block_new() -> GapPrelimEditBlock {
    let mut edit_block = GapPrelimEditBlock {
        edit_ops: Vec::new(),
        num_ops_allocated: 0,
        num_ops: 0,
        last_op: EGapAlignOpType::Invalid,
    };
    s_gap_prelim_edit_block_realloc(&mut edit_block, 100);
    edit_block
}

/// Matches C `GapPrelimEditBlockFree`.
pub fn gap_prelim_edit_block_free(
    edit_block: Option<GapPrelimEditBlock>,
) -> Option<GapPrelimEditBlock> {
    drop(edit_block);
    None
}

/// Matches C `GapPrelimEditBlockReset`; allocated storage is retained.
pub fn gap_prelim_edit_block_reset(edit_block: Option<&mut GapPrelimEditBlock>) {
    if let Some(edit_block) = edit_block {
        edit_block.num_ops = 0;
        edit_block.last_op = EGapAlignOpType::Invalid;
    }
}
