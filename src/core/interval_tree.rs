// NCBI BLAST+ 2.17.0+, blast_itree.c/.h.
// Author: Jason Papadopoulos. See LICENSE-NCBI.
use super::blast_hits::BlastHSP;
use super::blast_query_info::BlastQueryInfo;
use std::ptr::NonNull;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[repr(i32)]
pub enum EIntervalDirection {
    Left = 0,
    Right = 1,
    Neither = 2,
}
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[repr(i32)]
pub enum EITreeIndexMethod {
    QueryOnly = 0,
    QueryAndSubject = 1,
    QueryOnlyStrandIndifferent = 2,
}
#[derive(Clone, Copy, Debug, Default)]
#[repr(C)]
pub struct SIntervalNode {
    pub leftend: i32,
    pub rightend: i32,
    pub leftptr: i32,
    pub midptr: i32,
    pub rightptr: i32,
    pub hsp: Option<NonNull<BlastHSP>>,
}
/// Nodes borrow stable HSP allocations and never own/free them. Keep every
/// indexed HSP alive until reset/free; moving or sorting Box handles is fine.
/// Query/subject coordinates, context and QueryInfo coordinate grouping must
/// stay stable while indexed. Score/metric changes are read live, but no HSP
/// mutation may overlap a tree call. Removing/evicting an owning HSP invalidates
/// the tree unless its lifetime is separately retained until reset/free.
#[derive(Debug)]
pub struct BlastIntervalTree {
    pub nodes: Vec<SIntervalNode>,
    pub num_alloc: i32,
    pub num_used: i32,
    pub s_min: i32,
    pub s_max: i32,
}

pub fn s_interval_node_init(
    tree: &mut BlastIntervalTree,
    parent_index: i32,
    dir: EIntervalDirection,
    ret_status: &mut i16,
) -> i32 {
    *ret_status = 0;
    if tree.num_used == tree.num_alloc {
        tree.num_alloc = tree.num_alloc.wrapping_mul(2);
        let count = tree.num_alloc as usize;
        if count < tree.nodes.len()
            || tree
                .nodes
                .try_reserve_exact(count - tree.nodes.len())
                .is_err()
        {
            *ret_status = 50;
            return 0;
        }
        tree.nodes.resize(count, SIntervalNode::default());
    }
    let new_index = tree.num_used;
    tree.num_used += 1;
    if dir == EIntervalDirection::Neither {
        return new_index;
    }
    let parent = tree.nodes[parent_index as usize];
    let node = &mut tree.nodes[new_index as usize];
    node.leftptr = 0;
    node.midptr = 0;
    node.rightptr = 0;
    node.hsp = None;
    let midpt = (parent.leftend as i64 + parent.rightend as i64) / 2;
    if dir == EIntervalDirection::Left {
        node.leftend = parent.leftend;
        node.rightend = midpt as i32;
    } else {
        node.leftend = midpt as i32 + 1;
        node.rightend = parent.rightend;
    }
    new_index
}
pub fn s_interval_root_node_init(
    tree: &mut BlastIntervalTree,
    region_start: i32,
    region_end: i32,
    retval: &mut i16,
) -> i32 {
    let index = s_interval_node_init(tree, 0, EIntervalDirection::Neither, retval);
    if *retval != 0 {
        return 0;
    }
    let node = &mut tree.nodes[index as usize];
    node.leftptr = 0;
    node.midptr = 0;
    node.rightptr = 0;
    node.hsp = None;
    node.leftend = region_start;
    node.rightend = region_end;
    index
}
pub fn blast_interval_tree_init(
    q_start: i32,
    q_end: i32,
    s_start: i32,
    s_end: i32,
) -> Option<Box<BlastIntervalTree>> {
    let size = 100;
    let mut nodes = Vec::new();
    if nodes.try_reserve_exact(size).is_err() {
        return None;
    }
    nodes.resize(size, SIntervalNode::default());
    let mut tree = Box::new(BlastIntervalTree {
        nodes,
        num_alloc: size as i32,
        num_used: 0,
        s_min: s_start,
        s_max: s_end,
    });
    let mut retval = 0;
    s_interval_root_node_init(&mut tree, q_start, q_end, &mut retval);
    if retval != 0 {
        return None;
    }
    Some(tree)
}
pub fn blast_interval_tree_free(
    tree: Option<Box<BlastIntervalTree>>,
) -> Option<Box<BlastIntervalTree>> {
    drop(tree);
    None
}
pub fn blast_interval_tree_reset(tree: &mut BlastIntervalTree) {
    tree.num_used = 1;
    let root = &mut tree.nodes[0];
    root.leftptr = 0;
    root.midptr = 0;
    root.rightptr = 0;
    root.hsp = None;
}
pub fn s_get_query_strand_offset(query_info: &BlastQueryInfo, context: i32) -> i32 {
    let mut c = context;
    while c != 0 {
        let frame = query_info.contexts[c as usize].frame;
        if frame == 0 || frame.signum() != query_info.contexts[c as usize - 1].frame.signum() {
            break;
        }
        c -= 1;
    }
    query_info.contexts[c as usize].query_offset
}
pub fn s_hs_ps_have_common_endpoint(
    in_hsp: &BlastHSP,
    in_q_start: i32,
    tree_hsp: &BlastHSP,
    tree_q_start: i32,
    which_end: EIntervalDirection,
) -> Option<NonNull<BlastHSP>> {
    if in_q_start != tree_q_start {
        return None;
    }
    if in_hsp.subject.frame.signum() != tree_hsp.subject.frame.signum() {
        return None;
    }
    let matches = if which_end == EIntervalDirection::Left {
        in_hsp.query.offset == tree_hsp.query.offset
            && in_hsp.subject.offset == tree_hsp.subject.offset
    } else {
        in_hsp.query.end == tree_hsp.query.end && in_hsp.subject.end == tree_hsp.subject.end
    };
    if matches {
        if in_hsp.score > tree_hsp.score {
            return Some(NonNull::from(in_hsp));
        }
        if in_hsp.score < tree_hsp.score {
            return Some(NonNull::from(tree_hsp));
        }
        let in_q_length = in_hsp.query.end - in_hsp.query.offset;
        let tree_q_length = tree_hsp.query.end - tree_hsp.query.offset;
        if in_q_length > tree_q_length {
            return Some(NonNull::from(tree_hsp));
        }
        if in_q_length < tree_q_length {
            return Some(NonNull::from(in_hsp));
        }
        let in_s_length = in_hsp.subject.end - in_hsp.subject.offset;
        let tree_s_length = tree_hsp.subject.end - tree_hsp.subject.offset;
        if in_s_length > tree_s_length {
            return Some(NonNull::from(tree_hsp));
        }
        if in_s_length < tree_s_length {
            return Some(NonNull::from(in_hsp));
        }
        return Some(NonNull::from(tree_hsp));
    }
    None
}
/// # Safety
/// All HSP handles reachable in tree must obey BlastIntervalTree's lifetime
/// and coordinate-stability contract.
pub unsafe fn s_midpoint_tree_has_hsp_endpoint(
    tree: &mut BlastIntervalTree,
    root_index: i32,
    in_hsp: &BlastHSP,
    in_q_start: i32,
    which_end: EIntervalDirection,
) -> bool {
    // SAFETY: caller retains all indexed allocations for this traversal.
    unsafe {
        let mut root = root_index as usize;
        let target_offset = if which_end == EIntervalDirection::Left {
            in_hsp.subject.offset
        } else {
            in_hsp.subject.end
        };
        loop {
            assert!(
                target_offset >= tree.nodes[root].leftend
                    && target_offset <= tree.nodes[root].rightend
            );
            let mut tmp = tree.nodes[root].midptr;
            let mut list = root;
            while tmp != 0 {
                let next = tmp as usize;
                let tree_hsp = tree.nodes[next].hsp.unwrap();
                let best = s_hs_ps_have_common_endpoint(
                    in_hsp,
                    in_q_start,
                    tree_hsp.as_ref(),
                    tree.nodes[next].leftptr,
                    which_end,
                );
                tmp = tree.nodes[next].midptr;
                if best == Some(tree_hsp) {
                    return true;
                } else if best == Some(NonNull::from(in_hsp)) {
                    tree.nodes[list].midptr = tmp;
                }
                list = next;
            }
            tmp = 0;
            let midpt = (tree.nodes[root].leftend as i64 + tree.nodes[root].rightend as i64) / 2;
            if (target_offset as i64) < midpt {
                tmp = tree.nodes[root].leftptr;
            } else if (target_offset as i64) > midpt {
                tmp = tree.nodes[root].rightptr;
            }
            if tmp == 0 {
                return false;
            }
            let next = tmp as usize;
            if let Some(tree_hsp) = tree.nodes[next].hsp {
                let best = s_hs_ps_have_common_endpoint(
                    in_hsp,
                    in_q_start,
                    tree_hsp.as_ref(),
                    tree.nodes[next].leftptr,
                    which_end,
                );
                if best == Some(tree_hsp) {
                    return true;
                } else if best == Some(NonNull::from(in_hsp)) {
                    if (target_offset as i64) < midpt {
                        tree.nodes[root].leftptr = 0;
                    } else if (target_offset as i64) > midpt {
                        tree.nodes[root].rightptr = 0;
                    }
                    return false;
                }
                break;
            }
            root = next;
        }
        false
    }
}
/// # Safety
/// All HSP handles reachable in tree must obey its lifetime/stability contract.
pub unsafe fn s_interval_tree_has_hsp_endpoint(
    tree: &mut BlastIntervalTree,
    in_hsp: &BlastHSP,
    in_q_start: i32,
    which_end: EIntervalDirection,
) -> bool {
    // SAFETY: caller retains all indexed HSP allocations.
    unsafe {
        let mut root = 0;
        let target_offset = in_q_start
            + if which_end == EIntervalDirection::Left {
                in_hsp.query.offset
            } else {
                in_hsp.query.end
            };
        loop {
            assert!(
                target_offset >= tree.nodes[root].leftend
                    && target_offset <= tree.nodes[root].rightend
            );
            let mut tmp = tree.nodes[root].midptr;
            if tmp != 0
                && s_midpoint_tree_has_hsp_endpoint(tree, tmp, in_hsp, in_q_start, which_end)
            {
                return true;
            }
            tmp = 0;
            let midpt = (tree.nodes[root].leftend as i64 + tree.nodes[root].rightend as i64) / 2;
            if (target_offset as i64) < midpt {
                tmp = tree.nodes[root].leftptr;
            } else if (target_offset as i64) > midpt {
                tmp = tree.nodes[root].rightptr;
            }
            if tmp == 0 {
                return false;
            }
            let next = tmp as usize;
            if let Some(tree_hsp) = tree.nodes[next].hsp {
                let best = s_hs_ps_have_common_endpoint(
                    in_hsp,
                    in_q_start,
                    tree_hsp.as_ref(),
                    tree.nodes[next].leftptr,
                    which_end,
                );
                if best == Some(tree_hsp) {
                    return true;
                } else if best == Some(NonNull::from(in_hsp)) {
                    if (target_offset as i64) < midpt {
                        tree.nodes[root].leftptr = 0;
                    } else if (target_offset as i64) > midpt {
                        tree.nodes[root].rightptr = 0;
                    }
                    return false;
                }
                break;
            }
            root = next;
        }
        false
    }
}
/// # Safety
/// hsp must point to a live stable BlastHSP allocation, retained until reset or
/// free. All other reachable handles must satisfy the same contract. Coordinate
/// and context keys must remain consistent with QueryInfo while indexed.
pub unsafe fn blast_interval_tree_add_hsp(
    hsp: NonNull<BlastHSP>,
    tree: &mut BlastIntervalTree,
    query_info: &BlastQueryInfo,
    index_method: EITreeIndexMethod,
) -> i16 {
    // SAFETY: supplied/previously inserted HSP allocations remain live.
    unsafe {
        let input = hsp.as_ref();
        let mut query_start = s_get_query_strand_offset(query_info, input.context);
        let (mut region_start, mut region_end);
        if index_method == EITreeIndexMethod::QueryOnlyStrandIndifferent
            && query_info.contexts[input.context as usize].frame == -1
        {
            region_end = query_start - input.query.offset;
            region_start = query_start - input.query.end;
            query_start =
                query_start - query_info.contexts[input.context as usize].query_length - 1;
        } else {
            region_start = query_start + input.query.offset;
            region_end = query_start + input.query.end;
        }
        assert!(region_start >= tree.nodes[0].leftend && region_end <= tree.nodes[0].rightend);
        assert!(input.query.offset <= input.query.end && input.subject.offset <= input.subject.end);
        let mut retval = 0;
        if index_method == EITreeIndexMethod::QueryAndSubject {
            assert!(input.subject.offset >= tree.s_min && input.subject.end <= tree.s_max);
            if s_interval_tree_has_hsp_endpoint(tree, input, query_start, EIntervalDirection::Left)
            {
                return retval;
            }
            if s_interval_tree_has_hsp_endpoint(tree, input, query_start, EIntervalDirection::Right)
            {
                return retval;
            }
        }
        let mut index_subject_range = false;
        let mut root = 0usize;
        let new_index =
            s_interval_node_init(tree, 0, EIntervalDirection::Neither, &mut retval) as usize;
        if retval != 0 {
            return retval;
        }
        tree.nodes[new_index].leftptr = query_start;
        tree.nodes[new_index].midptr = 0;
        tree.nodes[new_index].hsp = Some(hsp);
        loop {
            assert!(
                region_start >= tree.nodes[root].leftend && region_end <= tree.nodes[root].rightend
            );
            let middle = (tree.nodes[root].leftend as i64 + tree.nodes[root].rightend as i64) / 2;
            let old_index;
            let which_half;
            if (region_end as i64) < middle {
                if tree.nodes[root].leftptr == 0 {
                    tree.nodes[root].leftptr = new_index as i32;
                    return retval;
                }
                old_index = tree.nodes[root].leftptr as usize;
                if tree.nodes[old_index].hsp.is_none() {
                    root = old_index;
                    continue;
                } else {
                    which_half = EIntervalDirection::Left;
                }
            } else if (region_start as i64) > middle {
                if tree.nodes[root].rightptr == 0 {
                    tree.nodes[root].rightptr = new_index as i32;
                    return retval;
                }
                old_index = tree.nodes[root].rightptr as usize;
                if tree.nodes[old_index].hsp.is_none() {
                    root = old_index;
                    continue;
                } else {
                    which_half = EIntervalDirection::Right;
                }
            } else {
                if index_subject_range
                    || index_method == EITreeIndexMethod::QueryOnly
                    || index_method == EITreeIndexMethod::QueryOnlyStrandIndifferent
                {
                    tree.nodes[new_index].midptr = tree.nodes[root].midptr;
                    tree.nodes[root].midptr = new_index as i32;
                    return retval;
                } else {
                    index_subject_range = true;
                    if tree.nodes[root].midptr == 0 {
                        let mid_index =
                            s_interval_root_node_init(tree, tree.s_min, tree.s_max, &mut retval);
                        if retval != 0 {
                            return retval;
                        }
                        tree.nodes[root].midptr = mid_index;
                    }
                    root = tree.nodes[root].midptr as usize;
                    region_start = input.subject.offset;
                    region_end = input.subject.end;
                    continue;
                }
            }
            let mid_index =
                s_interval_node_init(tree, root as i32, which_half, &mut retval) as usize;
            if retval != 0 {
                return retval;
            }
            let old_hsp = tree.nodes[old_index].hsp.unwrap().as_ref();
            if which_half == EIntervalDirection::Left {
                tree.nodes[root].leftptr = mid_index as i32;
            } else {
                tree.nodes[root].rightptr = mid_index as i32;
            }
            let (mut old_region_start, mut old_region_end);
            if index_subject_range {
                old_region_start = old_hsp.subject.offset;
                old_region_end = old_hsp.subject.end;
            } else if index_method == EITreeIndexMethod::QueryOnlyStrandIndifferent
                && query_info.contexts[old_hsp.context as usize].frame == -1
            {
                let q_start = s_get_query_strand_offset(query_info, old_hsp.context);
                old_region_end = q_start - old_hsp.query.offset;
                old_region_start = q_start - old_hsp.query.end;
            } else {
                old_region_start = tree.nodes[old_index].leftptr + old_hsp.query.offset;
                old_region_end = tree.nodes[old_index].leftptr + old_hsp.query.end;
            }
            root = mid_index;
            let middle = (tree.nodes[root].leftend as i64 + tree.nodes[root].rightend as i64) / 2;
            if (old_region_end as i64) < middle {
                tree.nodes[mid_index].leftptr = old_index as i32;
            } else if (old_region_start as i64) > middle {
                tree.nodes[mid_index].rightptr = old_index as i32;
            } else if index_subject_range
                || index_method == EITreeIndexMethod::QueryOnly
                || index_method == EITreeIndexMethod::QueryOnlyStrandIndifferent
            {
                tree.nodes[mid_index].midptr = old_index as i32;
            } else {
                let mid_index2 =
                    s_interval_root_node_init(tree, tree.s_min, tree.s_max, &mut retval) as usize;
                if retval != 0 {
                    return retval;
                }
                old_region_start = old_hsp.subject.offset;
                old_region_end = old_hsp.subject.end;
                tree.nodes[mid_index].midptr = mid_index2 as i32;
                let middle = (tree.nodes[mid_index2].leftend as i64
                    + tree.nodes[mid_index2].rightend as i64)
                    / 2;
                if (old_region_end as i64) < middle {
                    tree.nodes[mid_index2].leftptr = old_index as i32;
                } else if (old_region_start as i64) > middle {
                    tree.nodes[mid_index2].rightptr = old_index as i32;
                } else {
                    tree.nodes[mid_index2].midptr = old_index as i32;
                }
            }
        }
    }
}
pub fn s_hsp_is_contained(
    in_hsp: &BlastHSP,
    in_q_start: i32,
    tree_hsp: &BlastHSP,
    tree_q_start: i32,
    min_diag_separation: i32,
) -> bool {
    if in_q_start != tree_q_start {
        return false;
    }
    if in_hsp.score <= tree_hsp.score
        && in_hsp.subject.frame.signum() == tree_hsp.subject.frame.signum()
        && tree_hsp.query.offset <= in_hsp.query.offset
        && tree_hsp.query.end >= in_hsp.query.offset
        && tree_hsp.subject.offset <= in_hsp.subject.offset
        && tree_hsp.subject.end >= in_hsp.subject.offset
        && tree_hsp.query.offset <= in_hsp.query.end
        && tree_hsp.query.end >= in_hsp.query.end
        && tree_hsp.subject.offset <= in_hsp.subject.end
        && tree_hsp.subject.end >= in_hsp.subject.end
    {
        if min_diag_separation == 0 {
            return true;
        }
        let start = tree_hsp
            .query
            .offset
            .wrapping_sub(tree_hsp.subject.offset)
            .wrapping_sub(in_hsp.query.offset.wrapping_sub(in_hsp.subject.offset))
            .wrapping_abs();
        let end = tree_hsp
            .query
            .end
            .wrapping_sub(tree_hsp.subject.end)
            .wrapping_sub(in_hsp.query.end.wrapping_sub(in_hsp.subject.end))
            .wrapping_abs();
        if start < min_diag_separation || end < min_diag_separation {
            return true;
        }
    }
    false
}
/// # Safety
/// All reachable HSP handles must satisfy BlastIntervalTree's contract.
pub unsafe fn s_midpoint_tree_contains_hsp(
    tree: &BlastIntervalTree,
    root_index: i32,
    in_hsp: &BlastHSP,
    in_q_start: i32,
    min_diag_separation: i32,
) -> bool {
    // SAFETY: traversal only dereferences retained indexed HSP allocations.
    unsafe {
        let mut node = root_index as usize;
        let region_start = in_hsp.subject.offset;
        let region_end = in_hsp.subject.end;
        while tree.nodes[node].hsp.is_none() {
            assert!(
                region_start >= tree.nodes[node].leftend && region_end <= tree.nodes[node].rightend
            );
            let mut tmp = tree.nodes[node].midptr;
            while tmp != 0 {
                let item = tree.nodes[tmp as usize];
                if s_hsp_is_contained(
                    in_hsp,
                    in_q_start,
                    item.hsp.unwrap().as_ref(),
                    item.leftptr,
                    min_diag_separation,
                ) {
                    return true;
                }
                tmp = item.midptr;
            }
            tmp = 0;
            let middle = (tree.nodes[node].leftend as i64 + tree.nodes[node].rightend as i64) / 2;
            if (region_end as i64) < middle {
                tmp = tree.nodes[node].leftptr;
            } else if (region_start as i64) > middle {
                tmp = tree.nodes[node].rightptr;
            }
            if tmp == 0 {
                return false;
            }
            node = tmp as usize;
        }
        s_hsp_is_contained(
            in_hsp,
            in_q_start,
            tree.nodes[node].hsp.unwrap().as_ref(),
            tree.nodes[node].leftptr,
            min_diag_separation,
        )
    }
}
/// # Safety
/// All reachable HSP allocations must remain live and their coordinate/context
/// keys must remain consistent with the tree and QueryInfo.
pub unsafe fn blast_interval_tree_contains_hsp(
    tree: &BlastIntervalTree,
    hsp: &BlastHSP,
    query_info: &BlastQueryInfo,
    min_diag_separation: i32,
) -> bool {
    // SAFETY: caller retains indexed HSP allocations for this traversal.
    unsafe {
        let mut node = 0usize;
        let query_start = s_get_query_strand_offset(query_info, hsp.context);
        let region_start = query_start + hsp.query.offset;
        let region_end = query_start + hsp.query.end;
        assert!(
            region_start >= tree.nodes[node].leftend && region_end <= tree.nodes[node].rightend
        );
        assert!(hsp.subject.offset >= tree.s_min && hsp.subject.end <= tree.s_max);
        assert!(hsp.query.offset <= hsp.query.end && hsp.subject.offset <= hsp.subject.end);
        while tree.nodes[node].hsp.is_none() {
            assert!(
                region_start >= tree.nodes[node].leftend && region_end <= tree.nodes[node].rightend
            );
            let mut tmp = tree.nodes[node].midptr;
            if tmp > 0
                && s_midpoint_tree_contains_hsp(tree, tmp, hsp, query_start, min_diag_separation)
            {
                return true;
            }
            tmp = 0;
            let middle = (tree.nodes[node].leftend as i64 + tree.nodes[node].rightend as i64) / 2;
            if (region_end as i64) < middle {
                tmp = tree.nodes[node].leftptr;
            } else if (region_start as i64) > middle {
                tmp = tree.nodes[node].rightptr;
            }
            if tmp == 0 {
                return false;
            }
            node = tmp as usize;
        }
        s_hsp_is_contained(
            hsp,
            query_start,
            tree.nodes[node].hsp.unwrap().as_ref(),
            tree.nodes[node].leftptr,
            min_diag_separation,
        )
    }
}
pub fn s_hsp_query_range_is_contained(
    in_hsp: &BlastHSP,
    in_q_start: i32,
    tree_hsp: &BlastHSP,
    tree_q_start: i32,
) -> i32 {
    if in_q_start != tree_q_start || in_hsp.score >= tree_hsp.score {
        return 0;
    }
    if tree_hsp.query.offset <= in_hsp.query.offset && tree_hsp.query.end >= in_hsp.query.end {
        return 1;
    }
    0
}
/// # Safety
/// All indexed HSPs must remain live; this operation requires QueryOnly indexing.
pub unsafe fn blast_interval_tree_num_redundant(
    tree: &BlastIntervalTree,
    hsp: &BlastHSP,
    query_info: &BlastQueryInfo,
) -> i32 {
    // SAFETY: caller retains indexed HSPs and uses the original query-only mode.
    unsafe {
        let mut node = 0usize;
        let in_query_start = s_get_query_strand_offset(query_info, hsp.context);
        let region_start = in_query_start + hsp.query.offset;
        let region_end = in_query_start + hsp.query.end;
        let mut num_redundant = 0;
        assert!(
            region_start >= tree.nodes[node].leftend && region_end <= tree.nodes[node].rightend
        );
        assert!(hsp.query.offset <= hsp.query.end && hsp.subject.offset <= hsp.subject.end);
        while tree.nodes[node].hsp.is_none() {
            assert!(
                region_start >= tree.nodes[node].leftend && region_end <= tree.nodes[node].rightend
            );
            let mut tmp = tree.nodes[node].midptr;
            while tmp != 0 {
                let item = tree.nodes[tmp as usize];
                num_redundant += s_hsp_query_range_is_contained(
                    hsp,
                    in_query_start,
                    item.hsp.unwrap().as_ref(),
                    item.leftptr,
                );
                tmp = item.midptr;
            }
            tmp = 0;
            let middle = (tree.nodes[node].leftend as i64 + tree.nodes[node].rightend as i64) / 2;
            if (region_end as i64) < middle {
                tmp = tree.nodes[node].leftptr;
            } else if (region_start as i64) > middle {
                tmp = tree.nodes[node].rightptr;
            }
            if tmp == 0 {
                return num_redundant;
            }
            node = tmp as usize;
        }
        num_redundant
            + s_hsp_query_range_is_contained(
                hsp,
                in_query_start,
                tree.nodes[node].hsp.unwrap().as_ref(),
                tree.nodes[node].leftptr,
            )
    }
}
