// NCBI BLAST+ 2.17.0+, compo_heap.c.
// Authors: NCBI BLAST developers. See LICENSE-NCBI.
pub const HEAP_INITIAL_CAPACITY: i32 = 100;
pub const HEAP_RESIZE_FACTOR: f64 = 1.5;
pub const HEAP_MIN_RESIZE: i32 = 100;
#[derive(Debug)]
pub struct BlastCompoHeapRecord<T> {
    pub best_evalue: f64,
    pub best_score: i32,
    pub subject_index: i32,
    pub these_alignments: Option<T>,
}
impl<T> Default for BlastCompoHeapRecord<T> {
    fn default() -> Self {
        Self {
            best_evalue: 0.0,
            best_score: 0,
            subject_index: 0,
            these_alignments: None,
        }
    }
}
#[derive(Debug)]
pub struct BlastCompoHeap<T> {
    pub n: i32,
    pub capacity: i32,
    pub heap_threshold: i32,
    pub ecutoff: f64,
    pub worst_evalue: f64,
    pub array: Option<Vec<BlastCompoHeapRecord<T>>>,
    pub heap_array: Option<Vec<BlastCompoHeapRecord<T>>>,
}
impl<T> Default for BlastCompoHeap<T> {
    fn default() -> Self {
        Self {
            n: 0,
            capacity: 0,
            heap_threshold: 0,
            ecutoff: 0.0,
            worst_evalue: 0.0,
            array: None,
            heap_array: None,
        }
    }
}
pub fn s_compo_heap_record_compare<T>(
    place1: &BlastCompoHeapRecord<T>,
    place2: &BlastCompoHeapRecord<T>,
) -> bool {
    let mut result = if place1.best_evalue > place2.best_evalue {
        1
    } else if place1.best_evalue < place2.best_evalue {
        -1
    } else {
        0
    };
    if result == 0 {
        result = if place2.best_score > place1.best_score {
            1
        } else if place2.best_score < place1.best_score {
            -1
        } else {
            0
        };
        if result == 0 {
            result = if place2.subject_index > place1.subject_index {
                1
            } else if place2.subject_index < place1.subject_index {
                -1
            } else {
                0
            };
        }
    }
    result > 0
}
pub fn s_compo_heap_record_swap<T>(
    record1: &mut BlastCompoHeapRecord<T>,
    record2: &mut BlastCompoHeapRecord<T>,
) {
    std::mem::swap(record1, record2);
}
pub fn s_compo_heap_is_valid<T>(heap_array: &[BlastCompoHeapRecord<T>], i: i32, n: i32) -> bool {
    let left = 2 * i;
    let right = 2 * i + 1;
    if right <= n {
        return !s_compo_heap_record_compare(&heap_array[right as usize], &heap_array[i as usize])
            && s_compo_heap_is_valid(heap_array, right, n);
    }
    if left <= n {
        return !s_compo_heap_record_compare(&heap_array[left as usize], &heap_array[i as usize])
            && s_compo_heap_is_valid(heap_array, left, n);
    }
    true
}
fn swap_slots<T>(array: &mut [BlastCompoHeapRecord<T>], i: usize, j: usize) {
    if i < j {
        let (a, b) = array.split_at_mut(j);
        s_compo_heap_record_swap(&mut a[i], &mut b[0]);
    } else if i > j {
        let (a, b) = array.split_at_mut(i);
        s_compo_heap_record_swap(&mut b[0], &mut a[j]);
    }
}
pub fn s_compo_heapify_down<T>(heap_array: &mut [BlastCompoHeapRecord<T>], top: i32, n: i32) {
    let mut largest = top;
    loop {
        let i = largest;
        let left = 2 * i;
        let right = 2 * i + 1;
        if left <= n
            && s_compo_heap_record_compare(&heap_array[left as usize], &heap_array[i as usize])
        {
            largest = left;
        } else {
            largest = i;
        }
        if right <= n
            && s_compo_heap_record_compare(
                &heap_array[right as usize],
                &heap_array[largest as usize],
            )
        {
            largest = right;
        }
        if largest != i {
            swap_slots(heap_array, i as usize, largest as usize);
        } else {
            break;
        }
    }
}
pub fn s_compo_heapify_up<T>(heap_array: &mut [BlastCompoHeapRecord<T>], mut i: i32) {
    let mut parent = i / 2;
    while parent >= 1
        && s_compo_heap_record_compare(&heap_array[i as usize], &heap_array[parent as usize])
    {
        swap_slots(heap_array, i as usize, parent as usize);
        i = parent;
        parent /= 2;
    }
}
pub fn s_convert_to_heap<T>(heap: &mut BlastCompoHeap<T>) {
    if let Some(array) = heap.array.take() {
        heap.heap_array = Some(array);
        let mut i = heap.n / 2;
        while i >= 1 {
            s_compo_heapify_down(heap.heap_array.as_mut().unwrap(), i, heap.n);
            i -= 1;
        }
    }
}
pub fn blast_compo_heap_would_insert<T>(
    heap: &mut BlastCompoHeap<T>,
    evalue: f64,
    score: i32,
    subject_index: i32,
) -> bool {
    if heap.n < heap.heap_threshold || evalue <= heap.ecutoff || evalue < heap.worst_evalue {
        true
    } else {
        if heap.heap_array.is_none() {
            s_convert_to_heap(heap);
        }
        let candidate = BlastCompoHeapRecord {
            best_evalue: evalue,
            best_score: score,
            subject_index,
            these_alignments: None,
        };
        s_compo_heap_record_compare(&heap.heap_array.as_ref().unwrap()[1], &candidate)
    }
}
pub fn s_comp_heap_record_insert_at_end<T>(
    array: &mut Vec<BlastCompoHeapRecord<T>>,
    length: &mut i32,
    capacity: &mut i32,
    alignments: T,
    evalue: f64,
    score: i32,
    subject_index: i32,
) -> i32 {
    if *length >= *capacity {
        let new_capacity =
            (HEAP_MIN_RESIZE + *capacity).max((HEAP_RESIZE_FACTOR * *capacity as f64) as i32);
        array.resize_with(new_capacity as usize + 1, BlastCompoHeapRecord::default);
        *capacity = new_capacity;
    }
    *length += 1;
    array[*length as usize] = BlastCompoHeapRecord {
        best_evalue: evalue,
        best_score: score,
        subject_index,
        these_alignments: Some(alignments),
    };
    0
}
pub fn blast_compo_heap_insert<T>(
    heap: &mut BlastCompoHeap<T>,
    alignments: T,
    evalue: f64,
    score: i32,
    subject_index: i32,
    discarded_alignments: &mut Option<T>,
) -> i32 {
    *discarded_alignments = None;
    if heap.array.is_some() && heap.n >= heap.heap_threshold {
        s_convert_to_heap(heap);
    }
    if let Some(array) = heap.array.as_mut() {
        let status = s_comp_heap_record_insert_at_end(
            array,
            &mut heap.n,
            &mut heap.capacity,
            alignments,
            evalue,
            score,
            subject_index,
        );
        if status != 0 {
            return -1;
        }
        if heap.worst_evalue < evalue {
            heap.worst_evalue = evalue;
        }
    } else {
        let array = heap.heap_array.as_mut().unwrap();
        if heap.n < heap.heap_threshold
            || (evalue <= heap.ecutoff && heap.worst_evalue <= heap.ecutoff)
        {
            let status = s_comp_heap_record_insert_at_end(
                array,
                &mut heap.n,
                &mut heap.capacity,
                alignments,
                evalue,
                score,
                subject_index,
            );
            if status != 0 {
                return -1;
            }
            s_compo_heapify_up(array, heap.n);
        } else {
            let mut candidate = BlastCompoHeapRecord {
                best_evalue: evalue,
                best_score: score,
                subject_index,
                these_alignments: Some(alignments),
            };
            if s_compo_heap_record_compare(&array[1], &candidate) {
                *discarded_alignments = array[1].these_alignments.take();
                array[1] = candidate;
            } else {
                *discarded_alignments = candidate.these_alignments.take();
            }
            s_compo_heapify_down(array, 1, heap.n);
        }
        heap.worst_evalue = array[1].best_evalue;
    }
    0
}
pub fn blast_compo_heap_filled_to_cutoff<T>(heap: &BlastCompoHeap<T>) -> bool {
    heap.n >= heap.heap_threshold && heap.worst_evalue <= heap.ecutoff
}
pub fn blast_compo_heap_initialize<T>(
    heap: &mut BlastCompoHeap<T>,
    heap_threshold: i32,
    ecutoff: f64,
) -> i32 {
    assert!(heap_threshold >= 0);
    heap.n = 0;
    heap.heap_threshold = heap_threshold;
    heap.ecutoff = ecutoff;
    heap.heap_array = None;
    heap.capacity = HEAP_INITIAL_CAPACITY.min(heap_threshold);
    heap.worst_evalue = 0.0;
    heap.array = Some(
        (0..heap.capacity + 1)
            .map(|_| BlastCompoHeapRecord::default())
            .collect(),
    );
    0
}
pub fn blast_compo_heap_release<T>(heap: &mut BlastCompoHeap<T>) {
    heap.heap_array = None;
    heap.array = None;
    heap.n = 0;
    heap.capacity = 0;
    heap.heap_threshold = 0;
}
pub fn blast_compo_heap_pop<T>(heap: &mut BlastCompoHeap<T>) -> Option<T> {
    s_convert_to_heap(heap);
    if heap.n > 0 {
        let array = heap.heap_array.as_mut().unwrap();
        let result = array[1].these_alignments.take();
        heap.n -= 1;
        if heap.n > 0 {
            array[1] = std::mem::take(&mut array[heap.n as usize + 1]);
            s_compo_heapify_down(array, 1, heap.n);
        }
        result
    } else {
        None
    }
}
