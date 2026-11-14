use btyper_blast::core::blast_extend::*;
use btyper_blast::core::na_ungapped::BlastUngappedData;

fn verify(fixture: &str) {
    let mut e = BlastExtendWord::default();
    let mut list = blast_init_hit_list_new();
    let mut dst = blast_init_hit_list_new();
    let mut moved = false;
    let mut move_source_pending = false;
    let mut sorted = false;
    for line in fixture.lines() {
        let mut t = line.split_whitespace();
        let op = t.next().unwrap();
        let mut n = || t.next().unwrap().parse::<i32>().unwrap();
        match op {
            "new" => {
                let q = n();
                let hash = n() != 0;
                let window = n();
                e = blast_extend_word_new(q as u32, hash, window);
            }
            "dirty" => {
                let offset = n();
                if let Some(d) = &mut e.diag_table {
                    d.offset = offset;
                    for (i, c) in d.hit_level_array.iter_mut().enumerate() {
                        c.last_hit_and_flag =
                            ((i as i32 - 17) as u32 & 0x7fffffff) | ((i as u32 % 2) << 31);
                        if let Some(lengths) = &mut d.hit_len_array {
                            lengths[i] = (i + 129) as u8;
                        }
                    }
                } else {
                    let h = e.hash_table.as_mut().unwrap();
                    h.offset = offset;
                    h.occupancy = 3;
                    for (i, v) in h.backbone.iter_mut().enumerate() {
                        *v = i as u32 % 3;
                    }
                    for i in 0..3 {
                        h.chain[i] = DiagHashCell {
                            diag: -(i as i32),
                            level_and_saved: ((i as i32 - 13) as u32 & 0x7fffffff)
                                | ((i as u32 % 2) << 31),
                            hit_len: i as i32 + 23,
                            next: i as u32 % 3,
                        };
                    }
                }
            }
            "exit" => {
                let len = n();
                assert_eq!(blast_extend_word_exit(Some(&mut e), len) as i32, n());
            }
            "state" => {
                let d = e.diag_table.as_ref().unwrap();
                assert_eq!(
                    [
                        d.offset,
                        d.diag_array_length,
                        d.diag_mask,
                        d.multiple_hits as i32,
                        d.actual_window
                    ],
                    [n(), n(), n(), n(), n()]
                );
                for (i, c) in d.hit_level_array.iter().enumerate() {
                    assert_eq!(
                        [
                            (c.last_hit_and_flag.wrapping_shl(1) as i32) >> 1,
                            (c.last_hit_and_flag >> 31) as i32,
                            d.hit_len_array.as_ref().map_or(-1, |v| v[i] as i32)
                        ],
                        [n(), n(), n()]
                    );
                }
            }
            "hash" => {
                let h = e.hash_table.as_ref().unwrap();
                assert_eq!(
                    [
                        h.offset,
                        h.num_buckets as i32,
                        h.occupancy as i32,
                        h.capacity as i32
                    ],
                    [n(), n(), n(), n()]
                );
                for value in &h.backbone {
                    assert_eq!(*value as i32, n());
                }
                for c in &h.chain[..3] {
                    assert_eq!(
                        [
                            c.diag,
                            (c.level_and_saved.wrapping_shl(1) as i32) >> 1,
                            (c.level_and_saved >> 31) as i32,
                            c.hit_len,
                            c.next as i32
                        ],
                        [n(), n(), n(), n(), n()]
                    );
                }
            }
            "null" => {
                assert_eq!(blast_extend_word_exit(None, 2) as i32, n());
                assert_eq!(s_blast_diag_clear(None), n());
            }
            "list_new" => {
                list = blast_init_hit_list_new();
            }
            "lock" => list.do_not_reallocate = true,
            "unlock" => list.do_not_reallocate = false,
            "save" => {
                let q = n();
                let s = n();
                let present = n() != 0;
                let data = present.then_some(BlastUngappedData {
                    q_start: q % 17,
                    s_start: q % 13,
                    length: 1 + q % 31,
                    score: q % 19 - 9,
                });
                assert_eq!(blast_save_initial_hit(&mut list, q, s, data) as i32, n());
            }
            "sort" => {
                blast_init_hit_list_sort_by_score(&mut list);
                sorted = true;
            }
            "move" => {
                blast_init_hit_list_move(&mut dst, &mut list);
                moved = true;
                move_source_pending = true;
            }
            "reset" => {
                blast_init_hit_list_reset(&mut dst);
                sorted = false;
            }
            "hsp" => blast_save_init_hsp(&mut dst, 3, 5, 7, 9, 11, 13),
            "list" => {
                let l = if move_source_pending {
                    move_source_pending = false;
                    &list
                } else if moved {
                    &dst
                } else {
                    &list
                };
                assert_eq!(
                    [
                        l.total,
                        l.allocated,
                        l.do_not_reallocate as i32,
                        blast_init_hit_list_is_sorted_by_score(l) as i32
                    ],
                    [n(), n(), n(), n()],
                    "{line}"
                );
                let mut expected = Vec::new();
                let mut actual = Vec::new();
                for h in &l.init_hsp_array {
                    let q = n();
                    let s = n();
                    let present = n();
                    let d = if present != 0 {
                        [n(), n(), n(), n()]
                    } else {
                        [0; 4]
                    };
                    expected.push((q, s, present, d));
                    actual.push((
                        h.offsets.q_off,
                        h.offsets.s_off,
                        h.ungapped_data.is_some() as i32,
                        h.ungapped_data
                            .as_ref()
                            .map_or([0; 4], |d| [d.q_start, d.s_start, d.length, d.score]),
                    ));
                }
                if sorted {
                    // C qsort leaves equal-key order unspecified across platforms.
                    // Require the same key at every sorted position and exactly the same hits.
                    for (a, b) in actual.iter().zip(&expected) {
                        assert_eq!((a.2, a.3), (b.2, b.3));
                    }
                    actual.sort();
                    expected.sort();
                }
                assert_eq!(actual, expected);
            }
            _ => panic!("unknown record {line}"),
        }
        assert_eq!(t.next(), None, "extra tokens: {line}");
    }
}
#[test]
fn original_extension_state_and_hit_storage() {
    verify(include_str!("fixtures/extend-2.17.0.txt"));
}
#[test]
#[ignore = "requires original C output via BLAST_EXTEND_FIXTURE"]
fn freshly_compiled_extension_state_and_hit_storage() {
    verify(&std::fs::read_to_string(std::env::var("BLAST_EXTEND_FIXTURE").unwrap()).unwrap());
}
