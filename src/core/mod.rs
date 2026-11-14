//! Low-level kernels corresponding to the pinned BLAST C core.
//!
//! These retain upstream conventions and preconditions. They are exposed for
//! differential verification; the stable search API will be layered above them.

pub mod blast_diagnostics;
pub mod blast_effective_lengths;
pub mod blast_encoding;
pub mod blast_extend;
pub mod blast_filter;
pub mod blast_lookup;
pub mod blast_nalookup;
pub mod blast_nascan;
pub mod blast_options;
pub mod blast_parameters;
pub mod blast_program;
pub mod blast_protein_matrix;
pub mod blast_query_info;
pub mod blast_seg;
pub mod blast_setup;
pub mod blast_stat;
pub mod blast_stat_gapped;
pub mod blast_util;
pub mod gapinfo;
pub mod genetic_code;
pub mod greedy_align;
pub mod lookup_util;
pub mod masksubj;
pub mod na_ungapped;
pub mod ncbi_math;

pub mod blast_aalookup;
pub mod blast_aascan;
pub mod blast_dynamic_align;
pub mod blast_hits;

pub mod aa_ungapped;
pub mod blast_gapalign;
pub mod blast_prelim;
pub mod blast_protein_gapalign;
pub mod interval_tree;

pub mod blast_stat_spouge;

pub mod boost_erf;

pub mod compo_linear_algebra;

pub mod compo_optimize_target_freq;

pub mod blast_composition_adjustment;

pub mod blast_composition_redo;

pub mod compo_heap;

pub mod hsp_sumstat;

pub mod blast_hsp_link;

pub mod blast_hsp_evaluate;

pub mod blast_kappa;

pub mod blast_traceback;
pub mod hspfilter_collector;

pub mod blast_result_pipeline;

pub mod blast_prelim_evaluate;

pub mod blast_query_split;
pub mod blast_subject_chunk;
pub mod blast_target_translation;
