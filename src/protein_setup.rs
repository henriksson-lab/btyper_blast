//! Native standard tblastn setup and six-frame preliminary gapped search.
//! Composition-adjusted redo, final traceback, significance and hit collection
//! follow this stage; these preliminary HSPs are not final BLAST results.
//! Assembly follows pinned BlastSetup_ScoreBlkInit/BLAST_GapAlignSetUp and
//! BLAST_PreliminarySearchEngine setup; search kernels retain their separate
//! one-to-one upstream translations. Native owners replace SeqSrc callbacks.
use crate::core::{
    aa_ungapped::{blast_aa_word_finder, ProteinSequence},
    blast_aalookup::*,
    blast_aascan::{blast_choose_protein_scan_subject, AminoAcidSubject},
    blast_diagnostics::{BlastGappedStats, BlastUngappedStats},
    blast_dynamic_align::{BlastDynamicAlign, DynamicSequence},
    blast_effective_lengths::*,
    blast_encoding::{EBlastEncoding, BLASTAA_SEQ_CODE},
    blast_extend::*,
    blast_gapalign::{blast_gap_align_struct_new, BlastGapAlignStruct},
    blast_hits::{blast_hsp_list_new, BlastHSPList},
    blast_lookup::SSeqRange,
    blast_nascan::BlastOffsetPair,
    blast_options::*,
    blast_parameters::*,
    blast_prelim::{blast_get_gapped_score, GapSequence, GapSubject},
    blast_program::EBlastProgramType,
    blast_protein_matrix::*,
    blast_query_info::BlastQueryInfo,
    blast_stat::*,
    blast_stat_spouge::{blast_gumbel_blk_calc, BlastGumbelBlk},
    blast_util::blast_get_all_translations,
    genetic_code::STANDARD_GENETIC_CODE,
};
use crate::{NucleotideDatabase, NucleotideSequence, PreparedProteinBatch};
use std::fmt;
const PROGRAM: EBlastProgramType = EBlastProgramType::Tblastn;
const MATRIX: &[u8] = b"BLOSUM62";
// blast_seqsrc.h; the standard SeqDB fallback, including v5 databases.
const BLAST_SEQSRC_MINLENGTH: i32 = 10;
const FRAME_ORDER: [i16; 6] = [1, 2, 3, -1, -2, -3];
#[derive(Debug)]
pub enum ProteinSetupError {
    EmptyDatabase,
    InvalidSubjectLength { length: usize },
    DatabaseTooLarge,
    Upstream { stage: &'static str, status: i32 },
    MissingState { stage: &'static str },
}
impl fmt::Display for ProteinSetupError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyDatabase => write!(f, "tblastn requires a nonempty nucleotide database"),
            Self::InvalidSubjectLength { length } => {
                write!(f, "invalid tblastn subject length {length}")
            }
            Self::DatabaseTooLarge => {
                write!(f, "database dimensions exceed native BLAST integer limits")
            }
            Self::Upstream { stage, status } => {
                write!(f, "BLAST {stage} failed with status {status}")
            }
            Self::MissingState { stage } => {
                write!(f, "BLAST {stage} did not initialize required state")
            }
        }
    }
}
impl std::error::Error for ProteinSetupError {}
fn checked(status: i32, stage: &'static str) -> Result<(), ProteinSetupError> {
    if status == 0 {
        Ok(())
    } else {
        Err(ProteinSetupError::Upstream { stage, status })
    }
}
fn required<T>(value: Option<T>, stage: &'static str) -> Result<T, ProteinSetupError> {
    value.ok_or(ProteinSetupError::MissingState { stage })
}
#[derive(Debug)]
pub struct ProteinPreliminaryFrame {
    pub frame: i16,
    pub hsp_list: Box<BlastHSPList>,
}
#[derive(Debug)]
pub struct ProteinPreliminarySubject {
    pub frames: Vec<ProteinPreliminaryFrame>,
}
impl ProteinPreliminarySubject {
    pub fn hsp_count(&self) -> usize {
        self.frames.iter().map(|f| f.hsp_list.hspcnt as usize).sum()
    }
    fn empty() -> Self {
        Self {
            frames: FRAME_ORDER
                .into_iter()
                .map(|frame| ProteinPreliminaryFrame {
                    frame,
                    hsp_list: blast_hsp_list_new(i32::MAX),
                })
                .collect(),
        }
    }
}
/// Reusable native BLOSUM62/11/1, word3/threshold13, two-hit/window40 setup.
/// Masked contexts feed lookup/ungapped search. The batch retains original
/// contexts independently for subsequent composition-adjusted traceback.
pub struct ProteinInitialSearch<'a> {
    batch: &'a PreparedProteinBatch,
    query_info: BlastQueryInfo,
    pub score_block: BlastScoreBlk,
    pub scoring_options: BlastScoringOptions,
    pub extension_parameters: BlastExtensionParameters,
    pub word_parameters: BlastInitialWordParameters,
    pub hit_parameters: BlastHitSavingParameters,
    pub warnings: Vec<StatisticalWarning>,
    pub effective_length_warnings: Vec<BlastEffectiveLengthWarning>,
    lookup: BlastAaLookupTable,
    state: BlastExtendWord,
    pairs: Vec<BlastOffsetPair>,
    initial_hits: BlastInitHitList,
    ungapped_stats: BlastUngappedStats,
    gapped_stats: BlastGappedStats,
    gap_align: Option<BlastGapAlignStruct<'static>>,
}
impl<'a> ProteinInitialSearch<'a> {
    pub fn new(
        batch: &'a PreparedProteinBatch,
        database: &NucleotideDatabase,
    ) -> Result<Self, ProteinSetupError> {
        Self::new_with_effective_searchspaces(batch, database, &[])
    }
    /// Split-query setup supplies the full query's effective spaces before
    /// deriving chunk hit and word thresholds, as in the original splitter.
    pub(crate) fn new_with_effective_searchspaces(
        batch: &'a PreparedProteinBatch,
        database: &NucleotideDatabase,
        effective_searchspaces: &[i64],
    ) -> Result<Self, ProteinSetupError> {
        Self::new_internal(batch, database, effective_searchspaces, true)
    }
    /// Original full-query traceback setup does not construct a full lookup
    /// when the preliminary search uses split query factories.
    pub(crate) fn new_for_traceback(
        batch: &'a PreparedProteinBatch,
        database: &NucleotideDatabase,
    ) -> Result<Self, ProteinSetupError> {
        Self::new_internal(batch, database, &[], false)
    }
    fn new_internal(
        batch: &'a PreparedProteinBatch,
        database: &NucleotideDatabase,
        effective_searchspaces: &[i64],
        build_lookup: bool,
    ) -> Result<Self, ProteinSetupError> {
        if database.is_empty() {
            return Err(ProteinSetupError::EmptyDatabase);
        }
        if database.total_length() > i64::MAX as u64
            || database.len() > i32::MAX as usize
            || database.max_length() > (i32::MAX as usize - 4) / 2
        {
            return Err(ProteinSetupError::DatabaseTooLarge);
        }
        // Standard SeqDB volumes do not store their minimum sequence length:
        // CSeqDBIdxFile.m_MinLen remains zero and CSeqDBImpl substitutes
        // BLAST_SEQSRC_MINLENGTH (10nt). BLAST_GapAlignSetUp then divides by
        // three for translated subjects. Using the actual shortest input here
        // raises the preliminary ungapped threshold and loses weak CBS hits.
        let min_protein_length = BLAST_SEQSRC_MINLENGTH / 3;
        let scoring_options = blast_scoring_options_new(PROGRAM);
        let mut score_block = BlastScoreBlk {
            alphabet_code: BLASTAA_SEQ_CODE as u8,
            alphabet_size: 28,
            protein_alphabet: true,
            ambiguous_res: vec![21],
            scale_factor: 1.0,
            matrix: vec![vec![0; 28]; 28],
            ..Default::default()
        };
        checked(
            blast_score_blk_protein_matrix_load(&mut score_block, MATRIX) as i32,
            "protein matrix",
        )?;
        checked(
            blast_score_blk_max_score_set(&mut score_block) as i32,
            "score limits",
        )?;
        let mut query_info = batch.query_info().clone();
        let mut warnings = vec![];
        checked(
            blast_score_blk_kbp_ungapped_calc(
                PROGRAM,
                &mut score_block,
                batch.encoded_contexts(),
                &mut query_info,
                Some(&mut warnings),
            ) as i32,
            "ungapped statistics",
        )?;
        let mut gapped = vec![None; query_info.contexts.len()];
        for c in query_info.first_context..=query_info.last_context {
            if !query_info.contexts[c as usize].is_valid {
                continue;
            }
            let mut kbp = blast_karlin_blk_new();
            checked(
                blast_karlin_blk_gapped_load_from_tables(
                    Some(&mut kbp),
                    scoring_options.gap_open,
                    scoring_options.gap_extend,
                    Some(MATRIX),
                    true,
                ) as i32,
                "gapped statistics",
            )?;
            gapped[c as usize] = Some(kbp);
        }
        score_block.kbp_gap = Some(gapped);
        let mut gbp = BlastGumbelBlk::default();
        checked(
            blast_gumbel_blk_calc(
                Some(&mut gbp),
                scoring_options.gap_open,
                scoring_options.gap_extend,
                Some(MATRIX),
                None,
            ) as i32,
            "Gumbel statistics",
        )?;
        gbp.db_length = database.total_length() as i64 / 3;
        score_block.gbp = Some(gbp);
        let mut length_options = blast_effective_lengths_options_new();
        if !effective_searchspaces.is_empty() {
            // Original chunk setup indexes this full vector by LOCAL context.
            length_options.num_searchspaces = effective_searchspaces.len() as i32;
            length_options.searchsp_eff = Some(effective_searchspaces.to_vec());
        }
        let length_parameters = blast_effective_lengths_parameters_new(
            &length_options,
            database.total_length() as i64,
            database.len() as i32,
        );
        let mut effective_length_warnings = vec![];
        checked(
            blast_calc_eff_lengths(
                PROGRAM,
                &EffectiveLengthScoringOptions {
                    gapped_calculation: true,
                    reward: 0,
                    penalty: 0,
                    gap_open: scoring_options.gap_open,
                    gap_extend: scoring_options.gap_extend,
                },
                &length_parameters,
                Some(&score_block),
                Some(&mut query_info),
                Some(MATRIX),
                Some(&mut effective_length_warnings),
            ) as i32,
            "effective lengths",
        )?;
        let mut extension_options = blast_extension_options_new(PROGRAM, true);
        extension_options.composition_based_stats = 2;
        let mut extension_parameters = None;
        checked(
            blast_extension_parameters_new(
                PROGRAM,
                &extension_options,
                &score_block,
                &query_info,
                Some(&mut extension_parameters),
            ) as i32,
            "extension parameters",
        )?;
        let extension_parameters = required(extension_parameters, "extension parameters")?;
        let mut hit_parameters = None;
        checked(
            blast_hit_saving_parameters_new(
                PROGRAM,
                &blast_hit_saving_options_new(PROGRAM, true),
                &score_block,
                &query_info,
                min_protein_length,
                extension_options.composition_based_stats,
                Some(&mut hit_parameters),
            ) as i32,
            "hit-saving parameters",
        )?;
        let hit_parameters = required(hit_parameters, "hit-saving parameters")?;
        let mut word_parameters = None;
        checked(
            blast_initial_word_parameters_new(
                PROGRAM,
                &blast_initial_word_options_new(PROGRAM),
                &hit_parameters,
                &score_block,
                &query_info,
                (database.total_length() / database.len() as u64) as u32,
                Some(&mut word_parameters),
            ) as i32,
            "initial-word parameters",
        )?;
        let word_parameters = required(word_parameters, "initial-word parameters")?;
        let (lookup, state, pairs) = if build_lookup {
            let query_length = batch.encoded_contexts().len();
            let mut lookup =
                blast_aa_lookup_table_new(BLAST_WORDSIZE_PROT, BLAST_WORD_THRESHOLD_TBLASTN as f64);
            blast_aa_lookup_index_query(
                &mut lookup,
                &score_block.matrix,
                batch.encoded_contexts(),
                batch.lookup_regions(),
                0,
            );
            let bone = if query_length >= 65534 {
                EBoneType::Backbone
            } else {
                EBoneType::Smallbone
            };
            checked(
                blast_aa_lookup_finalize(&mut lookup, bone),
                "protein lookup finalization",
            )?;
            blast_choose_protein_scan_subject(&mut lookup);
            let state = blast_extend_word_new(
                query_length as u32,
                false,
                word_parameters.options.window_size,
            );
            let pairs = vec![BlastOffsetPair::default(); 4096 + lookup.longest_chain as usize];
            (lookup, state, pairs)
        } else {
            (
                BlastAaLookupTable::default(),
                BlastExtendWord::default(),
                Vec::new(),
            )
        };
        // BLAST_GapAlignSetUp creates this owner once per search, using the
        // database maximum length. Retain the source DP/traceback buffers.
        let gap_align = if build_lookup {
            let mut scoring = None;
            checked(
                blast_scoring_parameters_new(Some(&scoring_options), &score_block, &mut scoring)
                    as i32,
                "scoring parameters",
            )?;
            let scoring = required(scoring, "scoring parameters")?;
            let mut gap = None;
            checked(
                blast_gap_align_struct_new(
                    Some(&scoring),
                    Some(&extension_parameters),
                    database.max_length() as u32,
                    Some(&score_block),
                    Some(&mut gap),
                ) as i32,
                "gapped scratch",
            )?;
            Some(bind_gap_align(required(gap, "gapped scratch")?, &[]))
        } else {
            // Final-only setup delegates actual alignment scratch to the
            // independently owned composition-redo callbacks.
            None
        };
        Ok(Self {
            batch,
            query_info,
            score_block,
            scoring_options,
            extension_parameters,
            word_parameters,
            hit_parameters,
            warnings,
            effective_length_warnings,
            lookup,
            state,
            pairs,
            initial_hits: blast_init_hit_list_new(),
            ungapped_stats: BlastUngappedStats::default(),
            gapped_stats: BlastGappedStats::default(),
            gap_align,
        })
    }
    pub fn query_info(&self) -> &BlastQueryInfo {
        &self.query_info
    }
    pub fn batch(&self) -> &PreparedProteinBatch {
        self.batch
    }
    pub fn stats(&self) -> &BlastUngappedStats {
        &self.ungapped_stats
    }
    pub fn gapped_stats(&self) -> &BlastGappedStats {
        &self.gapped_stats
    }
    pub fn lookup(&self) -> &BlastAaLookupTable {
        &self.lookup
    }
    pub fn preliminary_subject(
        &mut self,
        subject: &NucleotideSequence,
    ) -> Result<ProteinPreliminarySubject, ProteinSetupError> {
        // The upstream translated search skips subjects shorter than word*3+2.
        if subject.len() < (self.lookup.word_length * 3 + 2) as usize {
            return Ok(ProteinPreliminarySubject::empty());
        }
        if subject.len() > (i32::MAX as usize - 4) / 2 {
            return Err(ProteinSetupError::InvalidSubjectLength {
                length: subject.len(),
            });
        }
        let mut translation = Vec::new();
        let mut offsets = [0u32; 7];
        checked(
            blast_get_all_translations(
                subject.packed_for_scanning(),
                EBlastEncoding::Ncbi2na,
                subject.len() as i32,
                &STANDARD_GENETIC_CODE,
                Some(&mut translation),
                Some(&mut offsets),
                None,
            ) as i32,
            "six-frame translation",
        )?;
        let mut scoring = None;
        checked(
            blast_scoring_parameters_new(
                Some(&self.scoring_options),
                &self.score_block,
                &mut scoring,
            ) as i32,
            "scoring parameters",
        )?;
        let scoring = required(scoring, "scoring parameters")?;
        let query = GapSequence {
            sequence: Some(DynamicSequence {
                data: self.batch.sequence_with_sentinels(),
                offset: 1,
            }),
            oof_sequence: None,
            length: self.batch.encoded_contexts().len() as i32,
        };
        let protein_query = ProteinSequence {
            sequence: self.batch.encoded_contexts(),
            length: query.length,
        };
        let mut frames = Vec::with_capacity(6);
        for (context, frame) in FRAME_ORDER.into_iter().enumerate() {
            blast_init_hit_list_reset(&mut self.initial_hits);
            let offset = offsets[context] as usize + 1;
            let length = (offsets[context + 1] - offsets[context] - 1) as i32;
            let mut list = None;
            if length >= self.lookup.word_length {
                let ranges = [SSeqRange {
                    left: 0,
                    right: length,
                }];
                let aa_subject = AminoAcidSubject {
                    sequence: &translation[offset..],
                    seq_ranges: &ranges,
                };
                let capacity = self.pairs.len() as i32;
                checked(
                    blast_aa_word_finder(
                        &aa_subject,
                        length,
                        &protein_query,
                        &self.query_info,
                        &self.lookup,
                        &self.score_block.matrix,
                        &self.word_parameters,
                        &mut self.state,
                        &mut self.pairs,
                        capacity,
                        &mut self.initial_hits,
                        Some(&mut self.ungapped_stats),
                    ) as i32,
                    "protein word finder",
                )?;
                let target = GapSubject {
                    sequence: DynamicSequence {
                        data: &translation,
                        offset: offset as i32,
                    },
                    length,
                    frame,
                };
                let mut gap = bind_gap_align(
                    required(self.gap_align.take(), "search-owned gapped scratch")?,
                    &self.score_block.matrix,
                );
                gap.dynamic.gap_x_dropoff = self.extension_parameters.gap_x_dropoff;
                let status = blast_get_gapped_score(
                    PROGRAM,
                    &query,
                    &self.query_info,
                    &target,
                    &mut gap,
                    &scoring,
                    &self.extension_parameters,
                    &self.hit_parameters,
                    &self.word_parameters,
                    &self.initial_hits,
                    Some(&mut list),
                    Some(&mut self.gapped_stats),
                    None,
                ) as i32;
                self.gap_align = Some(bind_gap_align(gap, &[]));
                checked(status, "protein preliminary gapped alignment")?;
            }
            frames.push(ProteinPreliminaryFrame {
                frame,
                hsp_list: list.unwrap_or_else(|| blast_hsp_list_new(i32::MAX)),
            });
        }
        Ok(ProteinPreliminarySubject { frames })
    }

    /// Search a chunk in an already translated physical frame allocation.
    /// The offset retains readable sentinels/adjacent residues for upstream DP.
    pub(crate) fn preliminary_translated_chunk(
        &mut self,
        translation: &[u8],
        physical_offset: i32,
        logical_length: i32,
        frame: i16,
    ) -> Result<Box<BlastHSPList>, ProteinSetupError> {
        if physical_offset <= 0
            || logical_length < 0
            || logical_length > 5_000_000
            || physical_offset as usize + logical_length as usize >= translation.len()
            || frame == 0
            || !(-3..=3).contains(&frame)
        {
            return Err(ProteinSetupError::InvalidSubjectLength {
                length: logical_length.max(0) as usize,
            });
        }
        blast_init_hit_list_reset(&mut self.initial_hits);
        if logical_length < self.lookup.word_length {
            return Ok(blast_hsp_list_new(i32::MAX));
        }
        let mut scoring = None;
        checked(
            blast_scoring_parameters_new(
                Some(&self.scoring_options),
                &self.score_block,
                &mut scoring,
            ) as i32,
            "scoring parameters",
        )?;
        let scoring = required(scoring, "scoring parameters")?;
        let query = GapSequence {
            sequence: Some(DynamicSequence {
                data: self.batch.sequence_with_sentinels(),
                offset: 1,
            }),
            oof_sequence: None,
            length: self.batch.encoded_contexts().len() as i32,
        };
        let protein_query = ProteinSequence {
            sequence: self.batch.encoded_contexts(),
            length: query.length,
        };
        let ranges = [SSeqRange {
            left: 0,
            right: logical_length,
        }];
        let aa_subject = AminoAcidSubject {
            sequence: &translation[physical_offset as usize..],
            seq_ranges: &ranges,
        };
        let capacity = self.pairs.len() as i32;
        checked(
            blast_aa_word_finder(
                &aa_subject,
                logical_length,
                &protein_query,
                &self.query_info,
                &self.lookup,
                &self.score_block.matrix,
                &self.word_parameters,
                &mut self.state,
                &mut self.pairs,
                capacity,
                &mut self.initial_hits,
                Some(&mut self.ungapped_stats),
            ) as i32,
            "protein chunk word finder",
        )?;
        let target = GapSubject {
            sequence: DynamicSequence {
                data: translation,
                offset: physical_offset,
            },
            length: logical_length,
            frame,
        };
        let mut list = None;
        let mut gap = bind_gap_align(
            required(self.gap_align.take(), "search-owned gapped scratch")?,
            &self.score_block.matrix,
        );
        gap.dynamic.gap_x_dropoff = self.extension_parameters.gap_x_dropoff;
        let status = blast_get_gapped_score(
            PROGRAM,
            &query,
            &self.query_info,
            &target,
            &mut gap,
            &scoring,
            &self.extension_parameters,
            &self.hit_parameters,
            &self.word_parameters,
            &self.initial_hits,
            Some(&mut list),
            Some(&mut self.gapped_stats),
            None,
        ) as i32;
        self.gap_align = Some(bind_gap_align(gap, &[]));
        checked(status, "protein chunk preliminary gapped alignment")?;
        Ok(list.unwrap_or_else(|| blast_hsp_list_new(i32::MAX)))
    }
}

// Short-lived matrix binding around search-owned stable scratch buffers.
fn bind_gap_align<'matrix>(
    gap: BlastGapAlignStruct<'_>,
    matrix: &'matrix [Vec<i32>],
) -> BlastGapAlignStruct<'matrix> {
    let dynamic = gap.dynamic;
    BlastGapAlignStruct {
        dynamic: BlastDynamicAlign {
            matrix,
            pssm: None,
            position_based: false,
            gap_x_dropoff: dynamic.gap_x_dropoff,
            dp_mem: dynamic.dp_mem,
            dp_mem_alloc: dynamic.dp_mem_alloc,
            state_struct: dynamic.state_struct,
            query_start: dynamic.query_start,
            query_stop: dynamic.query_stop,
            subject_start: dynamic.subject_start,
            subject_stop: dynamic.subject_stop,
            score: dynamic.score,
        },
        greedy_align_mem: gap.greedy_align_mem,
        edit_script: gap.edit_script,
        fwd_prelim_tback: gap.fwd_prelim_tback,
        rev_prelim_tback: gap.rev_prelim_tback,
        greedy_query_seed_start: gap.greedy_query_seed_start,
        greedy_subject_seed_start: gap.greedy_subject_seed_start,
    }
}
