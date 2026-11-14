//! Native binding of the pinned tblastn complete-frame subject chunk loop.
//! Translation precedes amino-acid chunking, preserving codon phase, both
//! strands, sentinels and the database's original packed ambiguity placeholders.
use crate::core::{
    blast_encoding::EBlastEncoding,
    blast_hits::{
        blast_hsp_list_adjust_offsets, blast_hsp_list_new, blast_hsp_num_max, BlastHSPList,
    },
    blast_prelim_evaluate::s_blast_search_engine_one_context,
    blast_program::EBlastProgramType,
    blast_subject_chunk::*,
    blast_util::blast_get_all_translations,
    genetic_code::STANDARD_GENETIC_CODE,
};
use crate::protein_setup::{
    ProteinInitialSearch, ProteinPreliminaryFrame, ProteinPreliminarySubject, ProteinSetupError,
};
use crate::NucleotideSequence;
const PROGRAM: EBlastProgramType = EBlastProgramType::Tblastn;
const FRAME_ORDER: [i16; 6] = [1, 2, 3, -1, -2, -3];
fn checked(status: i32, stage: &'static str) -> Result<(), ProteinSetupError> {
    if status == 0 {
        Ok(())
    } else {
        Err(ProteinSetupError::Upstream { stage, status })
    }
}

/// Results already include the source OneContext per-chunk endpoint purges,
/// coordinate shifts and overlap merges. Feed them to engine core as
/// PreliminaryContextResult::Finished, without another endpoint purge.
pub(crate) fn preliminary_chunked_protein_subject(
    search: &mut ProteinInitialSearch<'_>,
    subject: &NucleotideSequence,
    oid: i32,
) -> Result<ProteinPreliminarySubject, ProteinSetupError> {
    if subject.len() > (i32::MAX as usize - 4) / 2 {
        return Err(ProteinSetupError::InvalidSubjectLength {
            length: subject.len(),
        });
    }
    if subject.len() < 11 {
        return Ok(ProteinPreliminarySubject {
            frames: FRAME_ORDER
                .into_iter()
                .map(|frame| ProteinPreliminaryFrame {
                    frame,
                    hsp_list: blast_hsp_list_new(i32::MAX),
                })
                .collect(),
        });
    }
    let mut translation = Vec::new();
    let mut offsets = [0u32; 7];
    // Source engine requests the original packed NCBI2na subject for this
    // preliminary translation. Ambiguity-aware translation belongs to redo.
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
        "chunked six-frame translation",
    )?;
    let hsp_num_max = blast_hsp_num_max(
        true,
        search
            .hit_parameters
            .options
            .as_ref()
            .ok_or(ProteinSetupError::MissingState {
                stage: "hit-saving options",
            })?
            .hsp_num_max,
    );
    let mut frames = Vec::with_capacity(6);
    for (context, frame) in FRAME_ORDER.into_iter().enumerate() {
        let frame_start = offsets[context] as i32 + 1;
        let frame_length = (offsets[context + 1] - offsets[context] - 1) as i32;
        let mut backup = SubjectSplitStruct::default();
        let mut chunk = SubjectChunk::default();
        s_backup_subject(frame_length, &mut backup, &mut chunk);
        let mut combined: Option<Box<BlastHSPList>> = None;
        while s_get_next_subject_chunk(&mut chunk, &mut backup, false, DBSEQ_CHUNK_OVERLAP)
            == SUBJECT_SPLIT_OK
        {
            let raw = search.preliminary_translated_chunk(
                &translation,
                frame_start + chunk.sequence_offset,
                chunk.length,
                frame,
            )?;
            let mut list = None;
            checked(
                s_blast_search_engine_one_context(PROGRAM, oid, || Ok(Some(raw)), &mut list) as i32,
                "translated chunk endpoint evaluation",
            )?;
            if let Some(list) = list.as_mut() {
                blast_hsp_list_adjust_offsets(list, backup.offset);
            }
            let overlap = if backup.offset == 0 {
                0
            } else {
                DBSEQ_CHUNK_OVERLAP
            };
            checked(
                blast_hsp_lists_merge(
                    &mut list,
                    &mut combined,
                    hsp_num_max,
                    &[backup.offset],
                    i32::MIN,
                    overlap,
                    true,
                    false,
                ) as i32,
                "translated subject chunk merging",
            )?;
        }
        frames.push(ProteinPreliminaryFrame {
            frame,
            hsp_list: combined.unwrap_or_else(|| blast_hsp_list_new(i32::MAX)),
        });
    }
    Ok(ProteinPreliminarySubject { frames })
}
