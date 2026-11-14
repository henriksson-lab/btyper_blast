//! Small validation driver. Search/results stay in memory; tabular display is stdout.
use btyper_blast::{
    core::{
        blast_encoding::NCBI4NA_TO_BLASTNA,
        blast_program::EBlastProgramType,
        blast_target_translation::{
            blast_hsp_get_target_translation, blast_target_translation_new,
        },
        genetic_code::STANDARD_GENETIC_CODE,
    },
    nucleotide_search::search_nucleotide_database,
    report::AlignmentReport,
    search_protein_database, NucleotideDatabase, NucleotideTask, PreparedNucleotideBatch,
    PreparedNucleotideQuery, PreparedProteinBatch, PreparedProteinQuery,
};
use std::{
    error::Error,
    fs::File,
    io::{self, BufRead, BufReader, Write},
};

fn read_queries(path: &str) -> Result<Vec<(String, Vec<u8>)>, Box<dyn Error>> {
    let mut records = Vec::new();
    let mut current: Option<(String, Vec<u8>)> = None;
    for line in BufReader::new(File::open(path)?).lines() {
        let line = line?;
        if let Some(header) = line.strip_prefix('>') {
            if let Some(record) = current.take() {
                records.push(record);
            }
            let id = header
                .split_whitespace()
                .next()
                .ok_or("empty FASTA identifier")?;
            current = Some((id.to_owned(), Vec::new()));
        } else if !line.trim().is_empty() {
            let record = current
                .as_mut()
                .ok_or("query residues before FASTA header")?;
            record
                .1
                .extend(line.bytes().filter(|b| !b.is_ascii_whitespace()));
        }
    }
    if let Some(record) = current {
        records.push(record);
    }
    if records.is_empty() {
        return Err("query FASTA is empty".into());
    }
    Ok(records)
}
fn write_hit(
    output: &mut impl Write,
    query_id: &str,
    subject_id: &str,
    h: AlignmentReport,
) -> io::Result<()> {
    writeln!(
        output,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.17e}\t{:.17e}\t{}\t{}\t{}",
        query_id,
        subject_id,
        h.qstart,
        h.qend,
        h.sstart,
        h.send,
        h.qframe,
        h.sframe,
        h.length,
        h.nident,
        h.mismatch,
        h.gaps,
        h.gapopen,
        h.score,
        h.bitscore,
        h.evalue,
        h.btop,
        h.qseq,
        h.sseq
    )
}
fn main() -> Result<(), Box<dyn Error>> {
    let mut args: Vec<String> = std::env::args().collect();
    let raw = args.get(1).is_some_and(|arg| arg == "--raw-json");
    if raw {
        args.remove(1);
    }
    if args.len() != 4 {
        return Err(
            "usage: native_search [--raw-json] <megablast|blastn|tblastn> <database.fna> <queries.fasta>".into(),
        );
    }
    if args[1] == "tblastn" {
        return protein_search(&args[2], &args[3], raw);
    }
    let task = match args[1].as_str() {
        "megablast" => NucleotideTask::Megablast,
        "blastn" => NucleotideTask::Blastn,
        _ => return Err("expected megablast or blastn".into()),
    };
    let database = NucleotideDatabase::from_fasta(BufReader::new(File::open(&args[2])?))?;
    let queries = read_queries(&args[3])?
        .into_iter()
        .map(|(id, seq)| PreparedNucleotideQuery::new(id, &seq, task))
        .collect::<Result<Vec<_>, _>>()?;
    let batch = PreparedNucleotideBatch::new(&queries)?;
    let results = search_nucleotide_database(&batch, &database)?;
    let mut output = io::BufWriter::new(io::stdout().lock());
    if raw {
        write_raw_results(&mut output, &results.results)?;
        output.flush()?;
        return Ok(());
    }
    let mut subject_range = Vec::new();
    for (query_index, hitlist) in results.results.hitlist_array.iter().enumerate() {
        let Some(hitlist) = hitlist else { continue };
        for list in hitlist.hsplist_array[..hitlist.hsplist_count as usize]
            .iter()
            .flatten()
        {
            let subject = database
                .sequence(list.oid as usize)
                .ok_or("invalid result OID")?;
            for hsp in list.hsp_array[..list.hspcnt as usize].iter().flatten() {
                let ctx = &results.query_info.contexts[hsp.context as usize];
                let start = ctx.query_offset as usize;
                let query = &batch.blastna_contexts()[start..start + ctx.query_length as usize];
                let origin = usize::try_from(hsp.subject.offset)?;
                let end = usize::try_from(hsp.subject.end)?;
                let residues = subject
                    .ncbi4na()
                    .get(origin..end)
                    .ok_or("invalid HSP subject interval")?;
                subject_range.clear();
                subject_range.extend(
                    residues
                        .iter()
                        .map(|base| NCBI4NA_TO_BLASTNA[*base as usize]),
                );
                let hit = AlignmentReport::nucleotide_range(
                    hsp,
                    query,
                    &subject_range,
                    hsp.subject.offset,
                )?;
                write_hit(
                    &mut output,
                    &results.query_ids[query_index],
                    subject.id(),
                    hit,
                )?;
            }
        }
    }
    output.flush()?;
    Ok(())
}

fn protein_search(database_path: &str, query_path: &str, raw: bool) -> Result<(), Box<dyn Error>> {
    let database = NucleotideDatabase::from_fasta(BufReader::new(File::open(database_path)?))?;
    let queries = read_queries(query_path)?
        .into_iter()
        .map(|(id, seq)| PreparedProteinQuery::new(id, &seq))
        .collect::<Result<Vec<_>, _>>()?;
    let batch = PreparedProteinBatch::new(&queries)?;
    let results = search_protein_database(&batch, &database)?;
    let mut output = io::BufWriter::new(io::stdout().lock());
    if raw {
        write_raw_results(&mut output, &results.results)?;
        output.flush()?;
        return Ok(());
    }
    for (query_index, hitlist) in results.results.hitlist_array.iter().enumerate() {
        let Some(hitlist) = hitlist else { continue };
        for list in hitlist.hsplist_array[..hitlist.hsplist_count as usize]
            .iter()
            .flatten()
        {
            let subject = database
                .sequence(list.oid as usize)
                .ok_or("invalid result OID")?;
            let mut target = blast_target_translation_new(
                subject.ncbi4na(),
                &STANDARD_GENETIC_CODE,
                EBlastProgramType::Tblastn,
                false,
            );
            for hsp in list.hsp_array[..list.hspcnt as usize].iter().flatten() {
                let view = blast_hsp_get_target_translation(&mut target, Some(hsp))
                    .ok_or("missing subject translation for display")?;
                let ctx = &results.query_info.contexts[hsp.context as usize];
                let start = ctx.query_offset as usize;
                let query = &batch.original_contexts()[start..start + ctx.query_length as usize];
                let hit = AlignmentReport::translated_subject_range(
                    hsp,
                    query,
                    view.residues,
                    view.start_shift,
                    subject.len() as i32,
                )?;
                write_hit(
                    &mut output,
                    &results.query_ids[query_index],
                    subject.id(),
                    hit,
                )?;
            }
        }
    }
    output.flush()?;
    Ok(())
}

fn write_raw_results(
    output: &mut impl Write,
    results: &btyper_blast::core::blast_hits::BlastHSPResults,
) -> io::Result<()> {
    writeln!(
        output,
        "{{\"kind\":\"results\",\"num_queries\":{}}}",
        results.num_queries
    )?;
    for (query_index, hitlist) in results.hitlist_array.iter().enumerate() {
        let count = hitlist.as_ref().map_or(0, |list| list.hsplist_count);
        writeln!(
            output,
            "{{\"kind\":\"query\",\"query_index\":{},\"hsplist_count\":{}}}",
            query_index, count
        )?;
        let Some(hitlist) = hitlist else { continue };
        for (subject_list_index, list) in hitlist.hsplist_array[..count as usize].iter().enumerate()
        {
            let Some(list) = list else { continue };
            writeln!(output, "{{\"kind\":\"subject\",\"query_index\":{},\"subject_list_index\":{},\"oid\":{},\"list_query_index\":{},\"hspcnt\":{},\"best_evalue\":{:.17e},\"best_evalue_bits\":\"{:016x}\"}}", query_index, subject_list_index, list.oid,list.query_index,list.hspcnt,list.best_evalue,list.best_evalue.to_bits())?;
            for (hsp_index, hsp) in list.hsp_array[..list.hspcnt as usize].iter().enumerate() {
                let Some(hsp) = hsp else { continue };
                write!(output, "{{\"kind\":\"hsp\",\"query_index\":{},\"subject_list_index\":{},\"oid\":{},\"hsp_index\":{},\"score\":{},\"num_ident\":{},\"num_positives\":{},\"bit_score\":{:.17e},\"bit_score_bits\":\"{:016x}\",\"evalue\":{:.17e},\"evalue_bits\":\"{:016x}\",\"context\":{},\"num\":{},\"comp_adjustment_method\":{},\"query\":", query_index,subject_list_index,list.oid,hsp_index,hsp.score,hsp.num_ident,hsp.num_positives,hsp.bit_score,hsp.bit_score.to_bits(),hsp.evalue,hsp.evalue.to_bits(),hsp.context,hsp.num,hsp.comp_adjustment_method)?;
                write_raw_segment(output, &hsp.query)?;
                write!(output, ",\"subject\":")?;
                write_raw_segment(output, &hsp.subject)?;
                write!(output, ",\"gap_info\":")?;
                if let Some(script) = &hsp.gap_info {
                    write!(output, "{{\"size\":{},\"op_type\":[", script.size)?;
                    for i in 0..script.size as usize {
                        if i != 0 {
                            write!(output, ",")?;
                        }
                        write!(output, "{}", script.op_type[i] as i32)?;
                    }
                    write!(output, "],\"num\":[")?;
                    for i in 0..script.size as usize {
                        if i != 0 {
                            write!(output, ",")?;
                        }
                        write!(output, "{}", script.num[i])?;
                    }
                    write!(output, "]}}")?;
                } else {
                    write!(output, "null")?;
                }
                writeln!(output, "}}")?;
            }
        }
    }
    Ok(())
}
fn write_raw_segment(
    output: &mut impl Write,
    segment: &btyper_blast::core::blast_hits::BlastSeg,
) -> io::Result<()> {
    write!(
        output,
        "{{\"frame\":{},\"offset\":{},\"end\":{},\"gapped_start\":{}}}",
        segment.frame, segment.offset, segment.end, segment.gapped_start
    )
}
