#!/usr/bin/env python3
"""Capture pinned original BLAST core results before Seq-align clamps E-values.

GDB owns only the child BLAST process. A source-header-compiled ABI reporter
supplies every struct offset/field width. Each completed exported
BLAST_ComputeTraceback_MT call supplies one final query batch after sorting and
pruning. All batches are retained, with explicit global/local query mapping.
"""
import json
import os
from pathlib import Path
import sys


def aggregate_batches(batches, nucleotide=False):
    """Retain raw per-batch provenance while presenting global query/context IDs.

    list_query_index deliberately stays the original core field. In nucleotide
    batches it is local, and batch metadata supplies its global interpretation.
    No HSP score, float, segment, seed or edit script is normalized here.
    """
    ordered = sorted(batches, key=lambda batch: batch["entry"]["batch_index"])
    rows = [{"kind": "results", "num_queries": sum(batch["rows"][0]["num_queries"] for batch in ordered)}]
    query_offset = context_offset = 0
    mapping = []
    for batch in ordered:
        entry = batch["entry"]
        count = batch["rows"][0]["num_queries"]
        if count != entry["num_queries"]:
            raise ValueError("returned results and entry QueryInfo query counts differ")
        mapping.append(dict(entry, query_index_offset=query_offset,
                            context_offset=context_offset,
                            hsp_count=sum(row["kind"] == "hsp" for row in batch["rows"])))
        for original in batch["rows"][1:]:
            row = dict(original)
            row["batch_index"] = entry["batch_index"]
            row["local_query_index"] = row["query_index"]
            row["query_index"] += query_offset
            if row["kind"] == "subject" and nucleotide:
                row["global_list_query_index"] = row["list_query_index"] + query_offset
            if row["kind"] == "hsp":
                row["local_context"] = row["context"]
                row["context"] += context_offset
            rows.append(row)
        query_offset += count
        context_offset += entry["last_context"] + 1
    return rows, mapping


def debugger_main(gdb):
    import struct

    config = json.loads(Path(os.environ["BTYPER_CORE_ORACLE_CONFIG"]).read_text())
    layout = config["layout"]
    if layout["pointer_size"] != 8 or layout["double_size"] != 8 or layout["byte_order"] != "little":
        raise RuntimeError("oracle requires the pinned Linux x86-64 little-endian ABI")
    state = {"snapshots": 0, "errors": [], "sort_returns": [],
             "capture_boundary": "BLAST_ComputeTraceback_MT return, after final sort/prune"}
    completed_batches = []
    next_batch_index = 0
    types = layout["types"]

    def memory(address, length):
        return bytes(gdb.selected_inferior().read_memory(address, length))

    def field(address, typename, name):
        descriptor = types[typename]["fields"][name]
        return memory(address + descriptor["offset"], descriptor["size"])

    def integer(address, typename, name, signed=True):
        return int.from_bytes(field(address, typename, name), "little", signed=signed)

    def pointer(address, typename, name):
        return integer(address, typename, name, False)

    def pointer_at(array, index):
        return int.from_bytes(memory(array + index * layout["pointer_size"], layout["pointer_size"]), "little")

    def floating(address, typename, name):
        data = field(address, typename, name)
        return struct.unpack("<d", data)[0], f"{int.from_bytes(data, 'little'):016x}"

    def segment(address):
        return {name: integer(address, "BlastSeg", name) for name in ("frame", "offset", "end", "gapped_start")}

    def capture(address):
        num_queries = integer(address, "BlastHSPResults", "num_queries")
        if not 0 <= num_queries <= 100000:
            raise ValueError(f"invalid core query count {num_queries}")
        rows = [{"kind": "results", "num_queries": num_queries}]
        query_array = pointer(address, "BlastHSPResults", "hitlist_array")
        for query_index in range(num_queries):
            query = pointer_at(query_array, query_index)
            if not query:
                continue
            count = integer(query, "BlastHitList", "hsplist_count")
            maximum = integer(query, "BlastHitList", "hsplist_max")
            if not 0 <= count <= 1000000:
                raise ValueError(f"invalid subject count {count}")
            worst, worst_bits = floating(query, "BlastHitList", "worst_evalue")
            rows.append({"kind": "query", "query_index": query_index,
                         "hsplist_count": count, "hsplist_max": maximum,
                         "worst_evalue": worst, "worst_evalue_bits": worst_bits,
                         "low_score": integer(query, "BlastHitList", "low_score"),
                         "heapified": bool(integer(query, "BlastHitList", "heapified", False)),
                         "storage": {name: integer(query, "BlastHitList", name)
                                     for name in ("hsplist_current", "num_hits")}})
            subject_array = pointer(query, "BlastHitList", "hsplist_array")
            for subject_list_index in range(count):
                subject = pointer_at(subject_array, subject_list_index)
                if not subject:
                    raise ValueError("NULL subject within live result count")
                oid = integer(subject, "BlastHSPList", "oid")
                list_query_index = integer(subject, "BlastHSPList", "query_index")
                hspcnt = integer(subject, "BlastHSPList", "hspcnt")
                allocated = integer(subject, "BlastHSPList", "allocated")
                if not 0 <= hspcnt <= allocated <= 10000000:
                    raise ValueError(f"invalid HSP count/capacity {hspcnt}/{allocated}")
                best, best_bits = floating(subject, "BlastHSPList", "best_evalue")
                rows.append({"kind": "subject", "query_index": query_index,
                             "subject_list_index": subject_list_index, "oid": oid,
                             "list_query_index": list_query_index, "hspcnt": hspcnt,
                             "best_evalue": best, "best_evalue_bits": best_bits,
                             "storage": {"allocated": allocated,
                                         "hsp_max": integer(subject, "BlastHSPList", "hsp_max"),
                                         "do_not_reallocate": bool(integer(subject, "BlastHSPList", "do_not_reallocate", False))}})
                hsp_array = pointer(subject, "BlastHSPList", "hsp_array")
                for hsp_index in range(hspcnt):
                    hsp = pointer_at(hsp_array, hsp_index)
                    if not hsp:
                        raise ValueError("NULL HSP within live result count")
                    row = {"kind": "hsp", "query_index": query_index,
                           "subject_list_index": subject_list_index, "oid": oid,
                           "list_query_index": list_query_index, "hsp_index": hsp_index,
                           "best_evalue_bits": best_bits}
                    for name in ("score", "num_ident", "num_positives", "context", "num", "comp_adjustment_method"):
                        row[name] = integer(hsp, "BlastHSP", name)
                    for name in ("bit_score", "evalue"):
                        row[name], row[name + "_bits"] = floating(hsp, "BlastHSP", name)
                    for name in ("query", "subject"):
                        row[name] = segment(hsp + types["BlastHSP"]["fields"][name]["offset"])
                    script = pointer(hsp, "BlastHSP", "gap_info")
                    row["gap_info"] = None
                    if script:
                        size = integer(script, "GapEditScript", "size")
                        if not 0 <= size <= 10000000:
                            raise ValueError(f"invalid edit script size {size}")
                        operations = pointer(script, "GapEditScript", "op_type")
                        numbers = pointer(script, "GapEditScript", "num")
                        int_size = types["GapEditScript"]["fields"]["size"]["size"]
                        row["gap_info"] = {"size": size,
                                           "op_type": [int.from_bytes(memory(operations + i * layout["op_size"], layout["op_size"]), "little", signed=True) for i in range(size)],
                                           "num": [int.from_bytes(memory(numbers + i * int_size, int_size), "little", signed=True) for i in range(size)]}
                    rows.append(row)
        return rows

    def save_state():
        Path(config["status"]).write_text(json.dumps(state, indent=2) + "\n")

    def caller_names():
        result = []
        frame = gdb.newest_frame()
        for _ in range(12):
            if frame is None:
                break
            result.append(frame.name())
            frame = frame.older()
        return result

    def check_architecture():
        architecture = gdb.newest_frame().architecture().name()
        if architecture != "i386:x86-64":
            raise RuntimeError(f"unsupported inferior architecture {architecture}")

    class CompletedSort(gdb.FinishBreakpoint):
        def __init__(self, address, callers):
            super().__init__(gdb.newest_frame(), internal=True)
            self.address, self.callers = address, callers

        def stop(self):
            try:
                state["sort_returns"].append({"results_address": hex(self.address),
                                               "num_queries": integer(self.address, "BlastHSPResults", "num_queries"),
                                               "callers": self.callers})
            except Exception as error:
                state["errors"].append(str(error))
            save_state()
            return False

    class ResultsSort(gdb.Breakpoint):
        def stop(self):
            try:
                check_architecture()
                address = int(gdb.parse_and_eval("$rdi"))
                if address:
                    CompletedSort(address, caller_names())
            except Exception as error:
                state["errors"].append(str(error))
                save_state()
            return False

    class CompletedBatch(gdb.FinishBreakpoint):
        def __init__(self, output_pointer, entry):
            super().__init__(gdb.newest_frame(), internal=True)
            self.output_pointer, self.entry = output_pointer, entry

        def stop(self):
            try:
                status = int(gdb.parse_and_eval("$rax")) & 0xffff
                if status:
                    raise ValueError(f"BLAST_ComputeTraceback_MT returned status {status}")
                address = int.from_bytes(memory(self.output_pointer, layout["pointer_size"]), "little")
                if not address:
                    raise ValueError("successful traceback returned NULL results")
                batch_rows = capture(address)
                completed_batches.append({"entry": self.entry, "rows": batch_rows})
                rows, mapping = aggregate_batches(completed_batches, nucleotide=config["contexts_per_query"] == 2)
                temporary = Path(config["output"] + ".partial")
                temporary.write_text("".join(json.dumps(row, allow_nan=False, separators=(",", ":")) + "\n" for row in rows))
                temporary.replace(config["output"])
                state["snapshots"] = len(completed_batches)
                state["batches"] = mapping
                state["num_queries"] = rows[0]["num_queries"]
                state["hsp_count"] = sum(row["kind"] == "hsp" for row in rows)
                state["last_results_address"] = hex(address)
                state["last_num_queries"] = batch_rows[0]["num_queries"]
                state["last_hsp_count"] = sum(row["kind"] == "hsp" for row in batch_rows)
            except Exception as error:
                state["errors"].append(str(error))
            save_state()
            return False

        def out_of_scope(self):
            state["errors"].append(f"traceback batch {self.entry['batch_index']} did not return normally")
            save_state()

    class TracebackBatch(gdb.Breakpoint):
        def stop(self):
            nonlocal next_batch_index
            try:
                check_architecture()
                query_info = int(gdb.parse_and_eval("$rcx"))
                entry = {"batch_index": next_batch_index,
                         "program_number": int(gdb.parse_and_eval("$rdi")) & 0xffffffff,
                         "num_queries": integer(query_info, "BlastQueryInfo", "num_queries"),
                         "first_context": integer(query_info, "BlastQueryInfo", "first_context"),
                         "last_context": integer(query_info, "BlastQueryInfo", "last_context"),
                         "callers": caller_names()}
                next_batch_index += 1
                # SysV AMD64: first six args in registers, arg10 at rsp+32.
                # Pinned blast_traceback.c defines results_out as arg10.
                stack = int(gdb.parse_and_eval("$rsp"))
                output_pointer = int.from_bytes(memory(stack + 32, layout["pointer_size"]), "little")
                if entry["first_context"] != 0 or entry["last_context"] + 1 != entry["num_queries"] * config["contexts_per_query"]:
                    raise ValueError("nonstandard query contexts require an explicit batch mapping")
                CompletedBatch(output_pointer, entry)
            except Exception as error:
                state["errors"].append(str(error))
                save_state()
            return False

    class TracebackSeeds(gdb.Breakpoint):
        def stop(self):
            try:
                subject = int(gdb.parse_and_eval("$rsi"))
                if not subject:
                    return False
                oid = integer(subject, "BlastHSPList", "oid")
                count = integer(subject, "BlastHSPList", "hspcnt")
                if not 0 <= count <= 10000000:
                    raise ValueError(f"invalid traceback input count {count}")
                array = pointer(subject, "BlastHSPList", "hsp_array")
                rows = []
                for i in range(count):
                    hsp = pointer_at(array, i)
                    if not hsp:
                        raise ValueError("NULL traceback input HSP")
                    row = {"kind": "traceback_seed", "oid": oid, "hsp_index": i,
                           "list_query_index": integer(subject, "BlastHSPList", "query_index"),
                           "context": integer(hsp, "BlastHSP", "context"),
                           "score": integer(hsp, "BlastHSP", "score")}
                    for name in ("query", "subject"):
                        row[name] = segment(hsp + types["BlastHSP"]["fields"][name]["offset"])
                    rows.append(row)
                with Path(config["traceback_seeds"]).open("a") as output:
                    output.write("".join(json.dumps(row, separators=(",", ":")) + "\n" for row in rows))
                state["traceback_seed_lists"] = state.get("traceback_seed_lists", 0) + 1
            except Exception as error:
                state["errors"].append(str(error))
            save_state()
            return False

    class GappedInputs(gdb.Breakpoint):
        def stop(self):
            try:
                subject = int(gdb.parse_and_eval("$rcx"))
                oid = integer(subject, "BLAST_SequenceBlk", "oid")
                if config.get("subject_oid") is not None and oid != config["subject_oid"]:
                    return False
                stack = int(gdb.parse_and_eval("$rsp"))
                stack_arg = lambda offset: int.from_bytes(memory(stack + offset, layout["pointer_size"]), "little")
                word, initial, extension = stack_arg(24), stack_arg(32), stack_arg(8)
                query_info = int(gdb.parse_and_eval("$rdx"))
                options = pointer(word, "BlastInitialWordParameters", "options")
                total = integer(initial, "BlastInitHitList", "total")
                if not 0 <= total <= 10000000:
                    raise ValueError(f"invalid initial hit count {total}")
                rows = [{"kind": "gapped_input", "oid": oid, "total": total,
                         "subject_length": integer(subject, "BLAST_SequenceBlk", "length"),
                         "subject_frame": integer(subject, "BLAST_SequenceBlk", "frame"),
                         "word": {name: integer(word, "BlastInitialWordParameters", name) for name in ("x_dropoff_max", "cutoff_score_min", "container_type")},
                         "options": {name: integer(options, "BlastInitialWordOptions", name) for name in ("window_size", "scan_range")},
                         "extension": {name: integer(extension, "BlastExtensionParameters", name) for name in ("gap_x_dropoff", "gap_x_dropoff_final")}}]
                for name in ("gap_trigger", "x_dropoff"):
                    rows[0]["options"][name], rows[0]["options"][name + "_bits"] = floating(options, "BlastInitialWordOptions", name)
                contexts = pointer(query_info, "BlastQueryInfo", "contexts")
                cutoffs = pointer(word, "BlastInitialWordParameters", "cutoffs")
                first, last = [integer(query_info, "BlastQueryInfo", key) for key in ("first_context", "last_context")]
                for c in range(first, last + 1):
                    context = contexts + c * types["BlastContextInfo"]["size"]
                    cutoff = cutoffs + c * types["BlastUngappedCutoffs"]["size"]
                    rows.append({"kind": "gapped_context", "oid": oid, "context": c,
                                 "query": {name: integer(context, "BlastContextInfo", name) for name in ("query_offset", "query_length", "frame", "is_valid")},
                                 "cutoff": {name: integer(cutoff, "BlastUngappedCutoffs", name) for name in ("x_dropoff_init", "x_dropoff", "cutoff_score", "reduced_nucl_cutoff_score")}})
                array = pointer(initial, "BlastInitHitList", "init_hsp_array")
                for i in range(total):
                    hit = array + i * types["BlastInitHSP"]["size"]
                    pair = hit + types["BlastInitHSP"]["fields"]["offsets"]["offset"]
                    ungapped = pointer(hit, "BlastInitHSP", "ungapped_data")
                    rows.append({"kind": "initial_hit", "oid": oid, "index": i,
                                 "offsets": {name: integer(pair, "BlastOffsetPair", "qs_offsets." + name, False) for name in ("q_off", "s_off")},
                                 "ungapped": None if not ungapped else {name: integer(ungapped, "BlastUngappedData", name) for name in ("q_start", "s_start", "length", "score")}})
                with Path(config["gapped_inputs"]).open("a") as output:
                    output.write("".join(json.dumps(row, separators=(",", ":")) + "\n" for row in rows))
                state["gapped_input_lists"] = state.get("gapped_input_lists", 0) + 1
            except Exception as error:
                state["errors"].append(str(error))
            save_state()
            return False

    def exited(event):
        state["inferior_exit_code"] = getattr(event, "exit_code", None)
        save_state()

    gdb.events.exited.connect(exited)
    TracebackBatch("BLAST_ComputeTraceback_MT", internal=True)
    ResultsSort("Blast_HSPResultsSortByEvalue", internal=True)
    if config.get("traceback_seeds"):
        TracebackSeeds("Blast_TracebackFromHSPList", internal=True)
    if config.get("gapped_inputs"):
        GappedInputs("BLAST_GetGappedScore", internal=True)
    save_state()


def host_main():
    import argparse
    import hashlib
    import platform
    import subprocess
    import tempfile

    repo = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--program", required=True, choices=("blastn", "megablast", "tblastn"))
    parser.add_argument("--database", required=True, type=Path, help="input nucleotide FASTA")
    parser.add_argument("--queries", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path, help="ordered raw core JSONL")
    parser.add_argument("--metadata", type=Path)
    parser.add_argument("--work-dir", type=Path)
    parser.add_argument("--traceback-seeds", type=Path, help="optional original preliminary HSPs at traceback entry")
    parser.add_argument("--gapped-inputs", type=Path, help="optional original ungapped hits/options at GetGappedScore entry")
    parser.add_argument("--subject-oid", type=int, help="filter optional gapped-input capture to one subject")
    parser.add_argument("--cc", default="cc")
    parser.add_argument("--gdb", default="gdb")
    args = parser.parse_args()
    if platform.system() != "Linux" or platform.machine() != "x86_64":
        raise ValueError("official pinned oracle executable is Linux x86-64")
    subprocess.run([sys.executable, str(repo / "scripts/blast_tools.py"), "verify"], check=True, stdout=subprocess.DEVNULL)
    pin = json.loads((repo / "upstream/blast-tools.json").read_text())
    source = json.loads((repo / "upstream/blast-source.json").read_text())
    if pin["version"] != source["version"]:
        raise ValueError("source and executable pins differ")
    work = args.work_dir.resolve() if args.work_dir else Path(tempfile.mkdtemp(prefix="btyper-core-oracle-"))
    work.mkdir(parents=True, exist_ok=True)
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    metadata = args.metadata.resolve() if args.metadata else Path(str(output) + ".metadata.json")
    executable = repo / pin["root"] / ("tblastn" if args.program == "tblastn" else "blastn")
    env = dict(os.environ, BLAST_USAGE_REPORT="false", NCBI_DONT_USE_NCBIRC="1", NCBI_DONT_USE_LOCAL_CONFIG="1")
    version = subprocess.run([str(executable), "-version"], check=True, capture_output=True, text=True, env=env).stdout
    if pin["version"] not in version:
        raise ValueError("reference executable version differs from pin")
    reporter = work / "layout-reporter"
    compile_command = [args.cc, "-std=c11", "-I" + str(repo / "tests/reference/include"),
                       "-I" + str(repo / source["root"] / "include"), str(repo / "scripts/core_result_layout.c"), "-o", str(reporter)]
    subprocess.run(compile_command, check=True)
    layout = json.loads(subprocess.run([str(reporter)], check=True, capture_output=True, text=True).stdout)
    query_ids = [line[1:].split()[0] for line in args.queries.read_text().splitlines() if line.startswith(">")]
    config = {"layout": layout, "output": str(output), "status": str(work / "capture-status.json"),
              "contexts_per_query": 1 if args.program == "tblastn" else 2}
    if args.traceback_seeds:
        seeds = args.traceback_seeds.resolve()
        seeds.parent.mkdir(parents=True, exist_ok=True)
        seeds.write_text("")
        config["traceback_seeds"] = str(seeds)
    if args.gapped_inputs:
        path = args.gapped_inputs.resolve()
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("")
        config["gapped_inputs"] = str(path)
        config["subject_oid"] = args.subject_oid
    output.unlink(missing_ok=True)
    config_path = work / "capture-config.json"
    config_path.write_text(json.dumps(config, indent=2) + "\n")
    env["BTYPER_CORE_ORACLE_CONFIG"] = str(config_path)
    database = work / "db"
    build_command = [str(repo / pin["root"] / "makeblastdb"), "-in", str(args.database.resolve()), "-dbtype", "nucl", "-blastdb_version", "5", "-out", str(database)]
    build = subprocess.run(build_command, check=True, capture_output=True, text=True, cwd=work, env=env)
    (work / "makeblastdb.log").write_text(build.stdout + build.stderr)
    command = [str(executable), "-db", str(database), "-query", str(args.queries.resolve()), "-num_threads", "1", "-outfmt", "6"]
    if args.program == "blastn":
        command += ["-task", "blastn"]
    # A GDB script uses Python repr for the path, so spaces/quotes are harmless.
    debugger_script = work / "capture.gdb"
    debugger_script.write_text("set pagination off\nset confirm off\nset breakpoint pending on\n"
                              "set debuginfod enabled off\nset print thread-events off\n"
                              "python\nexec(compile(open(" + repr(str(Path(__file__).resolve())) + ").read(), " + repr(str(Path(__file__).resolve())) + ", 'exec'))\nend\nrun\n")
    debugger_command = [args.gdb, "--batch", "--nx", "-x", str(debugger_script), "--args"] + command
    completed = subprocess.run(debugger_command, cwd=work, env=env, capture_output=True, text=True)
    (work / "gdb.stdout").write_text(completed.stdout)
    (work / "gdb.stderr").write_text(completed.stderr)
    status_path = Path(config["status"])
    status = json.loads(status_path.read_text()) if status_path.exists() else {}
    result = {"capture_schema_version": 2, "program": args.program,
              "reference_version": pin["version"], "archive_sha256": pin["sha256"],
              "executable": str(executable), "executable_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
              "query_ids": query_ids, "query_batch_mapping": status.get("batches", []),
              "query_index_scope": "query_index/context are global; local_query_index/local_context preserve each batch; list_query_index remains the unmodified source field",
              "list_query_index_scope": "raw source field is unchanged; nucleotide subject record global_list_query_index explicitly adds recorded query_index_offset; composition-redo default0 is not rebased",
              "oracle_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "inputs": {name: {"path": str(path.resolve()), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()} for name, path in (("database", args.database), ("queries", args.queries))},
              "layout": layout, "layout_compile_command": compile_command,
              "capture_boundary_source": "src/algo/blast/core/blast_traceback.c:BLAST_ComputeTraceback_MT",
              "capture_boundary_source_sha256": hashlib.sha256((repo / source["root"] / "src/algo/blast/core/blast_traceback.c").read_bytes()).hexdigest(),
              "layout_reporter_sha256": hashlib.sha256((repo / "scripts/core_result_layout.c").read_bytes()).hexdigest(),
              "layout_header_sha256": {name: hashlib.sha256((repo / source["root"] / "include/algo/blast/core" / name).read_bytes()).hexdigest() for name in ("blast_hits.h", "gapinfo.h", "ncbi_std.h", "blast_extend.h", "blast_parameters.h", "blast_options.h", "blast_query_info.h", "blast_def.h")}, "debugger_command": debugger_command,
              "debugger_exit_status": completed.returncode, "work_dir": str(work), "capture": status,
              "precision": "original core double bits before Seq-align Evalue<1e-180 clamping; storage fields are diagnostic"}
    metadata.parent.mkdir(parents=True, exist_ok=True)
    metadata.write_text(json.dumps(result, indent=2) + "\n")
    if completed.returncode or status.get("inferior_exit_code") != 0 or not status.get("snapshots") or status.get("errors") or status.get("num_queries") != len(query_ids):
        raise RuntimeError(f"core capture failed; inspect {work / 'gdb.stderr'} and {metadata}")
    print(f"{status['hsp_count']} core HSPs across {status['num_queries']} queries captured from {status['snapshots']} final batch return(s): {output}")


try:
    import gdb
except ImportError:
    if __name__ == "__main__":
        host_main()
else:
    debugger_main(gdb)
