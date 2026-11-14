#!/usr/bin/env python3
"""Parity-gated paired FASTA-to-results benchmarks against pinned BLAST+.

This measures end-to-end user workflows, not isolated search engines. Original:
FASTA -> fresh on-disk v5 database -> BLAST CLI -> captured 19-field output.
Native: FASTA -> in-memory database -> native CLI -> captured same output.
Executable verification, parity comparisons, temporary-directory creation/deletion,
and report generation are outside measured intervals. No cache flushing occurs.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import random
import shutil
import statistics
import subprocess
import sys
import tempfile
import time

from compare_search import compare, PROGRAMS
from reference_search import FIELDS

try:
    import resource
except ImportError:  # Windows: CPU child accounting is unavailable in stdlib.
    resource = None

REPO = Path(__file__).resolve().parents[1]


def file_metadata(path):
    path = path.resolve()
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": digest.hexdigest()}


def cpu_usage():
    if resource is None:
        return None
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    return usage.ru_utime, usage.ru_stime


def command(command_line, env, cwd=None):
    cpu_before = cpu_usage()
    start = time.perf_counter_ns()
    process = subprocess.run(command_line, capture_output=True, check=False, env=env, cwd=cwd)
    elapsed = (time.perf_counter_ns() - start) / 1e9
    cpu_after = cpu_usage()
    if process.returncode:
        raise RuntimeError(f"command failed ({process.returncode}): {command_line!r}\n"
                           + process.stderr.decode("utf-8", errors="replace")[-12000:])
    return {"wall_seconds": elapsed,
            "cpu_user_seconds": None if cpu_before is None else cpu_after[0] - cpu_before[0],
            "cpu_system_seconds": None if cpu_before is None else cpu_after[1] - cpu_before[1],
            "stdout_bytes": len(process.stdout), "stderr_bytes": len(process.stderr)}, process.stdout, process.stderr


def rows(output):
    result = [line.split("\t") for line in output.decode("utf-8").splitlines()]
    for index, row in enumerate(result):
        if len(row) != len(FIELDS):
            raise ValueError(f"row {index}: {len(row)} fields, expected {len(FIELDS)}")
    return result


def original_hits(output):
    return [dict(zip(FIELDS, row)) for row in rows(output)]


def percentile(sorted_values, fraction):
    position = fraction * (len(sorted_values) - 1)
    lower = int(position)
    upper = min(lower + 1, len(sorted_values) - 1)
    return sorted_values[lower] + (position - lower) * (sorted_values[upper] - sorted_values[lower])


def paired_summary(pairs, seed, samples):
    original = [pair["reference"]["wall_seconds"] for pair in pairs]
    native = [pair["native"]["wall_seconds"] for pair in pairs]
    rng = random.Random(seed)
    ratios = []
    for _ in range(samples):
        indexes = [rng.randrange(len(pairs)) for _ in pairs]
        ratios.append(sum(native[i] for i in indexes) / sum(original[i] for i in indexes))
    ratios.sort()
    low, high = percentile(ratios, 0.025), percentile(ratios, 0.975)
    ratio = statistics.mean(native) / statistics.mean(original)
    return {"reference_mean_seconds": statistics.mean(original),
            "reference_median_seconds": statistics.median(original),
            "native_mean_seconds": statistics.mean(native),
            "native_median_seconds": statistics.median(native),
            "native_over_reference_ratio_of_means": ratio,
            "reference_over_native_speedup_ratio_of_means": 1.0 / ratio,
            "paired_bootstrap_95pct_native_over_reference": [low, high],
            "bootstrap_method": "percentile interval; resample whole repetition pairs with replacement; ratio of arithmetic mean wall times",
            "bootstrap_samples": samples, "bootstrap_seed": seed,
            "improvement_criterion": "at least 5 measured pairs and 95% CI upper bound of native/reference < 1",
            "statistically_supported_workflow_improvement": len(pairs) >= 5 and high < 1.0,
            "reference_build_mean_seconds": statistics.mean(pair["reference"]["build"]["wall_seconds"] for pair in pairs),
            "reference_search_mean_seconds": statistics.mean(pair["reference"]["search"]["wall_seconds"] for pair in pairs)}


def git_output(*args):
    try:
        return subprocess.run(["git", *args], cwd=REPO, capture_output=True,
                              text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def source_metadata():
    # Include untracked translation/native adapter files: a commit identifier
    # alone is insufficient while this translation is being assembled.
    paths = [REPO / "Cargo.toml", REPO / "Cargo.lock"]
    for directory in ("src", "examples"):
        paths.extend(sorted((REPO / directory).rglob("*.rs")))
    records = [{"path": str(path.relative_to(REPO)), "sha256": file_metadata(path)["sha256"]}
               for path in paths if path.is_file()]
    canonical = json.dumps(records, sort_keys=True, separators=(",", ":")).encode()
    return {"git_commit": git_output("rev-parse", "HEAD"),
            "git_status": git_output("status", "--porcelain"),
            "rust_source_and_cargo_sha256": hashlib.sha256(canonical).hexdigest(),
            "files": records}


def cpu_model():
    cpuinfo = Path("/proc/cpuinfo")
    if cpuinfo.is_file():
        for line in cpuinfo.read_text(errors="replace").splitlines():
            if line.startswith(("model name", "Hardware")) and ":" in line:
                return line.split(":", 1)[1].strip()
    return os.environ.get("PROCESSOR_IDENTIFIER") or platform.processor()


def tool_version(name):
    executable = shutil.which(name)
    if executable is None:
        return None
    return subprocess.run([executable, "--version"], capture_output=True,
                          text=True, check=True).stdout.strip()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native", required=True, type=Path)
    parser.add_argument("--database", required=True, type=Path)
    parser.add_argument("--nucleotide-queries", type=Path)
    parser.add_argument("--protein-queries", type=Path)
    parser.add_argument("--program", action="append", choices=PROGRAMS,
                        help="Repeat to select presets; default: all three")
    parser.add_argument("--repetitions", type=int, default=7)
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--seed", type=int, default=217)
    parser.add_argument("--bootstrap-samples", type=int, default=10000)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--manifest", type=Path, help="Optional workload provenance JSON; recorded verbatim, not treated as a representativeness claim")
    parser.add_argument("--build-metadata", type=Path, help="Optional native build provenance JSON, e.g. exact cargo command/profile/RUSTFLAGS/target")
    parser.add_argument("--temporary-directory", type=Path,
                        help="Parent directory for fresh reference v5 databases; record its filesystem/storage when interpreting results")
    args = parser.parse_args()
    if args.repetitions < 1 or args.warmups < 0 or args.bootstrap_samples < 100:
        parser.error("repetitions must be >=1, warmups >=0, bootstrap-samples >=100")
    programs = list(dict.fromkeys(args.program or PROGRAMS))
    for program in programs:
        if (args.protein_queries if program == "tblastn" else args.nucleotide_queries) is None:
            parser.error(f"{program} requires --{'protein' if program == 'tblastn' else 'nucleotide'}-queries")
    native = args.native.resolve()
    if not native.exists() and native.with_suffix(".exe").exists():
        native = native.with_suffix(".exe")
    database = args.database.resolve()
    queries = {program: (args.protein_queries if program == "tblastn" else args.nucleotide_queries).resolve()
               for program in programs}
    # Verify the official pinned binaries before any measurement.
    subprocess.run([sys.executable, str(REPO / "scripts/blast_tools.py"), "verify"],
                   check=True, stdout=subprocess.DEVNULL)
    pin = json.loads((REPO / "upstream/blast-tools.json").read_text())
    source_pin = json.loads((REPO / "upstream/blast-source.json").read_text())
    if pin["version"] != source_pin["version"]:
        raise ValueError("reference and translated source pins differ")
    binaries = (REPO / pin["root"]).resolve()
    env = dict(os.environ, BLAST_USAGE_REPORT="false", NCBI_DONT_USE_NCBIRC="1",
               NCBI_DONT_USE_LOCAL_CONFIG="1", OMP_NUM_THREADS="1",
               OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", RAYON_NUM_THREADS="1", LC_ALL="C")
    report = {"schema_version": 1, "started_utc": datetime.now(timezone.utc).isoformat(),
              "scope": "end-to-end FASTA-to-captured-results user workflow; not engine-only speed",
              "interpretation": "Results apply only to the supplied inputs, artifact, settings and host. No workload representativeness or engine-only speed is inferred from this workflow comparison.",
              "workflow_difference": "reference uses two processes and fresh on-disk v5 storage; native uses one process and in-memory storage",
              "timing_boundary": "workflow wall: immediately before build/native subprocess through complete search/native captured stdout+stderr; excludes temporary directory creation/deletion, pin verification, parity checks and reporting",
              "cache_policy": "warm-cache workflow measurements; no filesystem/page-cache eviction or flushing",
              "output_policy": "both sides capture all 19 fields into Python memory on every run; compare outside timed intervals",
              "precision": "exact ordered non-floating fields; compare_search.py BLAST CLI precision rules for Evalue/bit score",
              "reference_version": pin["version"], "reference_archive_sha256": pin["sha256"],
              "translated_source_pin": source_pin, "threads_per_search": 1,
              "warmups_per_program": args.warmups, "measured_pairs_per_program": args.repetitions,
              "random_seed": args.seed, "platform": platform.platform(),
              "python": sys.version, "cpu": cpu_model(), "machine": platform.machine(), "logical_cpus": os.cpu_count(),
              "load_average_at_start": list(os.getloadavg()) if hasattr(os, "getloadavg") else None,
              "clock": {"name": "perf_counter_ns", "resolution_seconds": time.get_clock_info("perf_counter").resolution},
              "harness": {name: file_metadata(REPO / "scripts" / name) for name in ("benchmark_search.py", "compare_search.py", "reference_search.py", "blast_tools.py")},
              "cpu_affinity": sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else None,
              "cpu_accounting": "POSIX RUSAGE_CHILDREN user/system deltas" if resource else "unavailable on this platform",
              "temporary_parent": str(args.temporary_directory.resolve()) if args.temporary_directory else tempfile.gettempdir(),
              "inputs": {"database": file_metadata(database),
                         **{program + "_queries": file_metadata(path) for program, path in queries.items()}},
              "executables": {"native": file_metadata(native),
                              **{name: file_metadata(binaries / name) for name in ("makeblastdb", "blastn", "tblastn", "blastdbcmd")}},
              "source": source_metadata(), "observed_rustc": tool_version("rustc"), "observed_cargo": tool_version("cargo"),
              "build_provenance_note": "observed toolchain/source describe this invocation, not necessarily the supplied native executable; exact build attribution requires --build-metadata",
              "build_metadata": None, "workload_manifest": None,
              "environment": {key: env.get(key) for key in ("BLAST_USAGE_REPORT", "NCBI_DONT_USE_NCBIRC", "NCBI_DONT_USE_LOCAL_CONFIG", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "RAYON_NUM_THREADS", "LC_ALL", "RUSTFLAGS", "CARGO_BUILD_TARGET")},
              "searches": {}}
    for argument, key in ((args.manifest, "workload_manifest"), (args.build_metadata, "build_metadata")):
        if argument:
            report[key] = {"file": file_metadata(argument), "content": json.loads(argument.read_text())}
    reference_commands = {}
    native_commands = {}
    def reference_workflow(program, verify_format=False):
        with tempfile.TemporaryDirectory(prefix="btyper-paired-", dir=args.temporary_directory) as temporary:
            db = Path(temporary) / "database"
            build_command = [str(binaries / "makeblastdb"), "-in", str(database), "-dbtype", "nucl", "-blastdb_version", "5", "-out", str(db)]
            executable = "tblastn" if program == "tblastn" else "blastn"
            search_command = [str(binaries / executable), "-query", str(queries[program]), "-db", str(db), "-num_threads", "1", "-outfmt", "6 " + " ".join(FIELDS)]
            if program == "blastn":
                search_command.extend(["-task", "blastn"])
            started = time.perf_counter_ns()
            build, _, build_stderr = command(build_command, env, temporary)
            search, output, search_stderr = command(search_command, env, temporary)
            total = (time.perf_counter_ns() - started) / 1e9
            db_bytes = sum(path.stat().st_size for path in Path(temporary).iterdir() if path.is_file())
            if verify_format:
                info = subprocess.run([str(binaries / "blastdbcmd"), "-db", str(db), "-info"], capture_output=True, check=True, env=env)
                if b"BLASTDB Version: 5" not in info.stdout:
                    raise ValueError("reference database is not v5")
            reference_commands[program] = {"build": [value.replace(str(db), "<temporary-db>") for value in build_command],
                                           "search": [value.replace(str(db), "<temporary-db>") for value in search_command]}
            return {"wall_seconds": total, "build": build, "search": search,
                    "cpu_user_seconds": None if build["cpu_user_seconds"] is None else build["cpu_user_seconds"] + search["cpu_user_seconds"],
                    "cpu_system_seconds": None if build["cpu_system_seconds"] is None else build["cpu_system_seconds"] + search["cpu_system_seconds"],
                    "database_bytes_on_disk": db_bytes,
                    "output_bytes": len(output), "output_sha256": hashlib.sha256(output).hexdigest(),
                    "diagnostics": (build_stderr + search_stderr).decode("utf-8", errors="replace")}, output
    def native_workflow(program):
        line = [str(native), program, str(database), str(queries[program])]
        native_commands[program] = line
        started = time.perf_counter_ns()
        metrics, output, diagnostics = command(line, env)
        total = (time.perf_counter_ns() - started) / 1e9
        return {**metrics, "process_wall_seconds": metrics["wall_seconds"], "wall_seconds": total, "output_bytes": len(output), "output_sha256": hashlib.sha256(output).hexdigest(),
                "diagnostics": diagnostics.decode("utf-8", errors="replace")}, output
    rng = random.Random(args.seed)
    def check_native_output(program, expected_output, native_output, label):
        try:
            compare(rows(native_output), original_hits(expected_output), program)
        except ValueError:
            directory = Path(str(args.report) + ".failure") / program / label
            directory.mkdir(parents=True, exist_ok=True)
            original_path = directory / "original.tsv"
            native_path = directory / "native.tsv"
            original_path.write_bytes(expected_output)
            native_path.write_bytes(native_output)
            report.setdefault("parity_failures", []).append({
                "program": program, "stage": label,
                "original_output": file_metadata(original_path),
                "native_output": file_metadata(native_path),
            })
            raise
    try:
        for program_index, program in enumerate(programs):
            print(f"{program}: parity preflight", file=sys.stderr, flush=True)
            reference, expected_output = reference_workflow(program, verify_format=True)
            native_result, native_output = native_workflow(program)
            check_native_output(program, expected_output, native_output, "preflight")
            preflight = {"passed": True, "hits": len(rows(expected_output)),
                         "reference_output_bytes": len(expected_output), "native_output_bytes": len(native_output),
                         "reference_output_sha256": reference["output_sha256"], "native_output_sha256": native_result["output_sha256"]}
            warmups, pairs = [], []
            for index in range(args.warmups + args.repetitions):
                order = ["reference", "native"]
                rng.shuffle(order)
                pair = {"index": index - args.warmups if index >= args.warmups else index,
                        "order": order, "warmup": index < args.warmups}
                for implementation in order:
                    metrics, output = reference_workflow(program) if implementation == "reference" else native_workflow(program)
                    # Recheck every output outside the timed interval. Native
                    # output floats are intentionally not compared as raw text.
                    if implementation == "reference":
                        if rows(output) != rows(expected_output):
                            raise ValueError(program + ": reference output changed between runs")
                    else:
                        check_native_output(program, expected_output, output, f"trial-{index}")
                    pair[implementation] = metrics
                (warmups if pair["warmup"] else pairs).append(pair)
                print(f"{program}: {'warmup' if pair['warmup'] else 'pair'} {pair['index'] + 1} complete", file=sys.stderr, flush=True)
            report["searches"][program] = {"parity": preflight, "reference_commands": reference_commands[program],
                                          "native_command": native_commands[program], "warmups": warmups,
                                          "pairs": pairs, "summary": paired_summary(pairs, args.seed + program_index + 1, args.bootstrap_samples)}
        report["status"] = "complete"
    except Exception as error:
        report["status"] = "failed"
        report["error"] = str(error)
        raise
    finally:
        report["finished_utc"] = datetime.now(timezone.utc).isoformat()
        report["rust_source_changed_during_run"] = source_metadata()["rust_source_and_cargo_sha256"] != report["source"]["rust_source_and_cargo_sha256"]
        report["native_executable_changed_during_run"] = file_metadata(native)["sha256"] != report["executables"]["native"]["sha256"]
        report["inputs_changed_during_run"] = any(
            file_metadata(Path(entry["path"]))["sha256"] != entry["sha256"]
            for entry in report["inputs"].values())
        build = report["build_metadata"]
        report["provided_build_artifact_matches"] = (
            build["content"].get("native_sha256") == report["executables"]["native"]["sha256"]
            if build is not None else None)
        report["measurement_integrity_passed"] = (
            report["status"] == "complete"
            and not report["rust_source_changed_during_run"]
            and not report["native_executable_changed_during_run"]
            and not report["inputs_changed_during_run"]
            and report["provided_build_artifact_matches"] is not False)
        if not report["measurement_integrity_passed"]:
            for search in report["searches"].values():
                search["summary"]["statistically_supported_workflow_improvement"] = False
                search["summary"]["invalidated_reason"] = "Whole selected benchmark failed parity or measurement-integrity checks; timings are diagnostic only"
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(str(args.report.resolve()))


if __name__ == "__main__":
    main()
