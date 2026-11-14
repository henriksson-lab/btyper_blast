#!/usr/bin/env python3
"""Summarize the complete fixed 2-dataset x 3-preset paired workflow suite.

The suite observation is one invocation of every fixed workload: sum the six
native arithmetic mean wall times and divide by the sum of the six reference
means. It is not a mean of speedups, an engine-only benchmark, or a claim about
other inputs/hosts. Each reference workflow includes fresh v5 database creation
and a search process; each native workflow includes in-memory construction and
search in one process. Both capture all output. These intentionally different
storage/process workflows are preserved from benchmark_search.py.

A deterministic joint paired bootstrap (10,000 draws, seed217) independently
resamples repetition-pair indexes within each case. Each draw supplies every
case's ratio and the ratio of the sums of resampled case means. Percentile95%
intervals classify CIupper<1 as a workflow win, CIlower>1 as a regression, and
all other intervals as inconclusive. Every case needs at least five pairs.
Warmups are excluded as declared by the original experiment; no measured trial,
outlier, slow case or unfavorable result is removed. Complete original reports
and raw observations are embedded in the summary for audit/recalculation.

Only complete, integrity-passing reports with compatible artifact/build/source,
thread/environment/host settings and pinned input manifests are accepted. The
script cannot independently prove that metadata or underlying measurements are
truthful; it checks their recorded consistency and preserves their scope.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import random
import statistics

REPO = Path(__file__).resolve().parents[1]
DATASETS = ("staph_st14", "yeast_s288c")
PROGRAMS = ("megablast", "blastn", "tblastn")
BOOTSTRAP_SAMPLES = 10_000
BOOTSTRAP_SEED = 217
MINIMUM_PAIRS = 5


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(value, label):
    require(isinstance(value, str) and len(value) == 64
            and all(char in "0123456789abcdef" for char in value), f"invalid {label} SHA256")
    return value


def positive_time(value, label):
    require(type(value) in (int, float) and math.isfinite(value) and value > 0,
            f"invalid positive wall time: {label}={value!r}")
    return value


def metadata_identity(value, label):
    require(isinstance(value, dict), f"missing {label} metadata")
    require(type(value.get("bytes")) is int and value["bytes"] > 0,
            f"invalid {label} byte count")
    return {"bytes": value["bytes"], "sha256": digest(value.get("sha256"), label)}


def manifest_inputs(report, dataset):
    manifest = report.get("workload_manifest")
    require(isinstance(manifest, dict), "missing pinned workload manifest")
    content = manifest.get("content")
    require(isinstance(content, dict) and content.get("dataset") == dataset,
            f"invalid workload manifest dataset {dataset}")
    recorded = metadata_identity(manifest.get("file"), "workload manifest")
    # Workload identities are fixed before looking at performance. Validate
    # against the repository's original per-dataset manifest, not report names.
    manifest_path = REPO / "vendor/benchmark-inputs" / dataset / "workload.json"
    original_bytes = manifest_path.read_bytes()
    require(recorded == {"bytes": len(original_bytes),
                         "sha256": hashlib.sha256(original_bytes).hexdigest()},
            f"{dataset}: workload manifest does not match the pinned local manifest")
    require(content == json.loads(original_bytes), f"{dataset}: changed manifest content")
    files = content.get("files", {})
    expected_names = {"database": "subjects.fna", "megablast_queries": "nucleotide-queries.fna",
                      "blastn_queries": "nucleotide-queries.fna", "tblastn_queries": "protein-queries.faa"}
    inputs = report.get("inputs", {})
    require(set(inputs) == set(expected_names), f"{dataset}: unexpected/missing input set")
    checked = {}
    for key, name in expected_names.items():
        identity = metadata_identity(inputs[key], f"{dataset}/{key}")
        require(identity == metadata_identity(files.get(name), f"manifest/{name}"),
                f"{dataset}: {key} is not pinned by the workload manifest")
        checked[key] = identity
    provenance = content.get("provenance", {})
    require(provenance.get("canonical_sha256") == checked["database"]["sha256"],
            f"{dataset}: database and provenance SHA256 disagree")
    return recorded, checked


def validate_trials(search, report, case):
    parity = search.get("parity", {})
    require(parity.get("passed") is True, f"{case}: parity did not pass")
    for key in ("reference_output_sha256", "native_output_sha256"):
        digest(parity.get(key), f"{case}/{key}")
    pairs = search.get("pairs")
    require(isinstance(pairs, list) and len(pairs) >= MINIMUM_PAIRS,
            f"{case}: at least {MINIMUM_PAIRS} measured pairs required")
    require(type(report.get("measured_pairs_per_program")) is int
            and len(pairs) == report["measured_pairs_per_program"], f"{case}: pair count mismatch")
    warmups = search.get("warmups")
    require(isinstance(warmups, list) and type(report.get("warmups_per_program")) is int
            and len(warmups) == report["warmups_per_program"], f"{case}: warmup count mismatch")
    for trials, warmup in ((pairs, False), (warmups, True)):
        require(all(type(trial.get("index")) is int for trial in trials)
                and [trial.get("index") for trial in trials] == list(range(len(trials))),
                f"{case}: duplicate, missing or reordered trial indexes")
        for trial in trials:
            label = f"{case}/{'warmup' if warmup else 'pair'}{trial['index']}"
            require(trial.get("warmup") is warmup, f"{label}: incorrect warmup flag")
            require(trial.get("order") in (["reference", "native"], ["native", "reference"]),
                    f"{label}: invalid paired execution order")
            for implementation in ("native", "reference"):
                metrics = trial.get(implementation, {})
                positive_time(metrics.get("wall_seconds"), f"{label}/{implementation}")
                require(metrics.get("output_sha256") == parity[f"{implementation}_output_sha256"],
                        f"{label}: changed {implementation} output")
                expected_bytes = parity.get(f"{implementation}_output_bytes")
                require(type(expected_bytes) is int and expected_bytes >= 0
                        and metrics.get("output_bytes") == expected_bytes,
                        f"{label}: changed {implementation} output size")
    return pairs


def validate_report(path, tools_pin, source_pin):
    payload = path.read_bytes()
    report = json.loads(payload)
    require(type(report.get("schema_version")) is int and report["schema_version"] == 1, f"{path}: unsupported report schema")
    require(report.get("status") == "complete", f"{path}: report is not complete")
    require(report.get("measurement_integrity_passed") is True, f"{path}: measurement integrity failed")
    for flag in ("rust_source_changed_during_run", "native_executable_changed_during_run",
                 "inputs_changed_during_run"):
        require(report.get(flag) is False, f"{path}: {flag} must be false")
    require(report.get("provided_build_artifact_matches") is True,
            f"{path}: verified native build attribution required")
    require(report.get("reference_version") == tools_pin["version"] == source_pin["version"],
            f"{path}: wrong reference/source version")
    require(report.get("reference_archive_sha256") == tools_pin["sha256"],
            f"{path}: wrong reference archive")
    require(report.get("translated_source_pin") == source_pin, f"{path}: translated source pin mismatch")
    require(type(report.get("threads_per_search")) is int and report["threads_per_search"] == 1, f"{path}: expected single-thread search")
    for field in ("platform", "cpu", "machine", "cache_policy", "timing_boundary",
                  "workflow_difference", "output_policy", "precision", "temporary_parent"):
        require(isinstance(report.get(field), str) and report[field], f"{path}: missing {field}")
    require(isinstance(report.get("environment"), dict), f"{path}: missing execution environment")
    for field in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "RAYON_NUM_THREADS"):
        require(report["environment"].get(field) == "1", f"{path}: {field} is not pinned to1")
    searches = report.get("searches", {})
    require(set(searches) == set(PROGRAMS), f"{path}: all three presets are mandatory; no extra cases")
    manifest = report.get("workload_manifest", {})
    dataset = manifest.get("content", {}).get("dataset")
    require(dataset in DATASETS, f"{path}: unexpected/missing dataset {dataset!r}")
    manifest_identity, inputs = manifest_inputs(report, dataset)
    native = metadata_identity(report.get("executables", {}).get("native"), "native executable")
    source = report.get("source", {})
    source_hash = digest(source.get("rust_source_and_cargo_sha256"), "Rust source/Cargo")
    build = report.get("build_metadata")
    require(isinstance(build, dict) and isinstance(build.get("content"), dict), f"{path}: build metadata required")
    build_file = metadata_identity(build.get("file"), "build metadata")
    content = build["content"]
    require(content.get("native_sha256") == native["sha256"], f"{path}: build artifact hash mismatch")
    require(content.get("rust_source_and_cargo_sha256") == source_hash, f"{path}: build/source mismatch")
    # Some builder records additionally contain before/after source snapshots.
    # The final compact builder record instead pins its source/artifact hashes;
    # do not fabricate absent snapshots or confuse metadata with independent
    # proof of a build's execution. Validate every snapshot that is provided.
    snapshots = [value for key, value in content.items() if key.startswith("source_") and isinstance(value, dict)]
    require(all(snapshot.get("rust_source_and_cargo_sha256") == source_hash for snapshot in snapshots),
            f"{path}: build source changed")
    reference_executables = {name: metadata_identity(report.get("executables", {}).get(name), name)
                             for name in ("makeblastdb", "blastn", "tblastn", "blastdbcmd")}
    cases = {}
    for program in PROGRAMS:
        cases[program] = validate_trials(searches[program], report, f"{dataset}/{program}")
    shared = {"native": native, "source_sha256": source_hash, "build_metadata": build_file,
              "reference_executables": reference_executables,
              **{key: report.get(key) for key in (
                  "threads_per_search", "platform", "cpu", "machine", "cpu_affinity", "environment",
                  "cache_policy", "timing_boundary", "workflow_difference", "output_policy", "precision",
                  "clock", "temporary_parent", "observed_rustc", "observed_cargo")},
              "harness": {name: metadata_identity(value, name) for name, value in report.get("harness", {}).items()}}
    require(shared["harness"], f"{path}: benchmark harness attribution missing")
    require(report.get("random_seed") == BOOTSTRAP_SEED, f"{path}: unexpected experiment random seed")
    return {"dataset": dataset, "path": str(path.resolve()), "sha256": hashlib.sha256(payload).hexdigest(),
            "record": report, "shared_identity": shared, "manifest_identity": manifest_identity,
            "pinned_inputs": inputs, "cases": cases}


def percentile(values, fraction):
    position = fraction * (len(values) - 1)
    lower = int(position)
    upper = min(lower + 1, len(values) - 1)
    return values[lower] + (position - lower) * (values[upper] - values[lower])


def result_statistics(reference_mean, native_mean, bootstrap_ratios):
    bootstrap_ratios.sort()
    interval = [percentile(bootstrap_ratios, 0.025), percentile(bootstrap_ratios, 0.975)]
    classification = "win" if interval[1] < 1 else "regression" if interval[0] > 1 else "inconclusive"
    ratio = native_mean / reference_mean
    return {"reference_mean_wall_seconds": reference_mean, "native_mean_wall_seconds": native_mean,
            "native_over_reference_ratio_of_means": ratio,
            "reference_over_native_speedup_ratio_of_means": 1 / ratio,
            "paired_bootstrap_95pct_native_over_reference": interval,
            "paired_bootstrap_95pct_reference_over_native": [1 / interval[1], 1 / interval[0]],
            "classification": classification, "statistically_supported_workflow_improvement": classification == "win"}


def summarize(reports):
    require(len(reports) == len(DATASETS), "exactly two dataset reports are required")
    tools = json.loads((REPO / "upstream/blast-tools.json").read_text())
    source = json.loads((REPO / "upstream/blast-source.json").read_text())
    validated = [validate_report(path, tools, source) for path in reports]
    by_dataset = {item["dataset"]: item for item in validated}
    require(len(by_dataset) == len(validated), "duplicate dataset/cases are forbidden")
    require(set(by_dataset) == set(DATASETS), "the complete fixed two-dataset suite is required")
    ordered_reports = [by_dataset[dataset] for dataset in DATASETS]
    shared = ordered_reports[0]["shared_identity"]
    require(all(item["shared_identity"] == shared for item in ordered_reports),
            "artifact/build/source/reference/thread/environment/host/harness identities differ between reports")
    cases = []
    for item in ordered_reports:
        for program in PROGRAMS:
            pairs = item["cases"][program]
            cases.append({"dataset": item["dataset"], "program": program, "pairs": pairs,
                          "reference": [pair["reference"]["wall_seconds"] for pair in pairs],
                          "native": [pair["native"]["wall_seconds"] for pair in pairs],
                          "bootstrap": []})
    rng = random.Random(BOOTSTRAP_SEED)
    suite_bootstrap = []
    for _ in range(BOOTSTRAP_SAMPLES):
        reference_means, native_means = [], []
        for case in cases:
            count = len(case["pairs"])
            indexes = [rng.randrange(count) for _ in range(count)]
            reference = math.fsum(case["reference"][i] for i in indexes) / count
            native = math.fsum(case["native"][i] for i in indexes) / count
            reference_means.append(reference)
            native_means.append(native)
            case["bootstrap"].append(native / reference)
        suite_bootstrap.append(math.fsum(native_means) / math.fsum(reference_means))
    results = []
    for case in cases:
        results.append({"dataset": case["dataset"], "program": case["program"],
                        "measured_pairs": len(case["pairs"]),
                        **result_statistics(statistics.mean(case["reference"]), statistics.mean(case["native"]),
                                            case["bootstrap"])})
    suite = result_statistics(math.fsum(case["reference_mean_wall_seconds"] for case in results),
                              math.fsum(case["native_mean_wall_seconds"] for case in results), suite_bootstrap)
    suite["reference_total_mean_wall_seconds"] = suite.pop("reference_mean_wall_seconds")
    suite["native_total_mean_wall_seconds"] = suite.pop("native_mean_wall_seconds")
    suite["case_count"] = len(results)
    suite["case_classifications"] = {name: sum(case["classification"] == name for case in results)
                                      for name in ("win", "regression", "inconclusive")}
    return {"schema_version": 1, "status": "complete", "measurement_integrity_passed": True,
            "scope": "One fixed whole-workflow invocation per each of six cases; not engine-only or universal speed",
            "interpretation": "Applies only to these pinned inputs, artifact, settings and host. A suite win does not imply every case wins; all case results are reported without selection.",
            "method": {"suite_metric": "sum of case native arithmetic mean wall times / sum of case reference arithmetic mean wall times",
                       "weighting": "one complete fixed workload per dataset/preset; no mean of ratios or speedups",
                       "bootstrap": "percentile95%; independently resample complete paired repetition indexes within each case, then ratio of sums of case means; case intervals use the same joint draws",
                       "case_intervals": "pointwise95%, not simultaneous family-wise intervals; fixed-suite interval is primary",
                       "bootstrap_samples": BOOTSTRAP_SAMPLES, "bootstrap_seed": BOOTSTRAP_SEED,
                       "minimum_pairs_each_case": MINIMUM_PAIRS,
                       "classification": "win iff CIupper<1; regression iff CIlower>1; otherwise inconclusive",
                       "trial_exclusion": "only predeclared warmups excluded; all measured trials/outliers retained",
                       "workflow_difference": shared["workflow_difference"], "timing_boundary": shared["timing_boundary"],
                       "raw_observations": "all original reports, pairs and warmups embedded below without omissions"},
            "common_identity": shared,
            "pinned_workloads": {item["dataset"]: {"manifest": item["manifest_identity"], "inputs": item["pinned_inputs"]}
                                 for item in ordered_reports},
            "cases": results, "suite": suite,
            "source_reports": [{key: item[key] for key in ("path", "sha256", "record")} for item in ordered_reports]}


def print_table(summary):
    print("Dataset | Preset | Pairs | Reference s | Native s | Native/reference [95% CI] | Result")
    print("--- | --- | ---: | ---: | ---: | --- | ---")
    for case in summary["cases"]:
        low, high = case["paired_bootstrap_95pct_native_over_reference"]
        print(f"{case['dataset']} | {case['program']} | {case['measured_pairs']} | "
              f"{case['reference_mean_wall_seconds']:.6f} | {case['native_mean_wall_seconds']:.6f} | "
              f"{case['native_over_reference_ratio_of_means']:.4f} [{low:.4f}, {high:.4f}] | {case['classification']}")
    suite = summary["suite"]
    low, high = suite["paired_bootstrap_95pct_native_over_reference"]
    print(f"Suite total: reference {suite['reference_total_mean_wall_seconds']:.6f}s, "
          f"native {suite['native_total_mean_wall_seconds']:.6f}s; native/reference "
          f"{suite['native_over_reference_ratio_of_means']:.4f} [{low:.4f}, {high:.4f}]: {suite['classification']}.")
    print("Fixed supplied workflows only. Reference includes fresh v5 disk database/process overhead; "
          "native uses in-memory construction. No engine-only or universal speed claim.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reports", nargs=2, type=Path, help="Complete staph and yeast paired benchmark reports")
    parser.add_argument("--output", required=True, type=Path, help="JSON summary, embedding all original raw observations")
    args = parser.parse_args()
    summary = summarize(args.reports)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    print_table(summary)
    print(f"Full summary/raw observations: {args.output.resolve()}")


if __name__ == "__main__":
    main()
