#!/usr/bin/env python3
"""Supplement existing paired reports with a reference-search-phase view.

Compare each native COMPLETE workflow with its pair's recorded original SEARCH
phase, excluding original database construction. No new measurements are run.
This favors the reference: native still pays parsing/encoding/database creation,
while original uses the database just created by its omitted build phase. It is
not an independently measured prebuilt-database or engine-only benchmark.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import random
import statistics

from summarize_benchmarks import (
    BOOTSTRAP_SAMPLES, BOOTSTRAP_SEED, PROGRAMS, positive_time,
    result_statistics, summarize,
)


def supplementary(reports):
    # Reuse all primary suite/artifact/input/parity/integrity validation gates.
    # Do not alter the primary whole-workflow summary or discard slower cases.
    primary = summarize(reports)
    cases = []
    for source in primary["source_reports"]:
        report = source["record"]
        dataset = report["workload_manifest"]["content"]["dataset"]
        for program in PROGRAMS:
            pairs = report["searches"][program]["pairs"]
            reference = [positive_time(pair["reference"]["search"]["wall_seconds"],
                                       f"{dataset}/{program} reference search phase") for pair in pairs]
            native = [pair["native"]["wall_seconds"] for pair in pairs]
            cases.append(dict(dataset=dataset, program=program, pairs=pairs,
                              reference=reference, native=native, bootstrap=[]))
    rng = random.Random(BOOTSTRAP_SEED)
    suite_ratios = []
    for _ in range(BOOTSTRAP_SAMPLES):
        reference_means, native_means = [], []
        for case in cases:
            count = len(case["pairs"])
            indexes = [rng.randrange(count) for _ in range(count)]
            reference = math.fsum(case["reference"][i] for i in indexes) / count
            native = math.fsum(case["native"][i] for i in indexes) / count
            case["bootstrap"].append(native / reference)
            reference_means.append(reference)
            native_means.append(native)
        suite_ratios.append(math.fsum(native_means) / math.fsum(reference_means))

    def stats(reference, native, ratios):
        values = result_statistics(reference, native, ratios)
        values["reference_search_mean_wall_seconds"] = values.pop("reference_mean_wall_seconds")
        values["statistically_supported_lower_native_wall_time"] = values.pop("statistically_supported_workflow_improvement")
        return values

    results = [dict(dataset=case["dataset"], program=case["program"],
                    measured_pairs=len(case["pairs"]),
                    **stats(statistics.mean(case["reference"]), statistics.mean(case["native"]),
                            case["bootstrap"])) for case in cases]
    suite = stats(math.fsum(case["reference_search_mean_wall_seconds"] for case in results),
                  math.fsum(case["native_mean_wall_seconds"] for case in results), suite_ratios)
    suite["reference_search_total_mean_wall_seconds"] = suite.pop("reference_search_mean_wall_seconds")
    suite["native_total_mean_wall_seconds"] = suite.pop("native_mean_wall_seconds")
    suite["case_classifications"] = {name: sum(case["classification"] == name for case in results)
                                     for name in ("win", "regression", "inconclusive")}
    return dict(schema_version=1, status="complete", measurement_integrity_passed=True,
                scope=__doc__.strip(),
                method=dict(
                    numerator="native complete workflow wall_seconds, including in-memory DB creation/output",
                    denominator="reference.search.wall_seconds, excluding preceding original DB build phase",
                    suite_metric="ratio of sums of six case arithmetic means, not average ratios",
                    bootstrap="independently resample complete paired repetition indexes per case; joint percentile95% suite interval; pointwise case intervals",
                    bootstrap_samples=BOOTSTRAP_SAMPLES, bootstrap_seed=BOOTSTRAP_SEED,
                    trial_exclusion="predeclared warmups only; all measured pairs/outliers retained",
                    timing_scope="post-hoc view of the existing experiment, not new engine-only or prebuilt-DB measurement"),
                common_identity=primary["common_identity"],
                cases=results, suite=suite,
                original_primary_suite=primary["suite"],
                source_reports=primary["source_reports"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reports", nargs=2, type=Path, help="the complete pinned staph and yeast reports")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    summary = supplementary(args.reports)
    summary["supplementary_script_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(summary, indent=2) + "\n")
    print("Native COMPLETE workflow / recorded reference SEARCH phase (not engine-only)")
    print("Dataset | Preset | Reference search s | Native total s | Ratio [95% interval] | Result")
    for case in summary["cases"]:
        low, high = case["paired_bootstrap_95pct_native_over_reference"]
        print(f"{case['dataset']} | {case['program']} | {case['reference_search_mean_wall_seconds']:.6f} | "
              f"{case['native_mean_wall_seconds']:.6f} | {case['native_over_reference_ratio_of_means']:.5f} "
              f"[{low:.5f}, {high:.5f}] | {case['classification']}")
    suite = summary["suite"]
    low, high = suite["paired_bootstrap_95pct_native_over_reference"]
    print(f"Suite ratio {suite['native_over_reference_ratio_of_means']:.5f} [{low:.5f}, {high:.5f}]: "
          f"{suite['classification']}; all six cases retained")


if __name__ == "__main__":
    main()
