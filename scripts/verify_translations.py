#!/usr/bin/env python3
"""Check versioned function correspondence and recorded C-fixture provenance."""
import hashlib
import json
from pathlib import Path
import re

REPO = Path(__file__).resolve().parents[1]


def digest(data):
    return hashlib.sha256(data).hexdigest()


def kernel_functions(path, include_crate=False):
    """Rustfmt-formatted free functions and simple inherent impl methods."""
    enclosing = ""
    for line in path.read_text().splitlines():
        implementation = re.match(r"^impl (\w+) \{", line)
        if implementation:
            enclosing = implementation[1]
        elif line == "}":
            enclosing = ""
        visibility = r"pub(?:\(crate\))?" if include_crate else r"pub"
        function = re.match(r"\s*" + visibility + r" (?:unsafe )?fn (\w+)(?:<.*>)?\s*\(", line)
        if function:
            yield enclosing, function[1]


def main():
    pin = json.loads((REPO / "upstream/blast-source.json").read_text())
    ledger = json.loads((REPO / "upstream/translations.json").read_text())
    fixture = json.loads((REPO / "upstream/kernel-fixture.json").read_text())
    originals, targets = set(), set()
    for entry in ledger:
        assert entry["version"] == pin["version"], entry
        original = (entry["source"], entry["symbol"])
        target = (entry["rust_path"], entry.get("rust_type", ""), entry["rust_symbol"])
        assert original not in originals and target not in targets, entry
        originals.add(original)
        targets.add(target)
        assert (entry.get("rust_type", ""), entry["rust_symbol"]) in set(kernel_functions(REPO / entry["rust_path"], include_crate=True)), target
        reference = REPO / pin["root"] / entry["source"]
        if reference.exists():
            original_bytes = reference.read_bytes()[entry["byte_start"]:entry["byte_end"]]
            assert digest(original_bytes) == entry["source_sha256"], original
    # Parallel workers may have new, unexposed source files in progress.
    # Audit every integrated core module and every ledger target; exposing a
    # new module makes its public functions mandatory in this same check.
    core_modules = set(re.findall(r"^pub mod (\w+);", (REPO / "src/core/mod.rs").read_text(), re.MULTILINE))
    mapped_paths = {entry["rust_path"] for entry in ledger}
    for directory in ("src/core", "src/seqdb", "src/util", "src/dustmask"):
        for path in (REPO / directory).glob("*.rs"):
            relative = path.relative_to(REPO).as_posix()
            if directory == "src/core" and path.stem not in core_modules and relative not in mapped_paths:
                continue
            for enclosing, name in kernel_functions(path):
                assert (path.relative_to(REPO).as_posix(), enclosing, name) in targets, f"Unmapped kernel: {path}:{enclosing}:{name}"
    data = json.loads((REPO / "upstream/translation-data.json").read_text())
    for entry in data:
        assert entry["version"] == pin["version"]
        text = (REPO / entry["rust_path"]).read_text()
        symbols = entry.get("symbols")
        if symbols is None:
            # Range records preserve exact provenance for individual tables,
            # declarations, or grouped constants from the composition closure.
            symbols = {entry["symbol"]: entry["rust_symbol"]}
        for targets in symbols.values():
            for target in targets.split("/"):
                name = target.rsplit("::", 1)[-1].split("<", 1)[0]
                assert re.search(r"\b(?:const|static|enum|struct|trait|type) " + re.escape(name) + r"\b", text), target
        reference = REPO / pin["root"] / entry["source"]
        if reference.exists():
            original = reference.read_bytes()
            if "byte_start" in entry:
                original = original[entry["byte_start"]:entry["byte_end"]]
            assert digest(original) == entry["source_sha256"], entry["source"]
    assert fixture["version"] == pin["version"]
    assert fixture["archive_sha256"] == pin["sha256"]
    for key in ("fixture", "harness", "generator", "shim"):
        assert digest((REPO / fixture[key]).read_bytes()) == fixture[key + "_sha256"], key
    records = len((REPO / fixture["fixture"]).read_text().splitlines())
    assert records == fixture["records"]
    search = json.loads((REPO / "tests/fixtures/search-smoke-2.17.0.json").read_text())
    assert search["reference_version"] == pin["version"]
    for key, filename in (("database", "subjects.fna"), ("nucleotide_queries", "nucleotide-queries.fna"),
                          ("protein_queries", "protein-queries.faa")):
        assert digest((REPO / "tests/data/smoke" / filename).read_bytes()) == search["inputs"][key]["sha256"]
    cpp_records = 0
    database_subjects = 0
    dust_records = 0
    query_records = 0
    lookup_records = 0
    ordinary_records = 0
    mb_records = 0
    seed_records = 0
    context_records = 0
    ungapped_records = 0
    extend_records = 0
    initial_hit_records = 0
    mini_extension_records = 0
    word_finder_records = 0
    statistics_records = 0
    composition_records = 0
    for filename in ("database-fixture.json", "makeblastdb-fixture.json", "symdust-fixture.json", "query-fixture.json", "lookup-fixture.json", "ordinary-lookup-fixture.json", "mb-lookup-fixture.json", "seed-lookup-fixture.json", "context-fixture.json", "ungapped-fixture.json", "extend-fixture.json", "initial-hit-fixture.json", "mini-extension-fixture.json", "word-finder-fixture.json", "statistics-fixture.json", "composition-fixture.json"):
        reference = json.loads((REPO / "upstream" / filename).read_text())
        assert reference["version"] == pin["version"]
        assert reference["archive_sha256"] == pin["sha256"]
        for key in ("fixture", "harness", "generator", "shim", "rng_shim", "range_parser", "integer_shim", "dust_shim", "dust_generator"):
            if key in reference:
                assert digest((REPO / reference[key]).read_bytes()) == reference[key + "_sha256"], key
        for entry in reference.get("fragments", []):
            path = REPO / pin["root"] / entry["source"]
            if path.exists():
                assert digest(path.read_bytes()[entry["byte_start"]:entry["byte_end"]]) == entry["sha256"]
        lines = (REPO / reference["fixture"]).read_text().splitlines()
        if filename in ("database-fixture.json", "symdust-fixture.json", "query-fixture.json", "lookup-fixture.json", "ordinary-lookup-fixture.json", "mb-lookup-fixture.json", "seed-lookup-fixture.json", "context-fixture.json", "ungapped-fixture.json", "extend-fixture.json", "initial-hit-fixture.json", "mini-extension-fixture.json", "word-finder-fixture.json", "statistics-fixture.json", "composition-fixture.json"):
            assert len(lines) == reference["records"]
            if filename == "database-fixture.json": cpp_records = len(lines)
            elif filename == "symdust-fixture.json": dust_records = len(lines)
            elif filename == "lookup-fixture.json": lookup_records = len(lines)
            elif filename == "ordinary-lookup-fixture.json": ordinary_records = len(lines)
            elif filename == "mb-lookup-fixture.json": mb_records = len(lines)
            elif filename == "seed-lookup-fixture.json": seed_records = len(lines)
            elif filename == "context-fixture.json": context_records = len(lines)
            elif filename == "ungapped-fixture.json": ungapped_records = len(lines)
            elif filename == "extend-fixture.json": extend_records = len(lines)
            elif filename == "initial-hit-fixture.json": initial_hit_records = len(lines)
            elif filename == "mini-extension-fixture.json": mini_extension_records = len(lines)
            elif filename == "word-finder-fixture.json": word_finder_records = len(lines)
            elif filename == "statistics-fixture.json": statistics_records = len(lines)
            elif filename == "composition-fixture.json": composition_records = len(lines)
            else:
                query_records = len(lines)
                tools = json.loads((REPO / "upstream/blast-tools.json").read_text())
                assert reference["executable_archive_sha256"] == tools["sha256"]
        else:
            tools = json.loads((REPO / "upstream/blast-tools.json").read_text())
            assert tools["version"] == pin["version"]
            assert reference["executable_archive_sha256"] == tools["sha256"]
            assert lines[0].split()[0] == "stats"
            assert len(lines) - 1 == int(lines[0].split()[1]) == reference["records"]
            database_subjects = reference["records"]
    print(f"{len(ledger)} mappings; {records} C, {cpp_records} writer/RNG, {dust_records} DUST, {lookup_records} lookup, {ordinary_records} ordinary, {mb_records} megablast, {seed_records} seed-adapter, {context_records} context and {ungapped_records} extension and {extend_records} extension-state and {initial_hit_records} initial-hit and {mini_extension_records} mini-extension and {word_finder_records} word-finder and {statistics_records} statistics and {composition_records} composition records; {database_subjects} database subjects and {query_records} original query setups verified")


if __name__ == "__main__":
    main()
