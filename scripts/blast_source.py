#!/usr/bin/env python3
"""Fetch the pinned reference, or verify its archive and selected reference files."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import tarfile
import urllib.request

REPO = Path(__file__).resolve().parents[1]


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("fetch", "verify"))
    args = parser.parse_args()
    pin = json.loads((REPO / "upstream/blast-source.json").read_text())
    cargo = (REPO / "Cargo.toml").read_text().split("[package.metadata.upstream]", 1)[1]
    metadata_version = re.search(r'^version = "([^"]+)"$', cargo, re.MULTILINE)
    library_version = re.search(
        r'UPSTREAM_BLAST_VERSION: &str = "([^"]+)"',
        (REPO / "src/lib.rs").read_text(),
    )
    if (not metadata_version or not library_version or
            metadata_version[1] != pin["version"] or library_version[1] != pin["version"]):
        raise SystemExit("Upstream version differs between source pin and Rust crate")
    archive = REPO / "vendor/downloads" / pin["archive"]
    if args.action == "fetch" and not archive.exists():
        archive.parent.mkdir(parents=True, exist_ok=True)
        partial = archive.with_suffix(archive.suffix + ".part")
        with urllib.request.urlopen(pin["url"], timeout=120) as source:
            with partial.open("wb") as destination:
                while block := source.read(1024 * 1024):
                    destination.write(block)
        if sha256(partial) != pin["sha256"]:
            raise SystemExit("Downloaded archive does not match the pinned SHA256")
        partial.replace(archive)
    if not archive.exists() or sha256(archive) != pin["sha256"]:
        raise SystemExit("Reference archive missing or SHA256 mismatch")
    root = REPO / pin["root"]
    if args.action == "fetch" and not root.exists():
        with tarfile.open(archive) as source:
            source.extractall(REPO / "vendor", filter="data")
    files = json.loads((REPO / "upstream/reference-files.json").read_text())
    for entry in files:
        path = root / entry["path"]
        if not path.exists() or sha256(path) != entry["sha256"]:
            raise SystemExit(f"Reference missing or modified: {entry['path']}")
    print(f"BLAST+ {pin['version']}: archive and {len(files)} reference files verified")


if __name__ == "__main__":
    main()
