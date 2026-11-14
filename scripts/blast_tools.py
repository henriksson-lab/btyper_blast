#!/usr/bin/env python3
"""Fetch or verify the pinned original BLAST Linux x86_64 reference executables."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import subprocess
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
    if platform.system() != "Linux" or platform.machine().lower() not in {"x86_64", "amd64"}:
        raise SystemExit("This reference binary pin is for Linux x86_64 only; Rust library targets are separate")
    pin = json.loads((REPO / "upstream/blast-tools.json").read_text())
    source = json.loads((REPO / "upstream/blast-source.json").read_text())
    if pin["version"] != source["version"]:
        raise SystemExit("Executable/source reference versions differ")
    archive = REPO / "vendor/downloads" / pin["archive"]
    if args.action == "fetch" and not archive.exists():
        archive.parent.mkdir(parents=True, exist_ok=True)
        partial = archive.with_suffix(archive.suffix + ".part")
        with urllib.request.urlopen(pin["url"], timeout=300) as response, partial.open("wb") as output:
            while block := response.read(1024 * 1024):
                output.write(block)
        if sha256(partial) != pin["sha256"]:
            raise SystemExit("Downloaded reference executable archive SHA256 mismatch")
        partial.replace(archive)
    if not archive.exists() or sha256(archive) != pin["sha256"]:
        raise SystemExit("Reference executable archive missing or SHA256 mismatch")
    root = REPO / pin["root"]
    if args.action == "fetch":
        wanted = {f"ncbi-blast-{pin['version']}/bin/{name}" for name in pin["tools"]}
        with tarfile.open(archive) as reference:
            members = [member for member in reference.getmembers() if member.name in wanted]
            if {member.name for member in members} != wanted:
                raise SystemExit("Reference archive does not contain all required tools")
            reference.extractall(REPO / "vendor/reference-tools", members=members, filter="data")
    with tarfile.open(archive) as reference:
        for name in pin["tools"]:
            path = root / name
            original = reference.extractfile(f"ncbi-blast-{pin['version']}/bin/{name}")
            if original is None or not path.exists():
                raise SystemExit(f"Missing original reference executable: {name}")
            expected = hashlib.sha256()
            while block := original.read(1024 * 1024):
                expected.update(block)
            if sha256(path) != expected.hexdigest():
                raise SystemExit(f"Original executable modified: {name}")
            result = subprocess.run([str(path), "-version"], capture_output=True, text=True, check=True)
            first_line = result.stdout.splitlines()[0]
            version = pin.get("tool_versions", {}).get(name, pin["version"])
            if first_line != f"{name}: {version}":
                raise SystemExit(f"Wrong reference version: {first_line}")
            print(first_line)


if __name__ == "__main__":
    main()
