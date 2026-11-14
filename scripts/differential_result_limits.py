#!/usr/bin/env python3
"""Focused whole-pipeline default subject-limit oracle; not a benchmark.

Compact deterministic workloads exercise equal-score subject ties,
interleaved score grades, and coupled CBS heaps at the default limits.
Compare all 19 tabular fields in order against the pinned original executable.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / 'scripts'))
from compare_search import compare
from reference_search import FIELDS
from differential_protein import AMINO, encode


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_fasta(path, records):
    path.write_text(''.join(f'>{name}\n{sequence}\n' for name, sequence in records))


def workload(program, graded):
    rng = random.Random(217550 if program != 'tblastn' else 2171050)
    alphabet = AMINO if program == 'tblastn' else 'ACGT'
    length, count = (40, 1101) if program == 'tblastn' else (100, 551)
    query = ''.join(rng.choice(alphabet) for _ in range(length))
    subjects = []
    grades = (0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5)
    for index in range(count):
        changes = grades[(index * 73) % len(grades)] if graded else 0
        sequence = list(query)
        for position in (8, 13, 20, 27, 32)[:changes]:
            sequence[position] = alphabet[(alphabet.index(sequence[position]) + 1) % len(alphabet)]
        sequence = ''.join(sequence)
        subjects.append((f's{index:04d}_m{changes}', encode(sequence) if program == 'tblastn' else sequence))
    return [('query', query)], subjects


def weak_multi_heap_workload():
    # The original early-termination gate requires EVERY query heap to hold
    # 500 matches with worst E <= .002. Weak positives alone cannot fill that
    # condition. After the shared strong prefix, query A's weak tail must run
    # while query B's heap is empty; query B's weak tail can be skipped later.
    rng = random.Random(2170500)
    query = ''.join(rng.choice(AMINO) for _ in range(40))
    # Exact seven-residue BLOSUM62 segment scores42. Stop flanks prevent
    # incidental extension; full target length remains40aa for both classes.
    weak = 'TAA' * 12 + encode(query[1:8]) + 'TAA' * 21
    subjects = [(f's{i:04d}_{"strong" if i < 500 else "weak"}',
                 encode(query) if i < 500 else weak) for i in range(1101)]
    return [('qA', query), ('qB', query)], subjects


def run(command, directory, env):
    result = subprocess.run(list(map(str, command)), cwd=directory, env=env,
                            capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError(f'{command[0]} exited {result.returncode}: {result.stderr}')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--native', required=True, type=Path)
    parser.add_argument('--program', action='append', choices=('megablast', 'blastn', 'tblastn'))
    parser.add_argument('--case', action='append', choices=('equal', 'rank_boundary', 'weak_multi_heap'))
    parser.add_argument('--work-dir', type=Path)
    args = parser.parse_args()
    subprocess.run([sys.executable, str(REPO / 'scripts/blast_tools.py'), 'verify'],
                   check=True, stdout=subprocess.DEVNULL)
    pin = json.loads((REPO / 'upstream/blast-tools.json').read_text())
    source = json.loads((REPO / 'upstream/blast-source.json').read_text())
    if pin['version'] != source['version']:
        raise ValueError('Source/reference version pins differ')
    reference = REPO / pin['root']
    native = args.native.resolve()
    env = dict(os.environ, BLAST_USAGE_REPORT='false', NCBI_DONT_USE_NCBIRC='1',
               NCBI_DONT_USE_LOCAL_CONFIG='1')
    for name in ('CHUNK_SIZE', 'OVERLAP_CHUNK_SIZE', 'ADAPTIVE_CBS'):
        env.pop(name, None)
    work = args.work_dir.resolve() if args.work_dir else Path(tempfile.mkdtemp(prefix='btyper-result-limits-'))
    work.mkdir(parents=True, exist_ok=True)
    report = {'version': pin['version'], 'source_pin': source, 'reference_archive_sha256': pin['sha256'],
              'native': str(native), 'native_sha256': digest(native), 'work_dir': str(work),
              'scope': 'All 19 CLI fields and order; strict 500-subject default-cap assertions. Not raw float-bit parity or speed evidence.',
              'runs': []}
    for program in args.program or ('megablast', 'blastn', 'tblastn'):
        tool = reference / ('tblastn' if program == 'tblastn' else 'blastn')
        version = run([tool, '-version'], REPO, env)
        if pin['version'] not in version.stdout:
            raise ValueError('Executable version differs from pin')
        cases = ('equal', 'rank_boundary', 'weak_multi_heap') if program == 'tblastn' else ('equal', 'rank_boundary')
        for name in cases:
            if args.case and name not in args.case:
                continue
            directory = work / f'{program}-{name}'
            directory.mkdir(exist_ok=True)
            queries, subjects = (weak_multi_heap_workload() if name == 'weak_multi_heap'
                                 else workload(program, name == 'rank_boundary'))
            query_path = directory / ('queries.faa' if program == 'tblastn' else 'queries.fna')
            db_path = directory / 'subjects.fna'
            write_fasta(query_path, queries)
            write_fasta(db_path, subjects)
            run([reference / 'makeblastdb', '-in', db_path, '-dbtype', 'nucl',
                 '-blastdb_version', '5', '-out', directory / 'db'], directory, env)
            original_command = [tool, '-query', query_path, '-db', directory / 'db',
                                '-num_threads', '1', '-outfmt', '6 ' + ' '.join(FIELDS)]
            if program != 'tblastn':
                original_command += ['-task', program]
            native_command = [native, program, db_path, query_path]
            original = run(original_command, directory, env)
            actual = run(native_command, directory, env)
            for label, output in (('original', original), ('native', actual)):
                (directory / f'{label}.tsv').write_text(output.stdout)
                (directory / f'{label}.stderr').write_text(output.stderr)
            entry = {'program': program, 'case': name, 'input_subjects': len(subjects),
                     'input_queries': len(queries),
                     'preliminary_limit': 1050 if program == 'tblastn' else 550,
                     'final_limit': 500, 'query_sha256': digest(query_path),
                     'subjects_sha256': digest(db_path), 'original_executable_sha256': digest(tool),
                     'original_command': list(map(str, original_command)), 'native_command': list(map(str, native_command))}
            try:
                expected = [dict(zip(FIELDS, row.split('\t'))) for row in original.stdout.splitlines()]
                rows = [row.split('\t') for row in actual.stdout.splitlines()]
                compare(rows, expected, program)
                original_ids = {row['sseqid'] for row in expected}
                actual_ids = {row[1] for row in rows}
                if len(original_ids) != 500 or actual_ids != original_ids:
                    raise ValueError(f'default-cap exercise incomplete: original {len(original_ids)}, native {len(actual_ids)} subjects')
                subjects_by_query = {query: {row['sseqid'] for row in expected if row['qseqid'] == query}
                                     for query, _ in queries}
                if any(len(values) != 500 for values in subjects_by_query.values()):
                    raise ValueError('Each query must independently exercise the final500-subject cap')
                entry.update(status='match', hsps=len(rows), distinct_subjects=500,
                             subjects_by_query={query: len(values) for query, values in subjects_by_query.items()},
                             original_output_sha256=digest(directory / 'original.tsv'),
                             native_output_sha256=digest(directory / 'native.tsv'))
            except Exception as error:
                entry.update(status='mismatch', error=str(error))
                report['runs'].append(entry)
                (work / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
                raise
            report['runs'].append(entry)
            (work / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
            print(f'{program} {name}: {len(rows)} HSPs / 500 subjects match all 19 fields and order', flush=True)
    if not report['runs']:
        raise ValueError('No cases selected for the requested programs')
    if digest(native) != report['native_sha256']:
        raise ValueError('Native executable changed during validation')
    report['status'] = 'match'
    (work / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(work / 'report.json')


if __name__ == '__main__':
    main()
